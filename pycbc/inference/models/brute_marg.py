 # Copyright (C) 2020 Alex Nitz
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

"""This module provides model classes that do brute force marginalization
using at the likelihood level.
"""
import math
import logging
import numpy

from pycbc.pool import BroadcastPool as Pool
from scipy.special import logsumexp

from .gaussian_noise import BaseGaussianNoise
from .tools import draw_sample, DistMarg

from pycbc.conversions import mass1_from_mtotal_q, mass2_from_mtotal_q
from pycbc.distributions import JointDistribution

_model = None
class likelihood_wrapper(object):
    def __init__(self, model):
        global _model
        _model = model

    def __call__(self, params):
        global _model

        _model.update(**params)

        loglr = _model.loglr
        lognl = _model.lognl
        loglikelihood = lognl + loglr
        logjacobian = _model.logjacobian
        logprior = _model.logprior

        return loglr, _model.current_stats


class BruteTotalMassMarginalize(BaseGaussianNoise, DistMarg):
    name = "marginalized_mtotal"

    def __init__(self, variable_params,
                 cores=1,
                 base_model=None,
                 marginalize_mtotal=None,
                 mtotal_grid=None,
                 mtotal_grid_num=None,
                 fiducial_mtotal=None,
                 **kwds):
        from pycbc.inference.models import models

        self.marginalize_vector_params = []
        if 'marginalize_vector_params' in kwds:
            if isinstance(kwds['marginalize_vector_params'], str) and kwds['marginalize_vector_params']:
                self.marginalize_vector_params = [p.strip() for p in kwds['marginalize_vector_params'].split(',') if p.strip()]
            elif isinstance(kwds['marginalize_vector_params'], (list, tuple)):
                self.marginalize_vector_params = [p.strip() for p in kwds['marginalize_vector_params'] if p.strip()]

        # --- Save the original mtotal prior distribution and bounds before deletion ---
        self.mtotal_prior_dist = None
        self.mtotal_bounds = None
        prior = kwds.get('prior', None)
        if prior and 'mtotal' in prior.bounds:
            for dist in prior.distributions:
                if dist.name == 'mtotal':
                    self.mtotal_prior_dist = dist
                    break
            self.mtotal_bounds = (prior.bounds['mtotal'].min, prior.bounds['mtotal'].max)

        # --- Setup bounds for mtotal marginalization using the saved bounds ---
        if marginalize_mtotal:
            prior = kwds.get('prior', {})
            if prior:
                prior_bounds = prior.bounds
                min_mtotal = prior_bounds['mtotal'].min
                max_mtotal = prior_bounds['mtotal'].max
            else:
                min_mtotal = max_mtotal = None

            if fiducial_mtotal is None and min_mtotal is not None and max_mtotal is not None:
                fiducial_mtotal = 0.5 * (min_mtotal + max_mtotal)
            else:
                fiducial_mtotal = float(fiducial_mtotal) if fiducial_mtotal else None

            if mtotal_grid is None and min_mtotal is not None and max_mtotal is not None:
                if mtotal_grid_num is None:
                    raise ValueError("Must specify mtotal_grid_num if mtotal_grid is not given")
                self.mtotal_grid = numpy.linspace(min_mtotal, max_mtotal, int(mtotal_grid_num))
            elif mtotal_grid is not None:
                self.mtotal_grid = numpy.array([float(m) for m in mtotal_grid])
            else:
                raise ValueError("Either provide 'mtotal_grid' (a list or array) or specify 'mtotal_grid_num' (the number of points to create a grid over the prior).")

        base_model_cls = models[base_model]
        self.model = base_model_cls(variable_params=variable_params, **kwds)
        self.call = likelihood_wrapper(self.model)

        marginalized_params = []

        if marginalize_mtotal:
            marginalized_params.append('mtotal')
        if 'marginalize_phase' in kwds and kwds['marginalize_phase']:
            marginalized_params.append('coa_phase')
        if 'marginalize_distance' in kwds and kwds['marginalize_distance']:
            marginalized_params.append(kwds.get('marginalize_distance_param', 'distance'))

        # ------ Reconstruct the marginalized prior ------
        marginalized_params.extend(self.marginalize_vector_params)

        variable_params = tuple(p for p in variable_params if p not in marginalized_params)

        prior_dists_by_param = {}
        for d in self.model.prior_distribution.distributions:
            if hasattr(d, 'params') and isinstance(d.params, (list, tuple)) and len(d.params) > 0:
                for p_name in d.params:
                    if isinstance(p_name, str):
                        prior_dists_by_param[p_name] = d
            else:
                if isinstance(d.name, str) and d.name in variable_params:
                    prior_dists_by_param[d.name] = d

        current_prior_dists = prior_dists_by_param.copy()
        for param in marginalized_params:
            if param in current_prior_dists:
                del current_prior_dists[param]

        variable_params = tuple(p for p in variable_params if p not in marginalized_params)

        self.model.prior_distribution = JointDistribution(variable_params, *current_prior_dists.values())
        kwds['prior'] = self.model.prior_distribution

        super().__init__(variable_params, **kwds)

        # Set up multiprocessing
        if cores > 1:
            self.pool = Pool(int(cores))
            self.mapfunc = self.pool.map
        else:
            self.pool = None
            self.mapfunc = map

        self.marginalize_mtotal = marginalize_mtotal
        self.fiducial_mtotal = fiducial_mtotal

    @property
    def _extra_stats(self):
        stats = self.model._extra_stats
        stats.append('maxl_mtotal')
        if 'maxl_loglr' not in stats:
            stats.append('maxl_loglr')
        return stats

    def scale_waveform(self, h_plus_ref, h_cross_ref, mtotal_ref, mtotal_new, q, f_lower, approximant):
        # Compute scaling factors
        time_scale = mtotal_new / mtotal_ref
        amp_scale = mtotal_new / mtotal_ref

        # Rescale amplitude
        h_plus_scaled = h_plus_ref * amp_scale
        h_cross_scaled = h_cross_ref * amp_scale

        return h_plus_scaled, h_cross_scaled

    def _loglr(self):
        if self.mtotal_grid is not None:
            ref_params = self.current_params.copy()
            for key in ['ra', 'dec', 'tc', 'polarization', 'distance','mtotal']:
                ref_params.pop(key, None)
            ref_params["mtotal"] = self.fiducial_mtotal
            ref_params["mass1"] = mass1_from_mtotal_q(ref_params["mtotal"], ref_params["q"])
            ref_params["mass2"] = mass2_from_mtotal_q(ref_params["mtotal"], ref_params["q"])

            wfs = self.model.waveform_generator.generate(**ref_params)

            self.reference_waveform = wfs

            approximant = ref_params['approximant']
            f_lower = ref_params['f_lower']
            q = self.current_params["q"]

            params = []
            self.scaled_waveforms = {}

            params = []

            for mtotal in self.mtotal_grid:
                scaled_waveform = {}
                for det, (hp, hc) in wfs.items():
                    scaled_hp, scaled_hc = self.scale_waveform(
                        hp, hc, self.fiducial_mtotal, mtotal, q, f_lower,
                        approximant=approximant
                    )
                    scaled_waveform[det] = (scaled_hp, scaled_hc)

                self.scaled_waveforms[mtotal] = scaled_waveform

                pmod = ref_params.copy()
                pmod.pop("mass1", None)
                pmod.pop("mass2", None)
                pmod.pop("mtotal", None)
                pmod["custom_waveform"] = self.scaled_waveforms[mtotal]
                pmod["rescale_mtotal"] = mtotal
                params.append(pmod)

            # Check first loglr value to decide whether to proceed scaling the waveform (skips calculation for loglr <0)
            first_val = self.call(params[0])
            first_loglr = first_val[0]

            if first_loglr < 0:
                print("Early exit: First loglr < 0, skipping full grid evaluation.")
                self._current_stats.maxl_mtotal = self.mtotal_grid[0]
                self._current_stats.maxl_loglr = first_loglr
                for stat in first_val[1]:
                    setattr(self._current_stats, stat, first_val[1][stat])
                return first_loglr

            # Proceed with full marginalization if the first is acceptable
            params_full = params
            vals = list(self.mapfunc(self.call, params_full))

            loglr = numpy.array([v[0] for v in vals])
            maxidx = loglr.argmax()
            maxstats = vals[maxidx][1]
            max_mtotal = self.mtotal_grid[maxidx]

            for stat in maxstats:
                setattr(self._current_stats, stat, maxstats[stat])
            self._current_stats.maxl_mtotal = max_mtotal
            self._current_stats.maxl_loglr = loglr[maxidx]
            print(logsumexp(loglr) - numpy.log(len(self.mtotal_grid)))
            return logsumexp(loglr) - numpy.log(len(self.mtotal_grid))
        else:
            return self.model.loglr


class BruteParallelGaussianMarginalize(BaseGaussianNoise):
    name = "brute_parallel_gaussian_marginalize"

    def __init__(self, variable_params,
                 cores=10,
                 base_model=None,
                 marginalize_phase=None,
                 **kwds):
        super().__init__(variable_params, **kwds)

        from pycbc.inference.models import models
        self.model = models[base_model](variable_params, **kwds)

        self.call = likelihood_wrapper(self.model)

        # size of pool for each likelihood call
        self.pool = Pool(int(cores))

        # Only one for now, but can be easily extended
        self.phase = None
        if marginalize_phase:
            samples = int(marginalize_phase)
            self.phase = numpy.linspace(0, 2.0 * numpy.pi, samples)

    @property
    def _extra_stats(self):
        stats = self.model._extra_stats
        stats.append('maxl_phase')
        if 'maxl_loglr' not in stats:
            stats.append('maxl_loglr')
        return stats

    def _loglr(self):
        if self.phase is not None:
            params = []
            for p in self.phase:
                pref = self.current_params.copy()
                pref['coa_phase'] = p
                params.append(pref)
            vals = list(self.pool.map(self.call, params))
            loglr = numpy.array([v[0] for v in vals])
            # get the maxl values
            if 'maxl_loglr' not in self.model._extra_stats:
                maxl_loglrs = loglr
            else:
                maxl_loglrs = numpy.array([v[1]['maxl_loglr'] for v in vals])
            maxidx = maxl_loglrs.argmax()
            maxstats = vals[maxidx][1]
            maxphase = self.phase[maxidx]
            # set the stats
            for stat in maxstats:
                setattr(self._current_stats, stat, maxstats[stat])
            self._current_stats.maxl_phase = maxphase
            self._current_stats.maxl_loglr = maxl_loglrs[maxidx]
            # calculate the marginal loglr and return
            return logsumexp(loglr) - numpy.log(len(self.phase))


class BruteLISASkyModesMarginalize(BaseGaussianNoise):
    name = "brute_lisa_sky_modes_marginalize"

    def __init__(self, variable_params,
                 cores=1,
                 loop_polarization=False,
                 base_model=None,
                 **kwds):
        super().__init__(variable_params, **kwds)

        from pycbc.inference.models import models
        kwds.update(models[base_model].extra_args_from_config(
            kwds['config_object'],
            "model",
            skip_args=[])
        )

        self.model = models[base_model](variable_params, **kwds)

        self.call = likelihood_wrapper(self.model)

        # size of pool for each likelihood call
        if cores > 1:
            self.pool = Pool(int(cores))
            self.mapfunc = self.pool.map
        else:
            self.pool = None
            self.mapfunc = map

        # Do I explicitly check the polarization + pi/2 points
        # We could also add other arguments here, ie only check longitude
        # or latitude symmetry points.
        if loop_polarization:
            self.num_sky_modes = 16
        else:
            self.num_sky_modes = 8

        self.reconstruct_sky_points = False

    @property
    def _extra_stats(self):
        stats = self.model._extra_stats
        return stats

    def _loglr(self):
        params = []
        for sym_num in range(self.num_sky_modes):
            pref = self.current_params.copy()
            self._apply_sky_point_rotation(pref, sym_num)
            params.append(pref)

        vals = list(self.mapfunc(self.call, params))
        loglr = numpy.array([v[0] for v in vals])

        if self.reconstruct_sky_points:
            return loglr

        max_llr_idx = loglr.argmax()
        max_llr = loglr[max_llr_idx]
        marg_lrfac = sum([math.exp(llr - max_llr) for llr in loglr])
        marg_llr = max_llr + math.log(marg_lrfac/self.num_sky_modes)

        # set the stats
        for sym_num in range(self.num_sky_modes):
            setattr(self._current_stats, f'llr_mode_{sym_num}', loglr[sym_num])

        return marg_llr

    def _apply_sky_point_rotation(self, pref, sky_num):
        """ Apply the sky point rotation for mode sky_num to parameters pref
        """
        lambdal = pref['eclipticlongitude']
        beta = pref['eclipticlatitude']
        psi = pref['polarization']
        inc = pref['inclination']

        pol_num = sky_num // 8
        sky_num = sky_num % 8
        long_num = sky_num % 4
        lat_num = sky_num // 4

        # Apply latitude symmetry mode
        if lat_num:
            beta = - beta
            inc = numpy.pi - inc
            psi = numpy.pi - psi

        # Apply longitudonal symmetry mode
        lambdal = (lambdal + long_num * 0.5 * numpy.pi) % (2*numpy.pi)
        psi = (psi + long_num * 0.5 * numpy.pi) % (2*numpy.pi)

        # Apply additional polarization mode (shouldn't be needed)
        if pol_num:
            psi = psi + (math.pi / 2.)

        pref['eclipticlongitude'] = lambdal
        pref['eclipticlatitude'] = beta
        pref['polarization'] = psi
        pref['inclination'] = inc

    @classmethod
    def from_config(cls, cp, **kwargs):
        kwargs['config_object'] = cp
        return super(BruteLISASkyModesMarginalize, cls).from_config(
            cp,
            **kwargs
        )

    def reconstruct(self, seed=None):
        """ Reconstruct a point from unwrapping the 8-fold sky symmetry
        """
        if seed:
            numpy.random.seed(seed)
        rec = {}

        logging.info('Reconstruct LISA sky mode symmetry')
        self.reconstruct_sky_points = True
        loglr = self.loglr
        xl = draw_sample(loglr)
        logging.info('Found point %d', xl)
        # Undo rotations
        pref = self.current_params.copy()
        self._apply_sky_point_rotation(pref, xl)

        for val in ['polarization', 'eclipticlongitude', 'eclipticlatitude',
                    'inclination']:
            rec[val] = pref[val]
        rec['loglr'] = loglr[xl]
        rec['loglikelihood'] = self.lognl + rec['loglr']
        self.reconstruct_sky_points = False
        return self.model.reconstruct(seed=seed, rec=rec)
