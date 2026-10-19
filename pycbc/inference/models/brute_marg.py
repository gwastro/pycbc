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

"""
This module provides model classes that do brute force marginalization
using at the likelihood level.
"""

import logging
import math

import numpy
from scipy.special import logsumexp

from pycbc.conversions import mass1_from_mtotal_q, mass2_from_mtotal_q
from pycbc.distributions import JointDistribution
from pycbc.pool import BroadcastPool as Pool

from .gaussian_noise import BaseGaussianNoise
from .tools import DistMarg, draw_sample

_model = None


class likelihood_wrapper:
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
    r"""Brute force marginalization over the total mass.

    This model wraps a base model and marginalizes the likelihood over a grid
    of total masses. Rather than generating a new waveform at every point on
    that grid, a single reference waveform is generated at
    ``fiducial_mtotal`` and rescaled to each total mass.

    At fixed mass ratio, spins, inclination and distance, a compact binary
    waveform depends on the total mass only through the dimensionless
    combination :math:`t / M`, with an overall amplitude proportional to
    :math:`M`. In the time domain this is

    .. math:: h(t; M) = \kappa\, h(t / \kappa; M_{\mathrm{ref}}),

    with :math:`\kappa = M / M_{\mathrm{ref}}`, and therefore in the frequency
    domain

    .. math:: \tilde{h}(f; M) = \kappa^{2}\,
              \tilde{h}(\kappa f; M_{\mathrm{ref}}).

    The reference waveform is generated from a correspondingly lower starting
    frequency, :math:`\min(\kappa) f_{\mathrm{low}}`, so that the lightest
    point on the grid still spans the full analysis band once rescaled.

    Because the rescaling is applied to the generated waveform, the base model
    must accept a pre-generated waveform through the ``custom_waveform``
    parameter; at present only :py:class:`MarginalizedTime
    <pycbc.inference.models.marginalized_gaussian_noise.MarginalizedTime>`
    does. The base model is also responsible for applying the detector
    response, so ``ra``, ``dec``, ``tc``, ``polarization`` and ``distance``
    must either be static or marginalized over by the base model.

    The marginalization assumes a uniform prior on the total mass evaluated
    on a uniformly spaced grid. The grid has to be fine enough to resolve the
    peak of the likelihood, which narrows as the signal gets louder; too
    coarse a grid biases the marginalized likelihood. Check convergence by
    increasing ``mtotal_grid_num`` until the likelihood stops moving.

    Parameters
    ----------
    variable_params : (tuple of) string(s)
        The model parameters.
    cores : int, optional
        Number of cores to use to evaluate the total mass grid in parallel.
        Default is 1.
    base_model : str
        Name of the model to marginalize on top of.
    marginalize_mtotal : bool, optional
        Turn on the total mass marginalization. If not set, the base model's
        likelihood is returned unchanged.
    mtotal_grid : list of float, optional
        Explicit grid of total masses to marginalize over. If not given, the
        grid is built from the ``mtotal`` prior bounds and
        ``mtotal_grid_num``.
    mtotal_grid_num : int, optional
        Number of uniformly spaced grid points to place across the ``mtotal``
        prior bounds. Required if ``mtotal_grid`` is not given.
    fiducial_mtotal : float, optional
        Total mass of the reference waveform. Defaults to the midpoint of the
        ``mtotal`` prior bounds.

    Notes
    -----
    The rescaling is exact only while the rescaled waveform still fits inside
    the analysis segment. For long signals whose duration approaches the
    segment length the rescaled waveform is mildly inaccurate (a match of
    roughly 0.99 rather than 1.0); use a longer segment in that regime.
    """

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
        vector_params = kwds.get('marginalize_vector_params', None)
        if isinstance(vector_params, str):
            vector_params = vector_params.split(',')
        if vector_params:
            self.marginalize_vector_params = [
                p.strip() for p in vector_params if p.strip()]

        self.mtotal_grid = None

        # --- Setup bounds for mtotal marginalization from the mtotal prior ---
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
                    raise ValueError("Must specify mtotal_grid_num if "
                                     "mtotal_grid is not given")
                self.mtotal_grid = numpy.linspace(
                    min_mtotal, max_mtotal, int(mtotal_grid_num))
            elif mtotal_grid is not None:
                self.mtotal_grid = numpy.array([float(m) for m in mtotal_grid])
            else:
                raise ValueError(
                    "Either provide 'mtotal_grid' (a list or array) or "
                    "specify 'mtotal_grid_num' (the number of points to "
                    "create a grid over the prior).")

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

        variable_params = tuple(
            p for p in variable_params if p not in marginalized_params)

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

        variable_params = tuple(
            p for p in variable_params if p not in marginalized_params)

        self.model.prior_distribution = JointDistribution(
            variable_params, *current_prior_dists.values())
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

    @staticmethod
    def _rescale_frequencyseries(htilde, kappa):
        r"""Evaluate :math:`\kappa^2 \tilde{h}(\kappa f)` on the grid of
        ``htilde``.

        The amplitude and phase are interpolated separately. The phase is
        unwrapped by smoothing its increments rather than the phase itself:
        for a chirp the phase can advance by more than :math:`\pi` between
        neighbouring frequency bins, which aliases a direct unwrap, but the
        *change* in that advance stays far below :math:`\pi`.
        """
        data = htilde.numpy()
        out = htilde.copy()
        out.data[:] = 0

        nonzero = numpy.flatnonzero(data)
        if len(nonzero) < 2:
            return out
        kmin, kmax = nonzero[0], nonzero[-1]

        f = numpy.arange(kmin, kmax + 1) * htilde.delta_f
        band = data[kmin:kmax + 1]
        amp = numpy.abs(band)

        phase = numpy.angle(band)
        increments = numpy.diff(phase)
        increments = (increments + numpy.pi) % (2 * numpy.pi) - numpy.pi
        increments = numpy.unwrap(increments)
        phase = phase[0] + numpy.concatenate(([0.0], numpy.cumsum(increments)))

        target = kappa * numpy.arange(len(out)) * htilde.delta_f
        inband = (target >= f[0]) & (target <= f[-1])
        out.data[inband] = kappa ** 2 * (
            numpy.interp(target[inband], f, amp)
            * numpy.exp(1j * numpy.interp(target[inband], f, phase))
        )
        return out

    def scale_waveform(self, h_plus_ref, h_cross_ref, mtotal_ref, mtotal_new):
        """Rescale a reference waveform from ``mtotal_ref`` to ``mtotal_new``.

        See the class documentation for the scaling relation used.
        """
        kappa = mtotal_new / mtotal_ref
        return (self._rescale_frequencyseries(h_plus_ref, kappa),
                self._rescale_frequencyseries(h_cross_ref, kappa))

    def _loglr(self):
        if self.mtotal_grid is None:
            self.model.update(**self.current_params)
            return self.model.loglr

        ref_params = self.current_params.copy()
        # These are applied by the base model, not by the waveform generator.
        for key in ['ra', 'dec', 'tc', 'polarization', 'distance', 'mtotal']:
            ref_params.pop(key, None)
        ref_params["mtotal"] = self.fiducial_mtotal
        ref_params["mass1"] = mass1_from_mtotal_q(ref_params["mtotal"],
                                                  ref_params["q"])
        ref_params["mass2"] = mass2_from_mtotal_q(ref_params["mtotal"],
                                                  ref_params["q"])

        # Rescaling maps f -> kappa * f, so the reference waveform has to
        # start low enough that the lightest grid point still covers the
        # analysis band once rescaled.
        kappa_min = self.mtotal_grid.min() / self.fiducial_mtotal
        ref_params["f_lower"] = ref_params["f_lower"] * min(1.0, kappa_min)

        wfs = self.model.waveform_generator.generate(**ref_params)
        self.reference_waveform = wfs

        params = []
        self.scaled_waveforms = {}
        for mtotal in self.mtotal_grid:
            scaled_waveform = {}
            for det, (hp, hc) in wfs.items():
                scaled_waveform[det] = self.scale_waveform(
                    hp, hc, self.fiducial_mtotal, mtotal)
            self.scaled_waveforms[mtotal] = scaled_waveform

            pmod = ref_params.copy()
            for key in ("mass1", "mass2", "mtotal"):
                pmod.pop(key, None)
            pmod["custom_waveform"] = scaled_waveform
            params.append(pmod)

        vals = list(self.mapfunc(self.call, params))
        loglr = numpy.array([v[0] for v in vals])

        # kept so that reconstruct() can draw a total mass from the grid
        self._grid_params = params
        self._grid_loglr = loglr

        # Record the base model's own diagnostics (including any parameters
        # it reconstructs while marginalizing) from the best fitting point on
        # the grid. The likelihood terms themselves are deliberately left
        # alone: the base model is only ever asked for its loglr, so its
        # remaining stats are unset, and copying them over would mask the
        # values this model computes for itself.
        reserved = ('loglikelihood', 'logjacobian', 'logprior', 'loglr',
                    'lognl')
        maxidx = loglr.argmax()
        for stat, value in vals[maxidx][1].items():
            if stat not in reserved:
                setattr(self._current_stats, stat, value)
        self._current_stats.maxl_mtotal = self.mtotal_grid[maxidx]
        self._current_stats.maxl_loglr = loglr[maxidx]

        return logsumexp(loglr) - numpy.log(len(self.mtotal_grid))

    def reconstruct(self, rec=None, seed=None):
        """Reconstruct the marginalized parameters.

        A total mass is drawn from the grid in proportion to its likelihood,
        and the base model is then asked to reconstruct whatever it
        marginalized over, holding that total mass.
        """
        if seed is not None:
            numpy.random.seed(seed)
        if rec is None:
            rec = {}

        if self.mtotal_grid is None:
            rec.update(self.model.reconstruct())
            return rec

        # evaluate the grid at the current parameters if that has not already
        # happened, then draw a total mass from it
        self.loglr
        xl = draw_sample(self._grid_loglr)
        rec["mtotal"] = self.mtotal_grid[xl]

        # put the base model at the drawn total mass, so that it reconstructs
        # its own marginalized parameters there
        self.model.update(**self._grid_params[xl])
        rec.update(self.model.reconstruct())
        return rec


class BruteParallelGaussianMarginalize(BaseGaussianNoise):
    name = "brute_parallel_gaussian_marginalize"

    def __init__(
        self, variable_params, cores=10, base_model=None, marginalize_phase=None, **kwds
    ):
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
        stats.append("maxl_phase")
        if "maxl_loglr" not in stats:
            stats.append("maxl_loglr")
        return stats

    def _loglr(self):
        if self.phase is not None:
            params = []
            for p in self.phase:
                pref = self.current_params.copy()
                pref["coa_phase"] = p
                params.append(pref)
            vals = list(self.pool.map(self.call, params))
            loglr = numpy.array([v[0] for v in vals])
            # get the maxl values
            if "maxl_loglr" not in self.model._extra_stats:
                maxl_loglrs = loglr
            else:
                maxl_loglrs = numpy.array([v[1]["maxl_loglr"] for v in vals])
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

    def __init__(
        self, variable_params, cores=1, loop_polarization=False, base_model=None, **kwds
    ):
        super().__init__(variable_params, **kwds)

        from pycbc.inference.models import models

        kwds.update(
            models[base_model].extra_args_from_config(
                kwds["config_object"], "model", skip_args=[]
            )
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
        marg_llr = max_llr + math.log(marg_lrfac / self.num_sky_modes)

        # set the stats
        for sym_num in range(self.num_sky_modes):
            setattr(self._current_stats, f"llr_mode_{sym_num}", loglr[sym_num])

        return marg_llr

    def _apply_sky_point_rotation(self, pref, sky_num):
        """Apply the sky point rotation for mode sky_num to parameters pref"""
        lambdal = pref["eclipticlongitude"]
        beta = pref["eclipticlatitude"]
        psi = pref["polarization"]
        inc = pref["inclination"]

        pol_num = sky_num // 8
        sky_num = sky_num % 8
        long_num = sky_num % 4
        lat_num = sky_num // 4

        # Apply latitude symmetry mode
        if lat_num:
            beta = -beta
            inc = numpy.pi - inc
            psi = numpy.pi - psi

        # Apply longitudonal symmetry mode
        lambdal = (lambdal + long_num * 0.5 * numpy.pi) % (2 * numpy.pi)
        psi = (psi + long_num * 0.5 * numpy.pi) % (2 * numpy.pi)

        # Apply additional polarization mode (shouldn't be needed)
        if pol_num:
            psi = psi + (math.pi / 2.0)

        pref["eclipticlongitude"] = lambdal
        pref["eclipticlatitude"] = beta
        pref["polarization"] = psi
        pref["inclination"] = inc

    @classmethod
    def from_config(cls, cp, **kwargs):
        kwargs["config_object"] = cp
        return super().from_config(cp, **kwargs)

    def reconstruct(self, seed=None):
        """Reconstruct a point from unwrapping the 8-fold sky symmetry"""
        if seed is not None:
            numpy.random.seed(seed)
        rec = {}

        logging.info("Reconstruct LISA sky mode symmetry")
        self.reconstruct_sky_points = True
        loglr = self.loglr
        xl = draw_sample(loglr)
        logging.info("Found point %d", xl)
        # Undo rotations
        pref = self.current_params.copy()
        self._apply_sky_point_rotation(pref, xl)

        for val in [
            "polarization",
            "eclipticlongitude",
            "eclipticlatitude",
            "inclination",
        ]:
            rec[val] = pref[val]
        rec["loglr"] = loglr[xl]
        rec["loglikelihood"] = self.lognl + rec["loglr"]
        self.reconstruct_sky_points = False
        return self.model.reconstruct(seed=seed, rec=rec)
