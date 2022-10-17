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
import numpy

from multiprocessing import Pool
from .gaussian_noise import BaseGaussianNoise
from scipy.special import logsumexp

_model = None
class likelihood_wrapper(object):
    def __init__(self, model):
        global _model
        _model = model

    def __call__(self, params):
        global _model
        _model.update(**params)
        loglr = _model.loglr
        return loglr, _model.current_stats

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


    @property
    def _extra_stats(self):
        stats = self.model._extra_stats
        return stats

    def _loglr(self):
        params = []
        for sym_num in range(8):
            # FIXME: Literature says this should be 8 points, but I always see
            #        two polarization modes if I use the indicated 8. Possibly
            #        there is some issue here, and I'm not using the *right*
            #        8 points. Note that what I think are the right 8 points
            #        are acheived by setting this to range(8).
            pref = self.current_params.copy()
            lambdal = pref['eclipticlongitude']
            beta = pref['eclipticlatitude']
            psi = pref['polarization']
            inc = pref['inclination']

            pol_num = sym_num // 8
            sym_num = sym_num % 8
            long_num = sym_num % 4
            lat_num = sym_num // 4

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
            params.append(pref)

        vals = list(self.mapfunc(self.call, params))
        loglr = numpy.array([v[0] for v in vals])

        max_llr_idx = loglr.argmax()
        max_llr = loglr[max_llr_idx]
        marg_lrfac = sum([math.exp(llr - max_llr) for llr in loglr])
        marg_llr = max_llr + math.log(marg_lrfac/8.)

        # set the stats
        for sym_num in range(8):
            setattr(self._current_stats, f'llr_mode_{sym_num}', loglr[sym_num])

        return marg_llr

    @classmethod
    def from_config(cls, cp, **kwargs):
        kwargs['config_object'] = cp
        return super(BruteLISASkyModesMarginalize, cls).from_config(cp, **kwargs)
