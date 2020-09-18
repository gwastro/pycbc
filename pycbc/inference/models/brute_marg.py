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
from multiprocessing import Pool

from .gaussian_noise import BaseGaussianNoise
from scipy.miscm import logsumexp

class BruteParallelGaussianMarginalize(BaseGaussianModel):
    name = "brute_parallel_gaussian_marginalize"

    def __init__(self, variable_params,
                 cores=10,
                 base_model=None,
                 **kwds):
        super(BruteParallelMarginalize, self).__init__(variable_params, **kwds)

        from pycbc.inference.models import models
        self.model = models[base_model](variable_params, **kwds)

        def likelihood_wrapper(self, params):
            self.model.update(**params)
            return self.model.loglr

        self.call = likelihood_wrapper

        # size of pool for each likelihood call
        self.pool = Pool(int(cores))

        # Only one for now, but can be easily extended
        self.phase = None
        if 'marginalize_phase' in kwds:
            samples = kwds['marginalize_phase_samples']
            self.phase = numpy.linspace(0, 2.0 * numpy.pi, int(samples))

    def _loglr(self):
        if self.phase:
            params = []
            for p in self.phase:
                pref = self.current_params.copy()
                pref['coa_phase'] = p
                params.append(pref)
            loglr = numpy.array(list(self.pool.map(self.call, params)))
            return logsumexp(loglr) - numpy.log(len(self.phase))
