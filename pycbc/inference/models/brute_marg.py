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
        return _model.loglr

class BruteParallelGaussianMarginalize(BaseGaussianNoise):
    name = "brute_parallel_gaussian_marginalize"

    def __init__(self, variable_params,
                 cores=10,
                 base_model=None,
                 marginalize_phase=None,
                 **kwds):
        super(BruteParallelGaussianMarginalize, self).__init__(variable_params,
                                     **kwds)

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

    def _loglr(self):
        if self.phase is not None:
            params = []
            for p in self.phase:
                pref = self.current_params.copy()
                pref['coa_phase'] = p
                params.append(pref)
            loglr = numpy.array(list(self.pool.map(self.call, params)))
            return logsumexp(loglr) - numpy.log(len(self.phase))
