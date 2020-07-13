# Copyright (C) 2020 Sumit Kumar, Alex Nitz
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


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""
Common utilities for samplers that rely on transforming between a unit cube
and the prior space. This is typical of many nested sampling algorithms.
"""
import numpy

from .. import models


def call_global_loglikelihood(cube):
    return models._global_instance.log_likelihood(cube)


def call_global_logprior(cube):
    return models._global_instance.prior_transform(cube)


def setup_calls(model, loglikelihood_function=None, copy_prior=False):
    """ Configure calls for MPI support
    """
    model_call = CubeModel(model, loglikelihood_function,
                           copy_prior=copy_prior)

    # these are used to help paralleize over multiple cores / MPI
    models._global_instance = model_call
    log_likelihood_call = call_global_loglikelihood
    prior_call = call_global_logprior
    return log_likelihood_call, prior_call


class CubeModel(object):
    """ Class for making PyCBC Inference 'model class'

    Parameters
    ----------
    model : inference.BaseModel instance
             A model instance from pycbc.
    """

    def __init__(self, model, loglikelihood_function=None, copy_prior=False):
        if model.sampling_transforms is not None:
            raise ValueError("Ultranest or dynesty do not support sampling transforms")
        self.model = model
        if loglikelihood_function is None:
            loglikelihood_function = 'loglikelihood'
        self.loglikelihood_function = loglikelihood_function
        self.copy_prior = copy_prior

    def log_likelihood(self, cube):
        """
        returns log likelihood function
        """
        params = dict(zip(self.model.sampling_params, cube))
        self.model.update(**params)
        if self.model.logprior == -numpy.inf:
            return -numpy.inf
        return getattr(self.model, self.loglikelihood_function)

    def prior_transform(self, cube):
        """
        prior transform function for ultranest sampler
        It takes unit cube as input parameter and apply
        prior transforms
        """
        if self.copy_prior:
            cube = cube.copy()

        # we preserve the type of cube to whatever we were given
        dict_cube = dict(zip(self.model.variable_params, cube))
        inv = self.model.prior_distribution.cdfinv(**dict_cube)
        for i, param in enumerate(self.model.variable_params):
            cube[i] = inv[param]
        return cube
