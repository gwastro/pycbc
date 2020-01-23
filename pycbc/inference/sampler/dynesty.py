# Copyright (C) 2019  Collin Capano, Sumit Kumar
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
This modules provides classes and functions for using the dynesty sampler
packages for parameter estimation.
"""


from __future__ import absolute_import

import logging
from pycbc.pool import choose_pool
import numpy
import dynesty
from dynesty.utils import resample_equal
from pycbc.inference.io import (DynestyFile, validate_checkpoint_files)
from pycbc.distributions import read_constraints_from_config
from .base import (BaseSampler, setup_output)
from .base_mcmc import get_optional_arg_from_config
from .. import models


#
# =============================================================================
#
#                                   Samplers
#
# =============================================================================
#

class DynestySampler(BaseSampler):
    """This class is used to construct an Dynesty sampler from the dynesty
    package.

    Parameters
    ----------
    model : model
        A model from ``pycbc.inference.models``.
    nlive : int
        Number of live points to use in sampler.
    dlogz: float
        Tolerance limit to the value of logz (also a convergence criteria)
    pool : function with map, Optional
        A provider of a map function that allows a function call to be run
        over multiple sets of arguments and possibly maps them to
        cores/nodes/etc.
    """
    name = "dynesty"
    _io = DynestyFile

    def __init__(self, model, nlive, dlogz, nprocesses=1,
                 loglikelihood_function=None, use_mpi=False, **kwargs):
        self.model = model
        # Set up the pool
        model_call = DynestyModel(model, loglikelihood_function)
        if nprocesses > 1:
            # these are used to help paralleize over multiple cores / MPI
            models._global_instance = model_call
            log_likelihood_call = _call_global_loglikelihood
            prior_call = _call_global_logprior
        else:
            prior_call = model_call.prior_transform
            log_likelihood_call = model_call.log_likelihood
        pool = choose_pool(mpi=use_mpi, processes=nprocesses)
        if pool is not None:
            pool.size = nprocesses

        self.nlive = nlive
        self.dlogz = dlogz
        self.names = model.sampling_params
        self.ndim = len(model.sampling_params)
        self.checkpoint_file = None
        self._sampler = dynesty.NestedSampler(log_likelihood_call,
                                              prior_call, self.ndim,
                                              nlive=self.nlive,
                                              dlogz=self.dlogz,
                                              pool=pool, **kwargs)

    def run(self):
        self._sampler.run_nested()

    @property
    def io(self):
        return self._io

    @property
    def niterations(self):
        return len(tuple(self.samples.values())[0])

    @classmethod
    def from_config(cls, cp, model, nprocesses=1, loglikelihood_function=None,
                    use_mpi=False):
        """
        Loads the sampler from the given config file.
        """
        section = "sampler"
        # check name
        assert cp.get(section, "name") == cls.name, (
            "name in section [sampler] must match mine")
        # get the number of live points to use
        nlive = int(cp.get(section, "nlive"))
        dlogz = float(cp.get(section, "dlogz"))
        loglikelihood_function = \
            get_optional_arg_from_config(cp, section, 'loglikelihood-function')

        # optional arguments for dynesty
        cargs = {'bound': str,
                 'bootstrap': int,
                 'enlarge': float,
                 'update_interval': float,
                 'sample': str}
        extra = {}
        for karg in cargs:
            if cp.has_option(section, karg):
                extra[karg] = cargs[karg](cp.get(section, karg))

        obj = cls(model, nlive=nlive, dlogz=dlogz, nprocesses=nprocesses,
                  loglikelihood_function=loglikelihood_function,
                  use_mpi=use_mpi, **extra)
        return obj

    def checkpoint(self):
        pass

    def finalize(self):
        logz = self._sampler.results.logz[-1:][0]
        dlogz = self._sampler.results.logzerr[-1:][0]
        logging.info("log Z, dlog Z: {}, {}".format(logz, dlogz))
        for fn in [self.checkpoint_file]:
            with self.io(fn, "a") as fp:
                fp.write_logevidence(logz, dlogz)
        logging.info("Writing samples to files")
        for fn in [self.checkpoint_file, self.backup_file]:
            self.write_results(fn)
        logging.info("Validating checkpoint and backup files")
        checkpoint_valid = validate_checkpoint_files(
            self.checkpoint_file, self.backup_file)
        if not checkpoint_valid:
            raise IOError("error writing to checkpoint file")

    @property
    def model_stats(self):
        logl = self._sampler.results.logl
        return {'loglikelihood': logl}

    @property
    def samples(self):
        samples_dict = {p: self.posterior_samples[:, i] for p, i in
                        zip(self.model.variable_params, range(self.ndim))}
        return samples_dict

    def set_initial_conditions(self, initial_distribution=None,
                               samples_file=None):
        """Sets up the starting point for the sampler.

        Should also set the sampler's random state.
        """
        pass

    def write_results(self, filename):
        """Writes samples, model stats, acceptance fraction, and random state
        to the given file.

        Parameters
        -----------
        filename : str
            The file to write to. The file is opened using the ``io`` class
            in an an append state.
        """
        with self.io(filename, 'a') as fp:
            # write samples
            fp.write_samples(self.samples, self.model.variable_params)
            # write stats
            fp.write_samples(self.model_stats)
            # write log evidence
            fp.write_logevidence(self._sampler.results.logz[-1:][0],
                                 self._sampler.results.logzerr[-1:][0])

    def setup_output(self, output_file, force=False):
        """Sets up the sampler's checkpoint and output files.

        The checkpoint file has the same name as the output file, but with
        ``.checkpoint`` appended to the name. A backup file will also be
        created.

        If the output file already exists, an ``OSError`` will be raised.
        This can be overridden by setting ``force`` to ``True``.

        Parameters
        ----------
        sampler : sampler instance
            Sampler
        output_file : str
            Name of the output file.
        force : bool, optional
            If the output file already exists, overwrite it.
        """
        setup_output(self, output_file, force=force)

    @property
    def posterior_samples(self):
        """
        Returns posterior samples from nested samples and weights
        given by dynsety sampler
        """

        dynesty_samples = self._sampler.results['samples']
        wt = numpy.exp(self._sampler.results['logwt'] - 
                       self._sampler.results['logz'][-1])
        # Make sure that sum of weights equal to 1
        weights = wt/numpy.sum(wt)
        posterior_dynesty = resample_equal(dynesty_samples, weights)
        return posterior_dynesty

    @property
    def logz(self):
        """
        return bayesian evidence estimated by
        dynesty sampler
        """

        return self._sampler.results.logz[-1:][0]

    @property
    def logz_err(self):
        """
        return error in bayesian evidence estimated by
        dynesty sampler
        """

        return self._sampler.results.logzerr[-1:][0]


def _call_global_loglikelihood(cube):
    return models._global_instance.log_likelihood(cube)


def _call_global_logprior(cube):
    return models._global_instance.prior_transform(cube)


class DynestyModel(object):
    """
    Class for making PyCBC Inference 'model class'
    Parameters
    ----------
    model : inference.BaseModel instance
             A model instance from pycbc.
    """

    def __init__(self, model, loglikelihood_function=None):
        if model.sampling_transforms is not None:
            raise ValueError("Dinesty does not support sampling transforms")
        self.model = model
        if loglikelihood_function is None:
            loglikelihood_function = 'loglikelihood'
        self.loglikelihood_function = loglikelihood_function

    def log_likelihood(self, cube):
        """
        returns log likelihood function
        """
        params = {p: v for p, v in zip(self.model.sampling_params, cube)}
        self.model.update(**params)
        if self.model.logprior == -numpy.inf:
            return -numpy.inf
        return getattr(self.model, self.loglikelihood_function)

    def prior_transform(self, cube):
        """
        prior transform function for dynesty sampler
        It takes unit cube as input parameter and apply
        prior transforms
        """
        prior_dists = self.model.prior_distribution.distributions
        dist_dict = {}
        for dist in prior_dists:
            dist_dict.update({param: dist for param in dist.params})
        for i, param in enumerate(self.model.variable_params):
            cube[i] = dist_dict[param].cdfinv(param, cube[i])
        return cube
