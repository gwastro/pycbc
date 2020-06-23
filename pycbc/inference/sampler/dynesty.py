# Copyright (C) 2019  Collin Capano, Sumit Kumar, Prayush Kumar
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
import os
import time
from pycbc.pool import choose_pool
import numpy
import dynesty
from dynesty.utils import resample_equal
from pycbc.inference.io import (DynestyFile, validate_checkpoint_files,
                                loadfile)
from pycbc.distributions import read_constraints_from_config
from .base import (BaseSampler, setup_output)
from .base_mcmc import get_optional_arg_from_config
from .base_cube import setup_calls
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
    pool : function with map, Optional
        A provider of a map function that allows a function call to be run
        over multiple sets of arguments and possibly maps them to
        cores/nodes/etc.
    """
    name = "dynesty"
    _io = DynestyFile

    def __init__(self, model, nlive, nprocesses=1, niter=0,
                 checkpoint_time_interval=None, maxcall=None,
                 loglikelihood_function=None, use_mpi=False,
                 run_kwds=None, **kwargs):
        self.model = model
        log_likelihood_call, prior_call = setup_calls(
            model,
            nprocesses=nprocesses,
            loglikelihood_function=loglikelihood_function)
        # Set up the pool
        pool = choose_pool(mpi=use_mpi, processes=nprocesses)
        if pool is not None:
            pool.size = nprocesses
        self.use_mpi = use_mpi
        self.nprocesses = nprocesses
        self.maxcall = maxcall
        self.checkpoint_time_interval = checkpoint_time_interval
        self.run_kwds = {} if run_kwds is None else run_kwds
        self.nlive = nlive
        self.names = model.sampling_params
        self.ndim = len(model.sampling_params)
        self.checkpoint_file = None
        self.__getstate__ = self.getstate
        self.niter = niter
        # Enable checkpointing if checkpoint_time_interval is set in config
        # file in sampler section
        if self.checkpoint_time_interval:
            self.run_with_checkpoint = True
        else:
            self.run_with_checkpoint = False

        if self.nlive < 0:
            # Interpret a negative input value for the number of live points
            # (which is clearly an invalid input in all senses)
            # as the desire to dynamically determine that number
            self._sampler = dynesty.DynamicNestedSampler(log_likelihood_call,
                                                         prior_call, self.ndim,
                                                         pool=pool, **kwargs)
        else:
            self._sampler = dynesty.NestedSampler(log_likelihood_call,
                                                  prior_call, self.ndim,
                                                  nlive=self.nlive,
                                                  pool=pool, **kwargs)
        self._sampler.pool = pool
        self._sampler.M = pool.map

    def run(self):
        diff_niter = 1
        if self.run_with_checkpoint is True:
            if self.checkpoint_file:
                try:
                    self.resume_from_checkpoint()
                    self.niter = self._sampler.results.niter
                    logging.info('Successfully read from the checkpoint file')
                except KeyError:
                    logging.info("Could not read checkpoint file")
            else:
                self.niter = 0
            n_checkpointing = 1
            t0 = time.time()
            while diff_niter != 0:
                delta_t = time.time()-t0
                self._sampler.run_nested(**self.run_kwds)
                diff_niter = self._sampler.results.niter - self.niter
                if delta_t >= self.checkpoint_time_interval:
                    logging.info('Checkpointing N={}'.format(n_checkpointing))
                    self.checkpoint()
                    n_checkpointing += 1
                    t0 = time.time()
                self.niter = self._sampler.results.niter
        else:
            self._sampler.run_nested(**self.run_kwds)


    @property
    def io(self):
        return self._io

    @property
    def niterations(self):
        return len(tuple(self.samples.values())[0])

    @classmethod
    def from_config(cls, cp, model, output_file=None, nprocesses=1,
                    use_mpi=False, loglikelihood_function=None):
        """
        Loads the sampler from the given config file.
        """
        section = "sampler"
        # check name
        assert cp.get(section, "name") == cls.name, (
            "name in section [sampler] must match mine")
        # get the number of live points to use
        nlive = int(cp.get(section, "nlive"))
        loglikelihood_function = \
            get_optional_arg_from_config(cp, section, 'loglikelihood-function')

        # optional run_nested arguments for dynesty
        rargs = {'maxiter': int,
                 'maxcall': int,
                 'dlogz': float,
                 'logl_max': float,
                 'n_effective': int,
                 }

        # optional arguments for dynesty
        cargs = {'bound': str,
                 'bootstrap': int,
                 'enlarge': float,
                 'update_interval': float,
                 'sample': str,
                 'checkpoint_time_interval': int
                 }
        extra = {}
        run_extra = {}
        for karg in cargs:
            if cp.has_option(section, karg):
                extra[karg] = cargs[karg](cp.get(section, karg))

        for karg in rargs:
            if cp.has_option(section, karg):
                run_extra[karg] = rargs[karg](cp.get(section, karg))

        obj = cls(model, nlive=nlive, nprocesses=nprocesses,
                  loglikelihood_function=loglikelihood_function,
                  use_mpi=use_mpi, run_kwds=run_extra, **extra)
        setup_output(obj, output_file, check_nsamples=False)
        if not obj.new_checkpoint:
            obj.resume_from_checkpoint()
        return obj

    def checkpoint(self):
        """Checkpoint function for dynesty sampler
        """
        self.__getstate__ = self.getstate()
        pickle_obj = self._sampler
        for fn in [self.checkpoint_file]:
            with loadfile(fn, 'a') as fp:
                fp.write_pickled_data_into_checkpoint_file(pickle_obj)
                fp.write_random_state()

    def resume_from_checkpoint(self):
        for fn in [self.checkpoint_file]:
            with loadfile(fn, 'r') as fp:
                self._sampler = fp.read_pickled_data_from_checkpoint_file()
            self._sampler.rstate = numpy.random
            self.set_state_from_file(fn)
            pool = choose_pool(mpi=self.use_mpi, processes=self.nprocesses)
            if pool is not None:
                pool.size = self.nprocesses
            self._sampler.M = pool.map
            self._sampler.pool = pool

    def set_state_from_file(self, filename):
        """Sets the state of the sampler back to the instance saved in a file.
        """
        with self.io(filename, 'r') as fp:
            rstate = fp.read_random_state()
        # set the numpy random state
        numpy.random.set_state(rstate)
        self._sampler.rstate.set_state(rstate)

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
            self.checkpoint_file, self.backup_file, check_nsamples=False)
        if not checkpoint_valid:
            raise IOError("error writing to checkpoint file")

    def getstate(self):
        """Strips unserializable attributes of sampler
        """
        self_dict = self._sampler.__dict__.copy()
        # Remove multiprocessing / mpi pool-related/ rstate attributes
        # as they cannot be pickled
        del self_dict['pool']
        del self_dict['M']
        del self_dict['rstate']
        del self_dict['loglikelihood']
        del self_dict['prior_transform']
        return self_dict

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
            fp.write_random_state()

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
