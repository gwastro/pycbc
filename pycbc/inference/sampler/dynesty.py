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
import copy
import os
import time
import numpy
import dynesty
from pycbc.pool import choose_pool
from dynesty import utils as dyfunc
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

    def __init__(self, model, nlive, nprocesses=1,
                 checkpoint_time_interval=None, maxcall=None,
                 loglikelihood_function=None, use_mpi=False,
                 run_kwds=None, **kwargs):

        self.model = model
        log_likelihood_call, prior_call = setup_calls(
            model,
            loglikelihood_function=loglikelihood_function)
        # Set up the pool
        self.pool = choose_pool(mpi=use_mpi, processes=nprocesses)

        self.maxcall = maxcall
        self.checkpoint_time_interval = checkpoint_time_interval
        self.run_kwds = {} if run_kwds is None else run_kwds
        self.nlive = nlive
        self.names = model.sampling_params
        self.ndim = len(model.sampling_params)
        self.checkpoint_file = None
        # Enable checkpointing if checkpoint_time_interval is set in config
        # file in sampler section
        if self.checkpoint_time_interval:
            self.run_with_checkpoint = True
            if self.maxcall is None:
                self.maxcall = 5000 * self.pool.size
            logging.info("Checkpointing enabled, will verify every %s calls"
                         " and try to checkpoint every %s seconds",
                         self.maxcall, self.checkpoint_time_interval)
        else:
            self.run_with_checkpoint = False

        # Check for cyclic boundaries
        periodic = []
        cyclic = self.model.prior_distribution.cyclic
        for i, param in enumerate(self.variable_params):
            if param in cyclic:
                logging.info('Param: %s will be cyclic', param)
                periodic.append(i)

        if len(periodic) == 0:
            periodic = None

        # Check for reflected boundaries. Dynesty only supports
        # reflection on both min and max of boundary.
        reflective = []
        reflect = self.model.prior_distribution.well_reflected
        for i, param in enumerate(self.variable_params):
            if param in reflect:
                logging.info("Param: %s will be well reflected", param)
                reflective.append(i)

        if len(reflective) == 0:
            reflective = None

        if self.nlive < 0:
            # Interpret a negative input value for the number of live points
            # (which is clearly an invalid input in all senses)
            # as the desire to dynamically determine that number
            self._sampler = dynesty.DynamicNestedSampler(log_likelihood_call,
                                                         prior_call, self.ndim,
                                                         pool=self.pool,
                                                         reflective=reflective,
                                                         periodic=periodic,
                                                         **kwargs)
            self.run_with_checkpoint = False
            logging.info("Checkpointing not currently supported with"
                         "DYNAMIC nested sampler")
        else:
            self._sampler = dynesty.NestedSampler(log_likelihood_call,
                                                  prior_call, self.ndim,
                                                  nlive=self.nlive,
                                                  reflective=reflective,
                                                  periodic=periodic,
                                                  pool=self.pool, **kwargs)

        # properties of the internal sampler which should not be pickled
        self.no_pickle = ['loglikelihood',
                          'prior_transform',
                          'propose_point',
                          'update_proposal',
                          '_UPDATE', '_PROPOSE',
                          'evolve_point']

    def run(self):
        diff_niter = 1
        if self.run_with_checkpoint is True:
            n_checkpointing = 1
            t0 = time.time()
            it = self._sampler.it

            logging.info('Starting from iteration: %s', it)
            while diff_niter != 0:
                self._sampler.run_nested(maxcall=self.maxcall, **self.run_kwds)

                delta_t = time.time() - t0
                diff_niter = self._sampler.it - it
                logging.info("Checking if we should checkpoint: %.2f s", delta_t)

                if delta_t >= self.checkpoint_time_interval:
                    logging.info('Checkpointing N={}'.format(n_checkpointing))
                    self.checkpoint()
                    n_checkpointing += 1
                    t0 = time.time()
                it = self._sampler.it
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
                 'dlogz': float,
                 'logl_max': float,
                 'n_effective': int,
                 }

        # optional arguments for dynesty
        cargs = {'bound': str,
                 'maxcall': int,
                 'bootstrap': int,
                 'enlarge': float,
                 'update_interval': float,
                 'sample': str,
                 'checkpoint_time_interval': float
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
        with loadfile(self.checkpoint_file, 'a') as fp:
            fp.write_random_state()

            # Dynesty has its own __getstate__ which deletes
            # random state information and the pool
            saved = {}
            for key in self.no_pickle:
                if hasattr(self._sampler, key):
                    saved[key] = getattr(self._sampler, key)
                    setattr(self._sampler, key, None)

            # For the dynamic sampler, we must also handle the internal
            # sampler object
            #saved_internal = {}
            #if self.nlive < 0:
            #    for key in self.no_pickle:
            #        if hasattr(self._sampler.sampler, key):
            #            saved[key] = getattr(self._sampler.sampler, key)
            #            setattr(self._sampler.sampler, key, None)

            #for key in self._sampler.__dict__:
            #    print(key, type(self._sampler.__dict__[key]))

            #for key in self._sampler.sampler.__dict__:
            #    print(key, type(self._sampler.sampler.__dict__[key]))

            fp.write_pickled_data_into_checkpoint_file(self._sampler)

        # Restore properties that couldn't be pickled if we are continuing
        for key in saved:
            setattr(self._sampler, key, saved[key])

        # Restore for dynamic nested sampler's internal sampler
        #for key in saved_internal:
        #    setattr(self._sampler.sampler, key, saved_internal[key])

    def resume_from_checkpoint(self):
        try:
            with loadfile(self.checkpoint_file, 'r') as fp:
                sampler = fp.read_pickled_data_from_checkpoint_file()

                for key in sampler.__dict__:
                    if key not in self.no_pickle:
                        value = getattr(sampler, key)
                        setattr(self._sampler, key, value)

                # If dynamic sampling, also restore internal sampler
                #if self.nlive < 0:
                #    for key in sampler.__dict__:
                #       if key not in self.no_pickle:
                #           value = getattr(sampler.sampler, key)
                #           setattr(self._sampler.sampler, key, value)

            self.set_state_from_file(self.checkpoint_file)
            logging.info("Found valid checkpoint file: %s",
                         self.checkpoint_file)
        except Exception as e:
            print(e)
            logging.info("Failed to load checkpoint file")

    def set_state_from_file(self, filename):
        """Sets the state of the sampler back to the instance saved in a file.
        """
        with self.io(filename, 'r') as fp:
            numpy.random.set_state(fp.read_random_state())

        self._sampler.rstate = numpy.random
        #if self.nlive < 0:
        #    self._sampler.sampler.rstate = numpy.random

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

    @property
    def samples(self):
        results = self._sampler.results
        samples = results.samples
        weights = numpy.exp(results.logwt - results.logz[-1])
        N = len(weights)
        positions = (numpy.random.random() + numpy.arange(N)) / N

        idx = numpy.zeros(N, dtype=numpy.int)
        cumulative_sum = numpy.cumsum(weights)
        i, j = 0, 0
        while i < N:
            if positions[i] < cumulative_sum[j]:
                idx[i] = j
                i += 1
            else:
                j += 1

        numpy.random.shuffle(idx)
        post = {'loglikelihood': self._sampler.results.logl[idx]}
        for i, param in enumerate(self.variable_params):
            post[param] = samples[:, i][idx]
        return post

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
            fp.write_samples(self.samples)

            # write log evidence
            fp.write_logevidence(self._sampler.results.logz[-1:][0],
                                 self._sampler.results.logzerr[-1:][0])

    @property
    def model_stats(self):
        pass

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
