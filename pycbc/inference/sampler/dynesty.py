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

import os
import time
import datetime
import signal
import logging
from pycbc.pool import (choose_pool, is_main_process)
import numpy
import dynesty
from dynesty.utils import resample_equal
from pycbc.inference.io import (DynestyFile, validate_checkpoint_files)
from pycbc.distributions import read_constraints_from_config
from .base import (BaseSampler, setup_output)
from .base_nested_sampler import BaseNestedSampler
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

class DynestySampler(BaseNestedSampler, BaseSampler):
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

    checkpoint : bool
        Enable checkpointing
    checkpoint_interval : int
        Checkpoint after every so many iterations
    checkpoint_period : int
        Checkpoint only after these many seconds since last checkpoint
    checkpoint_on_signal : bool
        Checkpoint if receiving external signals
    checkpoint_signal : str
        Signal to send when checkpointing

    nprocesses : int
        Number of processes to use in MPI pool
    use_mpi : bool
        Enable / disable MPI usage in multiprocessing pool
    run_kwds : dict
        Keyword options for dynesty's sampling routine
    kwargs : dict

    """
    name = "dynesty"
    _io = DynestyFile

    def __init__(self, model, nlive,
                 loglikelihood_function=None,
                 checkpoint=True,
                 checkpoint_interval=1000,
                 checkpoint_period=600,
                 checkpoint_on_signal=True,
                 checkpoint_signal=None,
                 nprocesses=1, use_mpi=False,
                 run_kwds=None, **kwargs):
        self.model = model
        log_likelihood_call, prior_call = setup_calls(
            model,
            nprocesses=nprocesses,
            loglikelihood_function=loglikelihood_function)
        # Set up the pool
        self.use_mpi = use_mpi
        self.nprocesses = nprocesses
        pool = choose_pool(mpi=use_mpi, processes=nprocesses)
        if pool is not None:
            pool.size = nprocesses

        self.run_kwds = {} if run_kwds is None else run_kwds
        self.nlive = nlive
        self.names = model.sampling_params
        self.ndim = len(model.sampling_params)

        self._checkpoint = checkpoint
        self._checkpoint_interval = checkpoint_interval
        self._checkpoint_signal = checkpoint_signal
        self._checkpoint_on_signal = checkpoint_on_signal
        self._checkpoint_period = checkpoint_period
        if self._checkpoint:
            logging.info("Checkpointing every {0} iters, {1} secs".format(
                self.checkpoint_interval, self.checkpoint_period))

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
        self._is_main_process = is_main_process()

    def set_sampler_specific_state_from_file(self, filename):
        """Set state of sampler object that could not be serialized"""
        self._sampler.rstate = numpy.random
        self._sampler.nqueue = -1
        logging.info("Creating a new multiprocessing pool")
        pool = choose_pool(mpi=self.use_mpi, processes=self.nprocesses)
        if pool is not None:
            pool.size = self.nprocesses
        self._sampler.pool = pool
        self._sampler.M = pool.map

    def getstate(self):
        """Strips unserializable attributes of sampler
        """
        self_dict = self._sampler.__dict__.copy()
        # Remove multiprocessing / mpi pool-related attributes
        # as they cannot be pickled
        del self_dict['pool']
        del self_dict['M']
        return self_dict

    def run(self):
        """Run the sampler."""
        if self._checkpoint_on_signal:
            # Enable on-signal checkpointing:
            # This purely guards against interruptions on compute clusters
            signal.signal(signal.SIGTERM, self.checkpoint_on_signal)
            signal.signal(signal.SIGALRM, self.checkpoint_on_signal)
            signal.signal(signal.SIGQUIT, self.checkpoint_on_signal)
            signal.signal(signal.SIGTSTP, self.checkpoint_on_signal)
            signal.signal(signal.SIGINT, self.checkpoint_on_signal)
            signal.signal(signal.SIGUSR1, self.checkpoint_on_signal)
            signal.signal(signal.SIGUSR2, self.checkpoint_on_signal)

        # Instead of calling the `run_nested` wrapper, we invoke
        # the core nested sampling loop for the dynesty sampler
        t0 = datetime.datetime.now()
        sampling_time = 0
        for it, res in enumerate(self._sampler.sample(**self.run_kwds)):
            i = it - 1
            dynesty.results.print_fn_fallback(res, i, self._sampler.ncall,
                                              dlogz=self.run_kwds['dlogz'])
            if not self._checkpoint or it == 0:
                continue
            if self.checkpoint_interval is not None:
                if it % self.checkpoint_interval != 0:
                    continue
            # Checkpoint periodically
            sampling_time += (datetime.datetime.now() - t0).total_seconds()
            t0 = datetime.datetime.now()
            if os.path.isfile(self.checkpoint_file):
                time_since_last_checkpoint_s =\
                    time.time() - os.path.getmtime(self.checkpoint_file)
            else:
                time_since_last_checkpoint_s = numpy.inf
            if time_since_last_checkpoint_s > self.checkpoint_period:
                self._sampler.kwargs['sampling_time'] = sampling_time
                if self.is_main_process:
                    self.checkpoint()

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
                 'sample': str}
        extra = {}
        run_extra = {}
        for karg in cargs:
            if cp.has_option(section, karg):
                extra[karg] = cargs[karg](cp.get(section, karg))

        for karg in rargs:
            if cp.has_option(section, karg):
                run_extra[karg] = rargs[karg](cp.get(section, karg))

        # get the checkpoint interval, if it's specified
        checkpoint_interval = cls.checkpoint_from_config(cp, section)
        checkpoint_signal = cls.ckpt_signal_from_config(cp, section)
        checkpoint_period = get_optional_arg_from_config(cp, section,
                                                         'checkpoint-period',
                                                         dtype=int)
        checkpoint_on_signal =\
            get_optional_arg_from_config(cp, section, 'checkpoint-on-signal',
                                         dtype=bool)
        checkpoint = get_optional_arg_from_config(cp, section,
                                                  'checkpoint', dtype=bool)
        if checkpoint is None and ((checkpoint_interval is not None) or
                                   (checkpoint_period is not None) or
                                   checkpoint_on_signal):
            checkpoint = True

        obj = cls(model, nlive=nlive, nprocesses=nprocesses,
                  loglikelihood_function=loglikelihood_function,
                  checkpoint=checkpoint,
                  checkpoint_interval=checkpoint_interval,
                  checkpoint_period=checkpoint_period,
                  checkpoint_signal=checkpoint_signal,
                  checkpoint_on_signal=checkpoint_on_signal,
                  use_mpi=use_mpi, run_kwds=run_extra, **extra)
        setup_output(obj, output_file)
        if not obj.new_checkpoint:
            obj.resume_from_checkpoint()
        return obj

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
