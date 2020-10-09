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
import dynesty, dynesty.dynesty, dynesty.nestedsamplers
from dynesty.utils import unitcheck, reflect
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

        if 'sample' in kwargs:
            if 'rwalk2' in kwargs['sample']:
                dynesty.dynesty._SAMPLING["rwalk"] = sample_rwalk_mod
                dynesty.nestedsamplers._SAMPLING["rwalk"] = sample_rwalk_mod
                kwargs['sample'] = 'rwalk'

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
        """Loads the sampler from the given config file. Many options are
        directly passed to the underlying dynesty sampler, see the official
        dynesty documentation for more details on these.

        The following options are retrieved in the ``[sampler]`` section:
        * ``name = STR`` :
            Required. This must match the sampler's name.
        * ``maxiter = INT``:
            The maximum number of iterations to run.
        * ``dlogz = FLOAT``:
            The target dlogz stopping condition.
        * ``logl_max = FLOAT``:
            The maximum logl stopping condition.
        * ``n_effective = INT``:
            Target effective number of samples stopping condition
        * ``samples = STR``:
            The method to sample the space. Should be one of 'uniform',
            'rwalk', 'rwalk2' (a modified version of rwalk), or 'slice'.
        * ``first_update_min_ncall = INT``:
            The minimum number of calls before updating the bounding region
            for the first time.
        * ``first_update_min_neff = FLOAT``:
            Don't update the the bounding region untill the efficiency drops
            below this value.
        * ``bound = STR``:
            The method of bounding of the prior volume.
            Should be one of 'single', 'balls', 'cubes', 'multi' or 'none'.
        * ``update_interval = INT``:
            Number of iterations between updating the bounding regions
        * ``enlarge = FLOAT``:
            Factor to enlarge the bonding region.
        * ``bootstrap = INT``:
            The number of bootstrap iterations to determine the enlargement
            factor.
        * ``maxcall = INT``:
            The maximum number of calls before checking if we should checkpoint
        * ``checkpoint_time_interval`` :
            Sets the time in seconds between checkpointing.
        * ``loglikelihood-function`` :
            The attribute of the model to use for the loglikelihood. If
            not provided, will default to ``loglikelihood``.

        Parameters
        ----------
        cp : WorkflowConfigParser instance
            Config file object to parse.
        model : pycbc.inference.model.BaseModel instance
            The model to use.
        output_file : str, optional
            The name of the output file to checkpoint and write results to.
        nprocesses : int, optional
            The number of parallel processes to use. Default is 1.
        use_mpi : bool, optional
            Use MPI for parallelization. Default is False.

        Returns
        -------
        DynestySampler :
            The sampler instance.
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
                 'checkpoint_time_interval': float,
                 'first_update_min_ncall': int,
                 'first_update_min_eff': float,
                 }
        extra = {}
        run_extra = {}
        for karg in cargs:
            if cp.has_option(section, karg):
                extra[karg] = cargs[karg](cp.get(section, karg))

        #This arg needs to be a dict
        first_update = {}
        if 'first_update_min_ncall' in extra:
            first_update['min_ncall'] = extra.pop('first_update_min_ncall')
            logging.info('First update: min_ncall:%s',
                         first_update['min_ncall'])
        if 'first_update_min_eff' in extra:
            first_update['min_eff'] = extra.pop('first_update_min_eff')
            logging.info('First update: min_eff:%s', first_update['min_eff'])
        extra['first_update'] = first_update

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
        # Dynesty has its own __getstate__ which deletes
        # random state information and the pool
        saved = {}
        for key in self.no_pickle:
            if hasattr(self._sampler, key):
                saved[key] = getattr(self._sampler, key)
                setattr(self._sampler, key, None)
        for fn in [self.checkpoint_file, self.backup_file]:
            with self.io(fn, "a") as fp:
                # Write random state
                fp.write_random_state()

                # Write pickled data
                fp.write_pickled_data_into_checkpoint_file(self._sampler)

                # Write nested samples
                fp.write_raw_samples(self.samples)

                # Write logz and dlogz
                logz = self._sampler.results.logz[-1:][0]
                dlogz = self._sampler.results.logzerr[-1:][0]
                fp.write_logevidence(logz, dlogz)

        # Restore properties that couldn't be pickled if we are continuing
        for key in saved:
            setattr(self._sampler, key, saved[key])

    def resume_from_checkpoint(self):
        try:
            with loadfile(self.checkpoint_file, 'r') as fp:
                sampler = fp.read_pickled_data_from_checkpoint_file()

                for key in sampler.__dict__:
                    if key not in self.no_pickle:
                        value = getattr(sampler, key)
                        setattr(self._sampler, key, value)

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
        self.checkpoint()
        logging.info("Validating checkpoint and backup files")
        checkpoint_valid = validate_checkpoint_files(
            self.checkpoint_file, self.backup_file, check_nsamples=False)
        if not checkpoint_valid:
            raise IOError("error writing to checkpoint file")

    @property
    def samples(self):
        """Returns raw nested samples
        """
        results = self._sampler.results
        samples = results.samples
        nest_samp = {}
        for i, param in enumerate(self.variable_params):
            nest_samp[param] = samples[:, i]
        nest_samp['logwt'] = results.logwt
        nest_samp['loglikelihood'] = results.logl
        return nest_samp

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
            fp.write_raw_samples(self.samples)

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


def sample_rwalk_mod(args):
    """ Modified version of dynesty.sampling.sample_rwalk

        Adapted from version used in bilby/dynesty
    """

    # Unzipping.
    (u, loglstar, axes, scale,
     prior_transform, loglikelihood, kwargs) = args
    rstate = numpy.random

    # Bounds
    nonbounded = kwargs.get('nonbounded', None)
    periodic = kwargs.get('periodic', None)
    reflective = kwargs.get('reflective', None)

    # Setup.
    n = len(u)
    walks = kwargs.get('walks', 10 * n)  # minimum number of steps
    maxmcmc = kwargs.get('maxmcmc', 2000)  # Maximum number of steps
    nact = kwargs.get('nact', 5)  # Number of ACT
    old_act = kwargs.get('old_act', walks)

    # Initialize internal variables
    accept = 0
    reject = 0
    nfail = 0
    act = numpy.inf
    u_list = []
    v_list = []
    logl_list = []

    ii = 0
    while ii < nact * act:
        ii += 1

        # Propose a direction on the unit n-sphere.
        drhat = rstate.randn(n)
        drhat /= numpy.linalg.norm(drhat)

        # Scale based on dimensionality.
        dr = drhat * rstate.rand() ** (1.0 / n)

        # Transform to proposal distribution.
        du = numpy.dot(axes, dr)
        u_prop = u + scale * du

        # Wrap periodic parameters
        if periodic is not None:
            u_prop[periodic] = numpy.mod(u_prop[periodic], 1)
        # Reflect
        if reflective is not None:
            u_prop[reflective] = reflect(u_prop[reflective])

        # Check unit cube constraints.
        if unitcheck(u_prop, nonbounded):
            pass
        else:
            nfail += 1
            # Only start appending to the chain once a single jump is made
            if accept > 0:
                u_list.append(u_list[-1])
                v_list.append(v_list[-1])
                logl_list.append(logl_list[-1])
            continue

        # Check proposed point.
        v_prop = prior_transform(numpy.array(u_prop))
        logl_prop = loglikelihood(numpy.array(v_prop))
        if logl_prop >= loglstar:
            u = u_prop
            v = v_prop
            logl = logl_prop
            accept += 1
            u_list.append(u)
            v_list.append(v)
            logl_list.append(logl)
        else:
            reject += 1
            # Only start appending to the chain once a single jump is made
            if accept > 0:
                u_list.append(u_list[-1])
                v_list.append(v_list[-1])
                logl_list.append(logl_list[-1])

        # If we've taken the minimum number of steps, calculate the ACT
        if accept + reject > walks:
            act = estimate_nmcmc(
                accept_ratio=accept / (accept + reject + nfail),
                old_act=old_act, maxmcmc=maxmcmc)

        # If we've taken too many likelihood evaluations then break
        if accept + reject > maxmcmc:
            logging.warning(
                "Hit maximum number of walks {} with accept={}, reject={}, "
                "and nfail={} try increasing maxmcmc"
                .format(maxmcmc, accept, reject, nfail))
            break

    # If the act is finite, pick randomly from within the chain
    if numpy.isfinite(act) and int(.5 * nact * act) < len(u_list):
        idx = numpy.random.randint(int(.5 * nact * act), len(u_list))
        u = u_list[idx]
        v = v_list[idx]
        logl = logl_list[idx]
    else:
        logging.debug("Unable to find a new point using walk: "
                      "returning a random point")
        u = numpy.random.uniform(size=n)
        v = prior_transform(u)
        logl = loglikelihood(v)

    blob = {'accept': accept, 'reject': reject, 'fail': nfail, 'scale': scale}
    kwargs["old_act"] = act

    ncall = accept + reject
    return u, v, logl, ncall, blob


def estimate_nmcmc(accept_ratio, old_act, maxmcmc, safety=5, tau=None):
    """ Estimate autocorrelation length of chain using acceptance fraction

    Using ACL = (2/acc) - 1 multiplied by a safety margin. Code adapated from
    CPNest:
        - https://github.com/johnveitch/cpnest/blob/master/cpnest/sampler.py
        - http://github.com/farr/Ensemble.jl

    Parameters
    ----------
    accept_ratio: float [0, 1]
        Ratio of the number of accepted points to the total number of points
    old_act: int
        The ACT of the last iteration
    maxmcmc: int
        The maximum length of the MCMC chain to use
    safety: int
        A safety factor applied in the calculation
    tau: int (optional)
        The ACT, if given, otherwise estimated.

    """
    if tau is None:
        tau = maxmcmc / safety

    if accept_ratio == 0.0:
        Nmcmc_exact = (1 + 1 / tau) * old_act
    else:
        Nmcmc_exact = (
            (1. - 1. / tau) * old_act +
            (safety / tau) * (2. / accept_ratio - 1.)
        )
        Nmcmc_exact = float(min(Nmcmc_exact, maxmcmc))
    return max(safety, int(Nmcmc_exact))

