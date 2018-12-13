# Copyright (C) 2018  Daniel Finstad
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
This modules provides classes and functions for using the Multinest sampler
packages for parameter estimation.
"""

from __future__ import absolute_import

import logging
import numpy
import pymultinest
from pycbc.pool import choose_pool

from .base import BaseSampler
from .base_mcmc import (BaseMCMC, raw_samples_to_dict,
                        blob_data_to_dict, get_optional_arg_from_config)
from pycbc.inference.io import (MultinestFile, validate_checkpoint_files)
from pycbc.distributions import read_constraints_from_config
from .. import models


#
# =============================================================================
#
#                                   Samplers
#
# =============================================================================
#

class MultinestSampler(BaseSampler):
    """This class is used to construct a nested sampler from
    the Multinest package.

    Parameters
    ----------
    model : model
        A model from ``pycbc.inference.models``.
    nlivepoints : int
        Number of live points to use in sampler.
    pool : function with map, Optional
        A provider of a map function that allows a function call to be run
        over multiple sets of arguments and possibly maps them to
        cores/nodes/etc.
    """
    name = "multinest"
    _io = MultinestFile

    def __init__(self, model, nlivepoints, checkpoint_interval=1000,
                 logpost_function=None, nprocesses=1, use_mpi=False,
                 importance_nested_sampling=False, evidence_tolerance=0.1,
                 sampling_efficiency=0.01, constraints=None):

        self.model = model
        # create a wrapper for calling the model
        if logpost_function is None:
            logpost_function = 'logposterior'
        model_call = models.CallModel(model, logpost_function)

        # Set up the pool
        if nprocesses > 1:
            # these are used to help paralleize over multiple cores / MPI
            models._global_instance = model_call
            model_call = models._call_global_model
        pool = choose_pool(mpi=use_mpi, processes=nprocesses)
        if pool is not None:
            pool.count = nprocesses

        # set up necessary attributes
        self._constraints = constraints
        self._nlivepoints = nlivepoints
        self._ndim = len(model.variable_params)
        rstate = numpy.random.get_state()
        self._random_state = rstate
        self._checkpoint_interval = checkpoint_interval
        self._ztol = evidence_tolerance
        self._eff = sampling_efficiency
        self._ins = importance_nested_sampling
        self._samples = None
        self._stats = None
        self._seed = 0 # FIXME need to pass in random seed?
        self._itercount = None
        self._lastclear = None
        self._logz = None
        self._dlogz = None
        self._importance_logz = None
        self._importance_dlogz = None

    @property
    def io(self):
        return self._io

    @property
    def niterations(self):
        """Get the current number of iterations."""
        itercount = self._itercount
        if itercount is None:
            itercount = 0
        lastclear = self._lastclear
        if lastclear is None:
            lastclear = 0
        return itercount + lastclear

    @property
    def checkpoint_interval(self):
        return self._checkpoint_interval

    @property
    def nlivepoints(self):
        return self._nlivepoints

    @property
    def logz(self):
        return self._logz

    @property
    def dlogz(self):
        return self._dlogz

    @property
    def importance_logz(self):
        return self._importance_logz

    @property
    def importance_dlogz(self):
        return self._importance_dlogz

    @property
    def samples(self):
        """A dict mapping ``variable_params`` to arrays of samples currently
        in memory.

        The arrays have shape ``nwalkers x niterations``.
        """
        samples_dict = {p: self._samples[:, i] for i, p in
                        enumerate(self.model.variable_params)}
        return samples_dict

    @property
    def model_stats(self):
        """A dict mapping the model's ``default_stats`` to arrays of values.

        The returned array has shape ``nwalkers x niterations``.
        """
        stats = self.model.default_stats
        callfuncs = [self.model._logjacobian, self.model._logprior,
                     self.model._loglikelihood]
        if 'loglr' in stats:
            callfuncs += [self.model._loglr]
        self._stats = {s: numpy.array([]) for s in stats}
        # calculate stats for each posterior sample
        for sample in self._samples:
            params = dict(zip(self.model.variable_params, sample))
            if self.model.sampling_transforms is not None:
                params = self.model.sampling_transforms.apply(params)
            self.model.update(**params)
            self.model.logposterior # this is necessary to avoid logprior being nan??
            for c in callfuncs:
                c()
            current_stats = self.model.get_current_stats(names=stats)
            for i, s in enumerate(stats):
                self._stats[s] = numpy.append(self._stats[s], current_stats[i])
        return self._stats

    def acceptance_fraction(self):
        # FIXME maybe get the real acceptance fraction at some point?
        return numpy.zeros(self.nlivepoints)

    def get_posterior_samples(self):
        # this is multinest's equal weighted posterior file
        post_file = self.backup_file[:-9]+'-post_equal_weights.dat'
        return numpy.loadtxt(post_file, ndmin=2)

    def check_if_finished(self):
        # do calculation that multinest does
        resume_file = self.backup_file[:-9] + '-resume.dat'
        current_vol, _, _ = numpy.loadtxt(
            resume_file, skiprows=6, unpack=True)
        maxloglike = max(self.get_posterior_samples()[:, -1])
        logz_remain = numpy.exp(maxloglike + numpy.log(current_vol) - self.logz)
        logging.info("Estimate of remaining logZ is {}".format(logz_remain))
        done = logz_remain < self._ztol
        return done

    def clear_samples(self):
        """Clears the samples and stats from memory.
        """
        # store the iteration that the clear is occuring on
        self._lastclear = self.niterations
        self._itercounter = 0
        # now clear the chain
        self._sampler.reset()
        self._sampler.clear_blobs()

    def set_initial_conditions(self, initial_distribution=None,
                               samples_file=None):
        """Sets the initial starting point for the MCMC.

        If a starting samples file is provided, will also load the random
        state from it.
        """
        #self.set_p0(samples_file=samples_file, prior=initial_distribution)
        # if a samples file was provided, use it to set the state of the
        # sampler
        if samples_file is not None:
            self.set_state_from_file(samples_file)

    def set_state_from_file(self, filename):
        """Sets the state of the sampler back to the instance saved in a file.
        """
        with self.io(filename, 'r') as fp:
            rstate = fp.read_random_state()
        # set the numpy random state
        numpy.random.set_state(rstate)
        # set emcee's generator to the same state
        self._random_state = rstate

    def ns_loglikelihood(self, cube):
        params = {p: v for p, v in zip(self.model.variable_params, cube)}
        # apply transforms
        if self.model.sampling_transforms is not None:
            params = self.model.sampling_transforms.apply(params)
        params = self.model._transform_params(**params) # waveform transforms
        # apply constraints
        if self._constraints is not None and not all([c(params) for c in self._constraints]):
            return -1e90
        else:
            # update model with current params
            self.model.update(**params)
            return self.model.logposterior

    def ns_prior(self, cube):
        transformed_cube = numpy.array(cube).copy()
        prior_dists = self.model.prior_distribution.distributions
        dist_dict = {d.params[0]: d for d in prior_dists}
        for i, p in enumerate(self.model.variable_params):
            bounds = dist_dict[p].bounds
            if dist_dict[p].name in ['uniform', 'uniform_angle']:
                scale = bounds[p].max - bounds[p].min
                transformed_cube[i] = cube[i] * scale + bounds[p].min
            else:
                transformed_cube[i] = dist_dict[p]._cdfinv(p, cube[i])
        return transformed_cube

    def run(self):
        if self.new_checkpoint:
            self._lastclear = 0
        else:
            with self.io(self.checkpoint_file, "r") as fp:
                self._lastclear = fp.niterations
        outputfiles_basename = self.backup_file[:-9] + '-'
        a = pymultinest.Analyzer(self._ndim,
                                 outputfiles_basename=outputfiles_basename)
        self._itercount = 0
        iterinterval = self.checkpoint_interval
        done = False
        while not done:
            logging.info("Running sampler for {} to {} iterations".format(
                self.niterations, self.niterations + iterinterval))
            # run multinest
            self.solve(self.ns_loglikelihood, self.ns_prior, self._ndim,
                       n_live_points=self.nlivepoints,
                       evidence_tolerance=self._ztol,
                       sampling_efficiency=self._eff,
                       importance_nested_sampling=self._ins,
                       max_iter=iterinterval,
                       outputfiles_basename=outputfiles_basename,
                       multimodal=False, verbose=True)
            # parse results from multinest output files
            nest_stats = a.get_mode_stats()
            self._logz = nest_stats['nested sampling global log-evidence']
            self._dlogz = nest_stats['nested sampling global log-evidence error']
            if self._ins:
                self._importance_logz = nest_stats['nested importance sampling global log-evidence']
                self._importance_dlogz = nest_stats['nested importance sampling global log-evidence error']
            self._samples = self.get_posterior_samples()[:,:-1]
            logging.info("Have {} posterior samples".format(self._samples.shape[0]))
            # update the itercounter
            self._itercount = self._itercount + iterinterval
            # dump the current results
            self.checkpoint()
            # check if we're finished
            done = self.check_if_finished()

    def solve(self, LogLikelihood, Prior, n_dims, **kwargs):
	kwargs['n_dims'] = n_dims
	files_temporary = False
	if 'outputfiles_basename' not in kwargs:
		files_temporary = True
		tempdir = tempfile.mkdtemp('pymultinest')
		kwargs['outputfiles_basename'] = tempdir + '/'
	outputfiles_basename = kwargs['outputfiles_basename']
	def SafePrior(cube, ndim, nparams):
	    try:
	        a = numpy.array([cube[i] for i in range(n_dims)])
		b = Prior(a)
		for i in range(n_dims):
		    cube[i] = b[i]
	    except Exception as e:
		import sys
		sys.stderr.write('ERROR in prior: %s\n' % e)
	        sys.exit(1)

	def SafeLoglikelihood(cube, ndim, nparams, lnew):
	    try:
	        a = numpy.array([cube[i] for i in range(n_dims)])
	        l = float(LogLikelihood(a))
		if not numpy.isfinite(l):
		    import sys
		    sys.stderr.write('WARNING: loglikelihood not finite: %f\n' % (l))
		    sys.stderr.write('         for parameters: %s\n' % a)
		    sys.stderr.write('         returned very low value instead\n')
		    return -1e100
		return l
	    except Exception as e:
		import sys
		sys.stderr.write('ERROR in loglikelihood: %s\n' % e)
		sys.exit(1)

	kwargs['LogLikelihood'] = SafeLoglikelihood
	kwargs['Prior'] = SafePrior
	pymultinest.run(**kwargs)

    def write_results(self, filename):
        """Writes samples, model stats, acceptance fraction, and random state
        to the given file.

        Parameters
        -----------
        filename : str
            The file to write to. The file is opened using the ``io`` class
            in an an append state.
        """
        # check to make sure there's at least 1 posterior sample
        if self._samples.shape[0] == 0:
            return
        with self.io(filename, 'a') as fp:
            # write samples
            fp.write_samples(self.samples, self.model.variable_params)
            # write stats
            fp.write_samples(self.model_stats)
            # write evidence
            fp.write_logevidence(self.logz, self.dlogz,
                                 self.importance_logz,
                                 self.importance_dlogz)
            # write acceptance
            fp.write_acceptance_fraction(self.acceptance_fraction())
            # write random state (use default numpy.random_state)
            fp.write_random_state()

    def checkpoint(self):
        """Dumps current samples to the checkpoint file."""
        logging.info("Writing samples to files")
        for fn in [self.checkpoint_file, self.backup_file]:
            self.write_results(fn)
            with self.io(fn, "a") as fp:
                fp.write_niterations(self.niterations)
        logging.info("Validating checkpoint and backup files")
        checkpoint_valid = validate_checkpoint_files(
            self.checkpoint_file, self.backup_file)
        if not checkpoint_valid:
            raise IOError("error writing to checkpoint file")
        # FIXME clear memory?

    def finalize(self):
        """All data is written by the last checkpoint in the run method, so
        this just passes."""
        pass

    @classmethod
    def from_config(cls, cp, model, nprocesses=1, use_mpi=False):
        """Loads the sampler from the given config file."""
        section = "sampler"
        # check name
        assert cp.get(section, "name") == cls.name, (
            "name in section [sampler] must match mine")
        # get the number of live points to use
        nlivepoints = int(cp.get(section, "nlivepoints"))
        # get the checkpoint interval, if it's specified
        checkpoint_interval = get_optional_arg_from_config(
            cp, section, 'checkpoint-interval', dtype=int)
        # get the logpost function
        lnpost = get_optional_arg_from_config(cp, section, 'logpost-function')
        # get the evidence tolerance, if specified
        ztol = get_optional_arg_from_config(cp, section, 'evidence-tolerance',
                                            dtype=float)
        # get the sampling efficiency, if specified
        eff = get_optional_arg_from_config(cp, section, 'sampling-efficiency',
                                           dtype=float)
        # get importance nested sampling setting, if specified
        ins = get_optional_arg_from_config(cp, section,
                                           'importance-nested-sampling',
                                           dtype=bool)
        # get constraints since we can't use the joint prior distribution
        constraints = read_constraints_from_config(cp)
        # build optional kwarg dict
        kwargnames = ['evidence_tolerance', 'sampling_efficiency',
                      'importance_nested_sampling']
        optional_kwargs = {k: v for k, v in zip(kwargnames, [ztol, eff, ins]) if
                           v is not None}
        obj = cls(model, nlivepoints, checkpoint_interval=checkpoint_interval,
                  logpost_function=lnpost, nprocesses=nprocesses,
                  use_mpi=use_mpi, constraints=constraints, **optional_kwargs)
        return obj
