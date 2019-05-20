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

from pycbc.inference.io import (MultinestFile, validate_checkpoint_files)
from pycbc.distributions import read_constraints_from_config
from pycbc.transforms import apply_transforms
from .base import BaseSampler
from .base_mcmc import get_optional_arg_from_config


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
    """
    name = "multinest"
    _io = MultinestFile

    def __init__(self, model, nlivepoints, checkpoint_interval=1000,
                 importance_nested_sampling=False,
                 evidence_tolerance=0.1, sampling_efficiency=0.01,
                 constraints=None):
        try:
            loglevel = logging.getLogger().getEffectiveLevel()
            logging.getLogger().setLevel(logging.WARNING)
            from pymultinest import Analyzer, run
            self.run_multinest = run
            self.analyzer = Analyzer
            logging.getLogger().setLevel(loglevel)
        except ImportError:
            raise ImportError("pymultinest is not installed.")

        super(MultinestSampler, self).__init__(model)

        self._constraints = constraints
        self._nlivepoints = nlivepoints
        self._ndim = len(model.variable_params)
        self._random_state = numpy.random.get_state()
        self._checkpoint_interval = checkpoint_interval
        self._ztol = evidence_tolerance
        self._eff = sampling_efficiency
        self._ins = importance_nested_sampling
        self._samples = None
        self._itercount = None
        self._logz = None
        self._dlogz = None
        self._importance_logz = None
        self._importance_dlogz = None

    @property
    def io(self):
        return self._io

    @property
    def niterations(self):
        """Get the current number of iterations.
        """
        itercount = self._itercount
        if itercount is None:
            itercount = 0
        return itercount

    @property
    def checkpoint_interval(self):
        """Get the number of iterations between checkpoints.
        """
        return self._checkpoint_interval

    @property
    def nlivepoints(self):
        """Get the number of live points used in sampling.
        """
        return self._nlivepoints

    @property
    def logz(self):
        """Get the current estimate of the log evidence.
        """
        return self._logz

    @property
    def dlogz(self):
        """Get the current error estimate of the log evidence.
        """
        return self._dlogz

    @property
    def importance_logz(self):
        """Get the current importance weighted estimate of the log
        evidence.
        """
        return self._importance_logz

    @property
    def importance_dlogz(self):
        """Get the current error estimate of the importance
        weighted log evidence.
        """
        return self._importance_dlogz

    @property
    def samples(self):
        """A dict mapping ``variable_params`` to arrays of samples currently
        in memory.
        """
        samples_dict = {p: self._samples[:, i] for i, p in
                        enumerate(self.model.variable_params)}
        return samples_dict

    @property
    def model_stats(self):
        """A dict mapping the model's ``default_stats`` to arrays of values.
        """
        stats = []
        for sample in self._samples:
            params = dict(zip(self.model.variable_params, sample))
            if self.model.sampling_transforms is not None:
                params = self.model.sampling_transforms.apply(params)
            self.model.update(**params)
            self.model.logposterior
            stats.append(self.model.get_current_stats())
        stats = numpy.array(stats)
        return {s: stats[:, i] for i, s in enumerate(self.model.default_stats)}

    def get_posterior_samples(self):
        """Read posterior samples from ASCII output file created by
        multinest.
        """
        post_file = self.backup_file[:-9]+'-post_equal_weights.dat'
        return numpy.loadtxt(post_file, ndmin=2)

    def check_if_finished(self):
        """Estimate remaining evidence to see if desired evidence-tolerance
        stopping criterion has been reached.
        """
        resume_file = self.backup_file[:-9] + '-resume.dat'
        current_vol, _, _ = numpy.loadtxt(
            resume_file, skiprows=6, unpack=True)
        maxloglike = max(self.get_posterior_samples()[:, -1])
        logz_remain = numpy.exp(maxloglike +
                                numpy.log(current_vol) - self.logz)
        logging.info("Estimate of remaining logZ is %s", logz_remain)
        done = logz_remain < self._ztol
        return done

    def set_initial_conditions(self, initial_distribution=None,
                               samples_file=None):
        """Sets the initial starting point for the sampler.

        If a starting samples file is provided, will also load the random
        state from it.
        """
        # use samples file to set the state of the sampler
        if samples_file is not None:
            self.set_state_from_file(samples_file)

    def set_state_from_file(self, filename):
        """Sets the state of the sampler back to the instance saved in a file.
        """
        with self.io(filename, 'r') as f_p:
            rstate = f_p.read_random_state()
        # set the numpy random state
        numpy.random.set_state(rstate)
        # set sampler's generator to the same state
        self._random_state = rstate

    def loglikelihood(self, cube, *extra_args):
        """Log likelihood evaluator that gets passed to multinest.
        """
        params = {p: v for p, v in zip(self.model.variable_params, cube)}
        # apply transforms
        if self.model.sampling_transforms is not None:
            params = self.model.sampling_transforms.apply(params)
        if self.model.waveform_transforms is not None:
            params = apply_transforms(params, self.model.waveform_transforms)
        # apply constraints
        if (self._constraints is not None and
                not all([c(params) for c in self._constraints])):
            return -numpy.inf
        self.model.update(**params)
        return self.model.loglikelihood

    def transform_prior(self, cube, *extra_args):
        """Transforms the unit hypercube that multinest makes its draws
        from, into the prior space defined in the config file.
        """
        prior_dists = self.model.prior_distribution.distributions
        dist_dict = {}
        for dist in prior_dists:
            dist_dict.update({param: dist for param in dist.params})
        for i, param in enumerate(self.model.variable_params):
            cube[i] = dist_dict[param].cdfinv(param, cube[i])
        return cube

    def run(self):
        """Runs the sampler until the specified evidence tolerance
        is reached.
        """
        if self.new_checkpoint:
            self._itercount = 0
        else:
            self.set_initial_conditions(samples_file=self.checkpoint_file)
            with self.io(self.checkpoint_file, "r") as f_p:
                self._itercount = f_p.niterations
        outputfiles_basename = self.backup_file[:-9] + '-'
        analyzer = self.analyzer(self._ndim,
                                 outputfiles_basename=outputfiles_basename)
        iterinterval = self.checkpoint_interval
        done = False
        while not done:
            logging.info("Running sampler for %s to %s iterations",
                         self.niterations, self.niterations + iterinterval)
            # run multinest
            self.run_multinest(self.loglikelihood, self.transform_prior,
                               self._ndim, n_live_points=self.nlivepoints,
                               evidence_tolerance=self._ztol,
                               sampling_efficiency=self._eff,
                               importance_nested_sampling=self._ins,
                               max_iter=iterinterval,
                               seed=numpy.random.randint(0, 1e6),
                               outputfiles_basename=outputfiles_basename,
                               multimodal=False, verbose=True)
            # parse results from multinest output files
            nest_stats = analyzer.get_mode_stats()
            self._logz = nest_stats["nested sampling global log-evidence"]
            self._dlogz = nest_stats[
                "nested sampling global log-evidence error"]
            if self._ins:
                self._importance_logz = nest_stats[
                    "nested importance sampling global log-evidence"]
                self._importance_dlogz = nest_stats[
                    "nested importance sampling global log-evidence error"]
            self._samples = self.get_posterior_samples()[:, :-1]
            logging.info("Have %s posterior samples", self._samples.shape[0])
            # update the itercounter
            self._itercount += iterinterval
            # make sure there's at least 1 posterior sample
            if self._samples.shape[0] == 0:
                continue
            # dump the current results
            self.checkpoint()
            # check if we're finished
            done = self.check_if_finished()

    def write_results(self, filename):
        """Writes samples, model stats, acceptance fraction, and random state
        to the given file.

        Parameters
        -----------
        filename : str
            The file to write to. The file is opened using the ``io`` class
            in an an append state.
        """
        with self.io(filename, 'a') as f_p:
            # write samples
            f_p.write_samples(self.samples, self.model.variable_params)
            # write stats
            f_p.write_samples(self.model_stats)
            # write evidence
            f_p.write_logevidence(self.logz, self.dlogz,
                                  self.importance_logz,
                                  self.importance_dlogz)
            # write random state (use default numpy.random_state)
            f_p.write_random_state()

    def checkpoint(self):
        """Dumps current samples to the checkpoint file."""
        logging.info("Writing samples to files")
        for f_n in [self.checkpoint_file, self.backup_file]:
            self.write_results(f_n)
            with self.io(f_n, "a") as f_p:
                f_p.write_niterations(self.niterations)
        logging.info("Validating checkpoint and backup files")
        checkpoint_valid = validate_checkpoint_files(
            self.checkpoint_file, self.backup_file)
        if not checkpoint_valid:
            raise IOError("error writing to checkpoint file")

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
        checkpoint = get_optional_arg_from_config(
            cp, section, 'checkpoint-interval', dtype=int)
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
        kwarg_names = ['evidence_tolerance', 'sampling_efficiency',
                       'importance_nested_sampling',
                       'checkpoint_interval']
        optional_kwargs = {k: v for k, v in
                           zip(kwarg_names, [ztol, eff, ins, checkpoint]) if
                           v is not None}
        obj = cls(model, nlivepoints, constraints=constraints,
                  **optional_kwargs)
        return obj
