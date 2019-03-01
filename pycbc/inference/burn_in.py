# Copyright (C) 2017  Collin Capano
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
This modules provides classes and functions for determining when Markov Chains
have burned in.
"""

from __future__ import division

import numpy
from scipy.stats import ks_2samp

from pycbc.io.record import get_vars_from_arg

# The value to use for a burn-in iteration if a chain is not burned in
NOT_BURNED_IN_ITER = -1


#
# =============================================================================
#
#                              Convenience functions
#
# =============================================================================
#


def ks_test(samples1, samples2, threshold=0.9):
    """Applies a KS test to determine if two sets of samples are the same.

    The ks test is applied parameter-by-parameter. If the two-tailed p-value
    returned by the test is greater than ``threshold``, the samples are
    considered to be the same.

    Parameters
    ----------
    samples1 : dict
        Dictionary of mapping parameters to the first set of samples.
    samples2 : dict
        Dictionary of mapping parameters to the second set of samples.
    threshold : float
        The thershold to use for the p-value. Default is 0.9.

    Returns
    -------
    dict :
        Dictionary mapping parameter names to booleans indicating whether the
        given parameter passes the KS test.
    """
    is_the_same = {}
    assert set(samples1.keys()) == set(samples2.keys()), (
        "samples1 and 2 must have the same parameters")
    # iterate over the parameters
    for param in samples1:
        s1 = samples1[param]
        s2 = samples2[param]
        _, p_value = ks_2samp(s1, s2)
        is_the_same[param] = p_value > threshold
    return is_the_same


def max_posterior(lnps_per_walker, dim):
    """Burn in based on samples being within dim/2 of maximum posterior.

    Parameters
    ----------
    lnps_per_walker : 2D array
        Array of values that are proportional to the log posterior values. Must
        have shape ``nwalkers x niterations``.
    dim : int
        The dimension of the parameter space.

    Returns
    -------
    burn_in_idx : array of int
        The burn in indices of each walker. If a walker is not burned in, its
        index will be be equal to the length of the chain.
    is_burned_in : array of bool
        Whether or not a walker is burned in.
    """
    if len(lnps_per_walker.shape) != 2:
        raise ValueError("lnps_per_walker must have shape "
                         "nwalkers x niterations")
    # find the value to compare against
    max_p = lnps_per_walker.max()
    criteria = max_p - dim/2.
    nwalkers, _ = lnps_per_walker.shape
    burn_in_idx = numpy.empty(nwalkers, dtype=int)
    is_burned_in = numpy.empty(nwalkers, dtype=bool)
    # find the first iteration in each chain where the logpost has exceeded
    # max_p - dim/2
    for ii in range(nwalkers):
        chain = lnps_per_walker[ii, :]
        passedidx = numpy.where(chain >= criteria)[0]
        is_burned_in[ii] = passedidx.size > 0
        if is_burned_in[ii]:
            burn_in_idx[ii] = passedidx[0]
        else:
            burn_in_idx[ii] = NOT_BURNED_IN_ITER
    return burn_in_idx, is_burned_in


def posterior_step(logposts, dim):
    """Finds the last time a chain made a jump > dim/2.

    Parameters
    ----------
    logposts : array
        1D array of values that are proportional to the log posterior values.
    dim : int
        The dimension of the parameter space.

    Returns
    -------
    int
        The index of the last time the logpost made a jump > dim/2. If that
        never happened, returns 0.
    """
    if logposts.ndim > 1:
        raise ValueError("logposts must be a 1D array")
    criteria = dim/2.
    dp = numpy.diff(logposts)
    indices = numpy.where(dp >= criteria)[0]
    if indices.size > 0:
        idx = indices[-1] + 1
    else:
        idx = 0
    return idx


#
# =============================================================================
#
#                              Burn in classes
#
# =============================================================================
#


class MCMCBurnInTests(object):
    """Provides methods for estimating burn-in of an ensemble MCMC."""

    available_tests = ('halfchain', 'min_iterations', 'max_posterior',
                       'posterior_step', 'nacl', 'ks_test',
                       )

    def __init__(self, sampler, burn_in_test, **kwargs):
        self.sampler = sampler
        # determine the burn-in tests that are going to be done
        self.do_tests = get_vars_from_arg(burn_in_test)
        self.burn_in_test = burn_in_test
        self.burn_in_data = {t: {} for t in self.do_tests}
        self.is_burned_in = False
        self.burn_in_iteration = NOT_BURNED_IN_ITER
        self.burn_in_index = NOT_BURNED_IN_ITER
        # Arguments specific to each test...
        # for nacl:
        self._nacls = int(kwargs.pop('nacls', 5))
        # for kstest:
        self._ksthreshold = float(kwargs.pop('ks_threshold', 0.9))
        # for max_posterior and posterior_step
        self._ndim = int(kwargs.pop('ndim', len(sampler.variable_params)))
        # for min iterations
        self._min_iterations = int(kwargs.pop('min_iterations', 0))

    def _getniters(self, filename):
        """Convenience function to get the number of iterations in the file.

        If `niterations` hasn't been written to the file yet, just returns 0.
        """
        with self.sampler.io(filename, 'r') as fp:
            try:
                niters = fp.niterations
            except KeyError:
                niters = 0
        return niters

    def _getnsamples(self, filename):
        """Convenience function to get the number of samples saved in the file.

        If no samples have been written to the file yet, just returns 0.
        """
        with self.sampler.io(filename, 'r') as fp:
            try:
                group = fp[fp.samples_group]
                # we'll just use the first parameter
                params = group.keys()
                nsamples = group[params[0]].shape[-1]
            except (KeyError, IndexError):
                nsamples = 0
        return nsamples

    def _index2iter(self, filename, index):
        """Converts the index in some samples at which burn in occurs to the
        iteration of the sampler that corresponds to.
        """
        with self.sampler.io(filename, 'r') as fp:
            thin_interval = fp.thinned_by
        return index * thin_interval

    def _iter2index(self, filename, iteration):
        """Converts an iteration to the index it corresponds to.
        """
        with self.sampler.io(filename, 'r') as fp:
            thin_interval = fp.thinned_by
        return iteration // thin_interval

    def _getlogposts(self, filename):
        """Convenience function for retrieving log posteriors.

        Parameters
        ----------
        filename : str
            The file to read.

        Returns
        -------
        array
            The log posterior values. They are not flattened, so have dimension
            nwalkers x niterations.
        """
        with self.sampler.io(filename, 'r') as fp:
            samples = fp.read_raw_samples(
                ['loglikelihood', 'logprior'], thin_start=0, thin_interval=1,
                flatten=False)
            logposts = samples['loglikelihood'] + samples['logprior']
        return logposts

    def _getacls(self, filename, start_index):
        """Convenience function for calculating acls for the given filename.

        Since we calculate the acls, this will also store it to the sampler.
        """
        acls = self.sampler.compute_acl(filename, start_index=start_index)
        # since we calculated it, save the acls to the sampler...
        # but only do this if this is the only burn in test
        if len(self.do_tests) == 1:
            self.sampler.acls = acls
        return acls

    def halfchain(self, filename):
        """Just uses half the chain as the burn-in iteration.
        """
        niters = self._getniters(filename)
        data = self.burn_in_data['halfchain']
        # this test cannot determine when something will burn in
        # only when it was not burned in in the past
        data['is_burned_in'] = True
        data['burn_in_iteration'] = niters/2

    def min_iterations(self, filename):
        """Just checks that the sampler has been run for the minimum number
        of iterations.
        """
        niters = self._getniters(filename)
        data = self.burn_in_data['min_iterations']
        data['is_burned_in'] = self._min_iterations < niters
        if data['is_burned_in']:
            data['burn_in_iteration'] = self._min_iterations
        else:
            data['burn_in_iteration'] = NOT_BURNED_IN_ITER

    def max_posterior(self, filename):
        """Applies max posterior test to self."""
        logposts = self._getlogposts(filename)
        burn_in_idx, is_burned_in = max_posterior(logposts, self._ndim)
        data = self.burn_in_data['max_posterior']
        # required things to store
        data['is_burned_in'] = is_burned_in.all()
        if data['is_burned_in']:
            data['burn_in_iteration'] = self._index2iter(
                filename, burn_in_idx.max())
        else:
            data['burn_in_iteration'] = NOT_BURNED_IN_ITER
        # additional info
        data['iteration_per_walker'] = self._index2iter(filename, burn_in_idx)
        data['status_per_walker'] = is_burned_in

    def posterior_step(self, filename):
        """Applies the posterior-step test."""
        logposts = self._getlogposts(filename)
        burn_in_idx = numpy.array([posterior_step(logps, self._ndim)
                                   for logps in logposts])
        data = self.burn_in_data['posterior_step']
        # this test cannot determine when something will burn in
        # only when it was not burned in in the past
        data['is_burned_in'] = True
        data['burn_in_iteration'] = self._index2iter(
            filename, burn_in_idx.max())
        # additional info
        data['iteration_per_walker'] = self._index2iter(filename, burn_in_idx)

    def nacl(self, filename):
        """Burn in based on ACL.

        This applies the following test to determine burn in:

        1. The first half of the chain is ignored.

        2. An ACL is calculated from the second half.

        3. If ``nacls`` times the ACL is < the length of the chain / 2,
           the chain is considered to be burned in at the half-way point.
        """
        nsamples = self._getnsamples(filename)
        kstart = int(nsamples / 2.)
        acls = self._getacls(filename, start_index=kstart)
        is_burned_in = {param: (self._nacls * acl) < kstart
                        for (param, acl) in acls.items()}
        data = self.burn_in_data['nacl']
        # required things to store
        data['is_burned_in'] = all(is_burned_in.values())
        if data['is_burned_in']:
            data['burn_in_iteration'] = self._index2iter(filename, kstart)
        else:
            data['burn_in_iteration'] = NOT_BURNED_IN_ITER
        # additional information
        data['status_per_parameter'] = is_burned_in

    def ks_test(self, filename):
        """Applies ks burn-in test."""
        nsamples = self._getnsamples(filename)
        with self.sampler.io(filename, 'r') as fp:
            # get the samples from the mid point
            samples1 = fp.read_raw_samples(
                ['loglikelihood', 'logprior'], iteration=int(nsamples/2.))
            # get the last samples
            samples2 = fp.read_raw_samples(
                ['loglikelihood', 'logprior'], iteration=-1)
        # do the test
        # is_the_same is a dictionary of params --> bool indicating whether or
        # not the 1D marginal is the same at the half way point
        is_the_same = ks_test(samples1, samples2, threshold=self._ksthreshold)
        data = self.burn_in_data['ks_test']
        # required things to store
        data['is_burned_in'] = all(is_the_same.values())
        if data['is_burned_in']:
            data['burn_in_iteration'] = self._index2iter(
                filename, int(nsamples/2.))
        else:
            data['burn_in_iteration'] = NOT_BURNED_IN_ITER
        # additional
        data['status_per_parameter'] = is_the_same

    def evaluate(self, filename):
        """Runs all of the burn-in tests."""
        for tst in self.do_tests:
            getattr(self, tst)(filename)
        # The iteration to use for burn-in depends on the logic in the burn-in
        # test string. For example, if the test was 'max_posterior | nacl' and
        # max_posterior burned-in at iteration 5000 while nacl burned in at
        # iteration 6000, we'd want to use 5000 as the burn-in iteration.
        # However, if the test was 'max_posterior & nacl', we'd want to use
        # 6000 as the burn-in iteration. The code below handles all cases by
        # doing the following: first, take the collection of burn in iterations
        # from all the burn in tests that were applied.  Next, cycle over the
        # iterations in increasing order, checking which tests have burned in
        # by that point. Then evaluate the burn-in string at that point to see
        # if it passes, and if so, what the iteration is. The first point that
        # the test passes is used as the burn-in iteration.
        data = self.burn_in_data
        burn_in_iters = numpy.unique([data[t]['burn_in_iteration']
                                      for t in self.do_tests])
        burn_in_iters.sort()
        for ii in burn_in_iters:
            test_results = {t: (data[t]['is_burned_in'] &
                                0 <= data[t]['burn_in_iteration'] <= ii)
                            for t in self.do_tests}
            is_burned_in = eval(self.burn_in_test, {"__builtins__": None},
                                test_results)
            if is_burned_in:
                break
        self.is_burned_in = is_burned_in
        if is_burned_in:
            self.burn_in_iteration = ii
            self.burn_in_index = self._iter2index(filename, ii)
        else:
            self.burn_in_iteration = NOT_BURNED_IN_ITER
            self.burn_in_index = NOT_BURNED_IN_ITER

    @classmethod
    def from_config(cls, cp, sampler):
        """Loads burn in from section [sampler-burn_in]."""
        section = 'sampler'
        tag = 'burn_in'
        burn_in_test = cp.get_opt_tag(section, 'burn-in-test', tag)
        kwargs = {}
        if cp.has_option_tag(section, 'nacl', tag):
            kwargs['nacl'] = int(cp.get_opt_tag(section, 'nacl', tag))
        if cp.has_option_tag(section, 'ks-threshold', tag):
            kwargs['ks_threshold'] = float(
                cp.get_opt_tag(section, 'ks-threshold', tag))
        if cp.has_option_tag(section, 'ndim', tag):
            kwargs['ndim'] = int(
                cp.get_opt_tag(section, 'ndim', tag))
        if cp.has_option_tag(section, 'min-iterations', tag):
            kwargs['min_iterations'] = int(
                cp.get_opt_tag(section, 'min-iterations', tag))
        return cls(sampler, burn_in_test, **kwargs)


class MultiTemperedMCMCBurnInTests(MCMCBurnInTests):
    """Adds support for multiple temperatures to the MCMCBurnInTests."""

    def _getacls(self, filename, start_index):
        """Convenience function for calculating acls for the given filename.

        This function is used by the ``n_acl`` burn-in test. That function
        expects the returned ``acls`` dict to just report a single ACL for
        each parameter. Since multi-tempered samplers return an array of ACLs
        for each parameter instead, this takes the max over the array before
        returning.

        Since we calculate the acls, this will also store it to the sampler.
        """
        acls = super(MultiTemperedMCMCBurnInTests, self)._getacls(
            filename, start_index)
        # return the max for each parameter
        return {param: vals.max() for (param, vals) in acls.items()}

    def _getlogposts(self, filename):
        """Convenience function for retrieving log posteriors.

        This just gets the coldest temperature chain, and returns arrays with
        shape nwalkers x niterations, so the parent class can run the same
        ``posterior_step`` function.
        """
        with self.sampler.io(filename, 'r') as fp:
            samples = fp.read_raw_samples(
                ['loglikelihood', 'logprior'], thin_start=0, thin_interval=1,
                temps=0, flatten=False)
            # reshape to drop the first dimension
            for (stat, arr) in samples.items():
                _, nwalkers, niterations = arr.shape
                samples[stat] = arr.reshape((nwalkers, niterations))
            logposts = samples['loglikelihood'] + samples['logprior']
        return logposts
