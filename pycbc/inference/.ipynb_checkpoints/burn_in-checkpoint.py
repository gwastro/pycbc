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


import logging
from abc import ABCMeta, abstractmethod
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


def nacl(nsamples, acls, nacls=5):
    """Burn in based on ACL.

    This applies the following test to determine burn in:

    1. The first half of the chain is ignored.

    2. An ACL is calculated from the second half.

    3. If ``nacls`` times the ACL is < the length of the chain / 2,
       the chain is considered to be burned in at the half-way point.

    Parameters
    ----------
    nsamples : int
        The number of samples of in the chain(s).
    acls : dict
        Dictionary of parameter -> ACL(s). The ACLs for each parameter may
        be an integer or an array of integers (for multiple chains).
    nacls : int, optional
        The number of ACLs the chain(s) must have gone past the halfway point
        in order to be considered burned in. Default is 5.

    Returns
    -------
    dict
        Dictionary of parameter -> boolean(s) indicating if the chain(s) pass
        the test. If an array of values was provided for the acls, the values
        will be arrays of booleans.
    """
    kstart = int(nsamples / 2.)
    return {param: (nacls * acl) < kstart for (param, acl) in acls.items()}


def evaluate_tests(burn_in_test, test_is_burned_in, test_burn_in_iter):
    """Evaluates burn in data from multiple tests.

    The iteration to use for burn-in depends on the logic in the burn-in
    test string. For example, if the test was 'max_posterior | nacl' and
    max_posterior burned-in at iteration 5000 while nacl burned in at
    iteration 6000, we'd want to use 5000 as the burn-in iteration.
    However, if the test was 'max_posterior & nacl', we'd want to use
    6000 as the burn-in iteration. This function handles all cases by
    doing the following: first, take the collection of burn in iterations
    from all the burn in tests that were applied.  Next, cycle over the
    iterations in increasing order, checking which tests have burned in
    by that point. Then evaluate the burn-in string at that point to see
    if it passes, and if so, what the iteration is. The first point that
    the test passes is used as the burn-in iteration.

    Parameters
    ----------
    burn_in_test : str
        The test to apply; e.g., ``'max_posterior & nacl'``.
    test_is_burned_in : dict
        Dictionary of test name -> boolean indicating whether a specific burn
        in test has passed.
    test_burn_in_iter : dict
        Dictionary of test name -> int indicating when a specific test burned
        in.

    Returns
    -------
    is_burned_in : bool
        Whether or not the data passes all burn in tests.
    burn_in_iteration :
        The iteration at which all the tests pass. If the tests did not all
        pass (``is_burned_in`` is false), then returns
        :py:data:`NOT_BURNED_IN_ITER`.
    """
    burn_in_iters = numpy.unique(list(test_burn_in_iter.values()))
    burn_in_iters.sort()
    for ii in burn_in_iters:
        test_results = {t: (test_is_burned_in[t] &
                            0 <= test_burn_in_iter[t] <= ii)
                        for t in test_is_burned_in}
        is_burned_in = eval(burn_in_test, {"__builtins__": None},
                            test_results)
        if is_burned_in:
            break
    if not is_burned_in:
        ii = NOT_BURNED_IN_ITER
    return is_burned_in, ii


#
# =============================================================================
#
#                              Burn in classes
#
# =============================================================================
#


class BaseBurnInTests(metaclass=ABCMeta):
    """Base class for burn in tests."""

    available_tests = ('halfchain', 'min_iterations', 'max_posterior',
                       'posterior_step', 'nacl',
                       )

    # pylint: disable=unnecessary-pass

    def __init__(self, sampler, burn_in_test, **kwargs):
        self.sampler = sampler
        # determine the burn-in tests that are going to be done
        self.do_tests = get_vars_from_arg(burn_in_test)
        self.burn_in_test = burn_in_test
        self.is_burned_in = False
        self.burn_in_iteration = NOT_BURNED_IN_ITER
        self.test_is_burned_in = {}  # burn in status per test
        self.test_burn_in_iteration = {}  # burn in iter per test
        self.test_aux_info = {}  # any additional information the test stores
        # Arguments specific to each test...
        # for nacl:
        self._nacls = int(kwargs.pop('nacls', 5))
        # for max_posterior and posterior_step
        self._ndim = int(kwargs.pop('ndim', len(sampler.variable_params)))
        # for min iterations
        self._min_iterations = int(kwargs.pop('min_iterations', 0))

    @abstractmethod
    def burn_in_index(self, filename):
        """The burn in index (retrieved from the iteration).

        This is an abstract method because how this is evaluated depends on
        if this is an ensemble MCMC or not.
        """
        pass

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
                params = list(group.keys())
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
        """
        return self.sampler.compute_acl(filename, start_index=start_index)

    def _getaux(self, test):
        """Convenience function for getting auxilary information.

        Parameters
        ----------
        test : str
            The name of the test to retrieve auxilary information about.

        Returns
        -------
        dict
            The ``test_aux_info[test]`` dictionary. If a dictionary does
            not exist yet for the given test, an empty dictionary will be
            created and saved to ``test_aux_info[test]``.
        """
        try:
            aux = self.test_aux_info[test]
        except KeyError:
            aux = self.test_aux_info[test] = {}
        return aux

    def halfchain(self, filename):
        """Just uses half the chain as the burn-in iteration.
        """
        niters = self._getniters(filename)
        # this test cannot determine when something will burn in
        # only when it was not burned in in the past
        self.test_is_burned_in['halfchain'] = True
        self.test_burn_in_iteration['halfchain'] = niters//2

    def min_iterations(self, filename):
        """Just checks that the sampler has been run for the minimum number
        of iterations.
        """
        niters = self._getniters(filename)
        is_burned_in = self._min_iterations < niters
        if is_burned_in:
            burn_in_iter = self._min_iterations
        else:
            burn_in_iter = NOT_BURNED_IN_ITER
        self.test_is_burned_in['min_iterations'] = is_burned_in
        self.test_burn_in_iteration['min_iterations'] = burn_in_iter

    @abstractmethod
    def max_posterior(self, filename):
        """Carries out the max posterior test and stores the results."""
        pass

    @abstractmethod
    def posterior_step(self, filename):
        """Carries out the posterior step test and stores the results."""
        pass

    @abstractmethod
    def nacl(self, filename):
        """Carries out the nacl test and stores the results."""
        pass

    @abstractmethod
    def evaluate(self, filename):
        """Performs all tests and evaluates the results to determine if and
        when all tests pass.
        """
        pass

    def write(self, fp, path=None):
        """Writes burn-in info to an open HDF file.

        Parameters
        ----------
        fp : pycbc.inference.io.base.BaseInferenceFile
            Open HDF file to write the data to. The HDF file should be an
            instance of a pycbc BaseInferenceFile.
        path : str, optional
            Path in the HDF file to write the data to. Default is (None) is
            to write to the path given by the file's ``sampler_group``
            attribute.
        """
        if path is None:
            path = fp.sampler_group
        fp.write_data('burn_in_test', self.burn_in_test, path)
        fp.write_data('is_burned_in', self.is_burned_in, path)
        fp.write_data('burn_in_iteration', self.burn_in_iteration, path)
        testgroup = 'burn_in_tests'
        # write individual test data
        for tst in self.do_tests:
            subpath = '/'.join([path, testgroup, tst])
            fp.write_data('is_burned_in', self.test_is_burned_in[tst], subpath)
            fp.write_data('burn_in_iteration',
                          self.test_burn_in_iteration[tst],
                          subpath)
            # write auxiliary info
            if tst in self.test_aux_info:
                for name, data in self.test_aux_info[tst].items():
                    fp.write_data(name, data, subpath)

    @staticmethod
    def _extra_tests_from_config(cp, section, tag):
        """For loading class-specific tests."""
        # pylint: disable=unused-argument
        return {}

    @classmethod
    def from_config(cls, cp, sampler):
        """Loads burn in from section [sampler-burn_in]."""
        section = 'sampler'
        tag = 'burn_in'
        burn_in_test = cp.get_opt_tag(section, 'burn-in-test', tag)
        kwargs = {}
        if cp.has_option_tag(section, 'nacl', tag):
            kwargs['nacl'] = int(cp.get_opt_tag(section, 'nacl', tag))
        if cp.has_option_tag(section, 'ndim', tag):
            kwargs['ndim'] = int(
                cp.get_opt_tag(section, 'ndim', tag))
        if cp.has_option_tag(section, 'min-iterations', tag):
            kwargs['min_iterations'] = int(
                cp.get_opt_tag(section, 'min-iterations', tag))
        # load any class specific tests
        kwargs.update(cls._extra_tests_from_config(cp, section, tag))
        return cls(sampler, burn_in_test, **kwargs)


class MCMCBurnInTests(BaseBurnInTests):
    """Burn-in tests for collections of independent MCMC chains.

    This differs from EnsembleMCMCBurnInTests in that chains are treated as
    being independent of each other. The ``is_burned_in`` attribute will be
    True if `any` chain passes the burn in tests (whereas in MCMCBurnInTests,
    all chains must pass the burn in tests). In other words, independent
    samples can be collected even if all of the chains are not burned in.
    """
    def __init__(self, sampler, burn_in_test, **kwargs):
        super(MCMCBurnInTests, self).__init__(sampler, burn_in_test, **kwargs)
        try:
            nchains = sampler.nchains
        except AttributeError:
            nchains = sampler.nwalkers
        self.nchains = nchains
        self.is_burned_in = numpy.zeros(self.nchains, dtype=bool)
        self.burn_in_iteration = numpy.repeat(NOT_BURNED_IN_ITER, self.nchains)

    def burn_in_index(self, filename):
        """The burn in index (retrieved from the iteration)."""
        burn_in_index = self._iter2index(filename, self.burn_in_iteration)
        # don't set if it isn't burned in
        burn_in_index[~self.is_burned_in] = NOT_BURNED_IN_ITER
        return burn_in_index

    def max_posterior(self, filename):
        """Applies max posterior test."""
        logposts = self._getlogposts(filename)
        burn_in_idx, is_burned_in = max_posterior(logposts, self._ndim)
        # convert index to iterations
        burn_in_iter = self._index2iter(filename, burn_in_idx)
        burn_in_iter[~is_burned_in] = NOT_BURNED_IN_ITER
        # save
        test = 'max_posterior'
        self.test_is_burned_in[test] = is_burned_in
        self.test_burn_in_iteration[test] = burn_in_iter

    def posterior_step(self, filename):
        """Applies the posterior-step test."""
        logposts = self._getlogposts(filename)
        burn_in_idx = numpy.array([posterior_step(logps, self._ndim)
                                   for logps in logposts])
        # this test cannot determine when something will burn in
        # only when it was not burned in in the past
        test = 'posterior_step'
        if test not in self.test_is_burned_in:
            self.test_is_burned_in[test] = numpy.ones(self.nchains, dtype=bool)
        # convert index to iterations
        self.test_burn_in_iteration[test] = self._index2iter(filename,
                                                             burn_in_idx)

    def nacl(self, filename):
        """Applies the :py:func:`nacl` test."""
        nsamples = self._getnsamples(filename)
        acls = self._getacls(filename, start_index=nsamples//2)
        is_burned_in = nacl(nsamples, acls, self._nacls)
        # stack the burn in results into an nparams x nchains array
        burn_in_per_chain = numpy.stack(list(is_burned_in.values())).all(
            axis=0)
        # store
        test = 'nacl'
        self.test_is_burned_in[test] = burn_in_per_chain
        try:
            burn_in_iter = self.test_burn_in_iteration[test]
        except KeyError:
            # hasn't been stored yet
            burn_in_iter = numpy.repeat(NOT_BURNED_IN_ITER, self.nchains)
            self.test_burn_in_iteration[test] = burn_in_iter
        burn_in_iter[burn_in_per_chain] = self._index2iter(filename,
                                                           nsamples//2)
        # add the status for each parameter as additional information
        self.test_aux_info[test] = is_burned_in

    def evaluate(self, filename):
        """Runs all of the burn-in tests."""
        # evaluate all the tests
        for tst in self.do_tests:
            logging.info("Evaluating %s burn-in test", tst)
            getattr(self, tst)(filename)
        # evaluate each chain at a time
        for ci in range(self.nchains):
            # some tests (like halfchain) just store a single bool for all
            # chains
            tibi = {t: r[ci] if isinstance(r, numpy.ndarray) else r
                    for t, r in self.test_is_burned_in.items()}
            tbi = {t: r[ci] if isinstance(r, numpy.ndarray) else r
                   for t, r in self.test_burn_in_iteration.items()}
            is_burned_in, burn_in_iter = evaluate_tests(self.burn_in_test,
                                                        tibi, tbi)
            self.is_burned_in[ci] = is_burned_in
            self.burn_in_iteration[ci] = burn_in_iter
        logging.info("Number of chains burned in: %i of %i",
                     self.is_burned_in.sum(), self.nchains)

    def write(self, fp, path=None):
        """Writes burn-in info to an open HDF file.

        Parameters
        ----------
        fp : pycbc.inference.io.base.BaseInferenceFile
            Open HDF file to write the data to. The HDF file should be an
            instance of a pycbc BaseInferenceFile.
        path : str, optional
            Path in the HDF file to write the data to. Default is (None) is
            to write to the path given by the file's ``sampler_group``
            attribute.
        """
        if path is None:
            path = fp.sampler_group
        super(MCMCBurnInTests, self).write(fp, path)
        # add number of chains burned in as additional metadata
        fp.write_data('nchains_burned_in', self.is_burned_in.sum(), path)


class MultiTemperedMCMCBurnInTests(MCMCBurnInTests):
    """Adds support for multiple temperatures to
    :py:class:`MCMCBurnInTests`.
    """

    def _getacls(self, filename, start_index):
        """Convenience function for calculating acls for the given filename.

        This function is used by the ``n_acl`` burn-in test. That function
        expects the returned ``acls`` dict to just report a single ACL for
        each parameter. Since multi-tempered samplers return an array of ACLs
        for each parameter instead, this takes the max over the array before
        returning.

        Since we calculate the acls, this will also store it to the sampler.

        Parameters
        ----------
        filename : str
            Name of the file to retrieve samples from.
        start_index : int
            Index to start calculating ACLs.

        Returns
        -------
        dict :
            Dictionary of parameter names -> array giving ACL for each chain.
        """
        acls = super(MultiTemperedMCMCBurnInTests, self)._getacls(
            filename, start_index)
        # acls will have shape ntemps x nchains, flatten to nchains
        return {param: vals.max(axis=0) for (param, vals) in acls.items()}

    def _getlogposts(self, filename):
        """Convenience function for retrieving log posteriors.

        This just gets the coldest temperature chain, and returns arrays with
        shape nwalkers x niterations, so the parent class can run the same
        ``posterior_step`` function.
        """
        return _multitemper_getlogposts(self.sampler, filename)


class EnsembleMCMCBurnInTests(BaseBurnInTests):
    """Provides methods for estimating burn-in of an ensemble MCMC."""

    available_tests = ('halfchain', 'min_iterations', 'max_posterior',
                       'posterior_step', 'nacl', 'ks_test',
                       )

    def __init__(self, sampler, burn_in_test, **kwargs):
        super(EnsembleMCMCBurnInTests, self).__init__(
            sampler, burn_in_test, **kwargs)
        # for kstest
        self._ksthreshold = float(kwargs.pop('ks_threshold', 0.9))

    def burn_in_index(self, filename):
        """The burn in index (retrieved from the iteration)."""
        if self.is_burned_in:
            index = self._iter2index(filename, self.burn_in_iteration)
        else:
            index = NOT_BURNED_IN_ITER
        return index

    def max_posterior(self, filename):
        """Applies max posterior test to self."""
        logposts = self._getlogposts(filename)
        burn_in_idx, is_burned_in = max_posterior(logposts, self._ndim)
        all_burned_in = is_burned_in.all()
        if all_burned_in:
            burn_in_iter = self._index2iter(filename, burn_in_idx.max())
        else:
            burn_in_iter = NOT_BURNED_IN_ITER
        # store
        test = 'max_posterior'
        self.test_is_burned_in[test] = all_burned_in
        self.test_burn_in_iteration[test] = burn_in_iter
        aux = self._getaux(test)
        # additional info
        aux['iteration_per_walker'] = self._index2iter(filename, burn_in_idx)
        aux['status_per_walker'] = is_burned_in

    def posterior_step(self, filename):
        """Applies the posterior-step test."""
        logposts = self._getlogposts(filename)
        burn_in_idx = numpy.array([posterior_step(logps, self._ndim)
                                   for logps in logposts])
        burn_in_iters = self._index2iter(filename, burn_in_idx)
        # this test cannot determine when something will burn in
        # only when it was not burned in in the past
        test = 'posterior_step'
        self.test_is_burned_in[test] = True
        self.test_burn_in_iteration[test] = burn_in_iters.max()
        # store the iteration per walker as additional info
        aux = self._getaux(test)
        aux['iteration_per_walker'] = burn_in_iters

    def nacl(self, filename):
        """Applies the :py:func:`nacl` test."""
        nsamples = self._getnsamples(filename)
        acls = self._getacls(filename, start_index=nsamples//2)
        is_burned_in = nacl(nsamples, acls, self._nacls)
        all_burned_in = all(is_burned_in.values())
        if all_burned_in:
            burn_in_iter = self._index2iter(filename, nsamples//2)
        else:
            burn_in_iter = NOT_BURNED_IN_ITER
        # store
        test = 'nacl'
        self.test_is_burned_in[test] = all_burned_in
        self.test_burn_in_iteration[test] = burn_in_iter
        # store the status per parameter as additional info
        aux = self._getaux(test)
        aux['status_per_parameter'] = is_burned_in

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
        is_burned_in = all(is_the_same.values())
        if is_burned_in:
            burn_in_iter = self._index2iter(filename, int(nsamples//2))
        else:
            burn_in_iter = NOT_BURNED_IN_ITER
        # store
        test = 'ks_test'
        self.test_is_burned_in[test] = is_burned_in
        self.test_burn_in_iteration[test] = burn_in_iter
        # store the test per parameter as additional info
        aux = self._getaux(test)
        aux['status_per_parameter'] = is_the_same

    def evaluate(self, filename):
        """Runs all of the burn-in tests."""
        # evaluate all the tests
        for tst in self.do_tests:
            logging.info("Evaluating %s burn-in test", tst)
            getattr(self, tst)(filename)
        is_burned_in, burn_in_iter = evaluate_tests(
            self.burn_in_test, self.test_is_burned_in,
            self.test_burn_in_iteration)
        self.is_burned_in = is_burned_in
        self.burn_in_iteration = burn_in_iter
        logging.info("Is burned in: %r", self.is_burned_in)
        if self.is_burned_in:
            logging.info("Burn-in iteration: %i",
                         int(self.burn_in_iteration))

    @staticmethod
    def _extra_tests_from_config(cp, section, tag):
        """Loads the ks test settings from the config file."""
        kwargs = {}
        if cp.has_option_tag(section, 'ks-threshold', tag):
            kwargs['ks_threshold'] = float(
                cp.get_opt_tag(section, 'ks-threshold', tag))
        return kwargs


class EnsembleMultiTemperedMCMCBurnInTests(EnsembleMCMCBurnInTests):
    """Adds support for multiple temperatures to
    :py:class:`EnsembleMCMCBurnInTests`.
    """

    def _getacls(self, filename, start_index):
        """Convenience function for calculating acls for the given filename.

        This function is used by the ``n_acl`` burn-in test. That function
        expects the returned ``acls`` dict to just report a single ACL for
        each parameter. Since multi-tempered samplers return an array of ACLs
        for each parameter instead, this takes the max over the array before
        returning.

        Since we calculate the acls, this will also store it to the sampler.
        """
        acls = super(EnsembleMultiTemperedMCMCBurnInTests, self)._getacls(
            filename, start_index)
        # return the max for each parameter
        return {param: vals.max() for (param, vals) in acls.items()}

    def _getlogposts(self, filename):
        """Convenience function for retrieving log posteriors.

        This just gets the coldest temperature chain, and returns arrays with
        shape nwalkers x niterations, so the parent class can run the same
        ``posterior_step`` function.
        """
        return _multitemper_getlogposts(self.sampler, filename)


def _multitemper_getlogposts(sampler, filename):
    """Retrieve log posteriors for multi tempered samplers."""
    with sampler.io(filename, 'r') as fp:
        samples = fp.read_raw_samples(
            ['loglikelihood', 'logprior'], thin_start=0, thin_interval=1,
            temps=0, flatten=False)
        # reshape to drop the first dimension
        for (stat, arr) in samples.items():
            _, nwalkers, niterations = arr.shape
            samples[stat] = arr.reshape((nwalkers, niterations))
        logposts = samples['loglikelihood'] + samples['logprior']
    return logposts
