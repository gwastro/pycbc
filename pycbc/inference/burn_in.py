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

import numpy
from scipy.stats import ks_2samp

def ks_test(sampler, fp, threshold=0.9):
    """Burn in based on whether the p-value of the KS test between the samples
    at the last iteration and the samples midway along the chain for each
    parameter is > ``threshold``.

    Parameters
    ----------
    sampler : pycbc.inference.sampler
        Sampler to determine burn in for. May be either an instance of a
        `inference.sampler`, or the class itself.
    fp : InferenceFile
        Open inference hdf file containing the samples to load for determing
        burn in.
    threshold : float
        The thershold to use for the p-value. Default is 0.9.

    Returns
    -------
    burn_in_idx : array
        Array of indices giving the burn-in index for each chain.
    is_burned_in : array
        Array of booleans indicating whether each chain is burned in.
    """
    nwalkers = fp.nwalkers
    niterations = fp.niterations
    # Create a dictionary which would have keys are the variable args and values
    # are booleans indicating whether the p-value for the parameters satisfies
    # the KS test
    is_burned_in_param = {}
    # iterate over the parameters
    for param in fp.variable_args:
        # read samples for the parameter from the last iteration of the chain
        samples_last_iter = sampler.read_samples(fp, param, iteration=-1,
                                                 flatten=True)[param]
        # read samples for the parameter from the iteration midway along the chain
        samples_chain_midpt = sampler.read_samples(fp, param, iteration=int(niterations/2),
                                                  flatten=True)[param]
        _, p_value = ks_2samp(samples_last_iter, samples_chain_midpt)
        # check if p_value is > than the desired range
        is_burned_in_param[param] = p_value > threshold
    # The chains are burned in if the p-value of the KS test lies in the range [0.1,0.9]
    # for all the parameters. If the KS test is passed, the chains have burned in at their
    # mid-way point. 
    if all(is_burned_in_param.values()):
        is_burned_in = numpy.ones(nwalkers, dtype=bool)
        burn_in_idx = numpy.repeat(niterations/2, nwalkers).astype(int)
    else:
        is_burned_in = numpy.zeros(nwalkers, dtype=bool)
        burn_in_idx = numpy.repeat(niterations, nwalkers).astype(int)
    return burn_in_idx, is_burned_in


def n_acl(sampler, fp, nacls=10):
    """Burn in based on ACL.

    The sampler is considered burned in if the number of itertions is >=
    ``nacls`` times the maximum ACL over all parameters, as measured from the
    first iteration.

    Parameters
    ----------
    sampler : pycbc.inference.sampler
        Sampler to determine burn in for. May be either an instance of a
        `inference.sampler`, or the class itself.
    fp : InferenceFile
        Open inference hdf file containing the samples to load for determing
        burn in.
    nacls : int
        Number of ACLs to use for burn in. Default is 10.

    Returns
    -------
    burn_in_idx : array
        Array of indices giving the burn-in index for each chain. By definition
        of this function, all chains reach burn in at the same iteration. Thus
        the returned array is the burn-in index repeated by the number of
        chains.
    is_burned_in : array
        Array of booleans indicating whether each chain is burned in. Since
        all chains obtain burn in at the same time, this is either an array
        of all False or True.
    """
    acl = numpy.array(sampler.compute_acls(fp, start_index=0).values()).max()
    burn_idx = nacls * acl
    is_burned_in = burn_idx < fp.niterations
    if not is_burned_in:
        burn_idx = fp.niterations
    nwalkers = fp.nwalkers
    return numpy.repeat(burn_idx, nwalkers).astype(int), \
           numpy.repeat(is_burned_in, nwalkers).astype(bool)


def max_posterior(sampler, fp):
    """Burn in based on samples being within dim/2 of maximum posterior.

    Parameters
    ----------
    sampler : pycbc.inference.sampler
        Sampler to determine burn in for. May be either an instance of a
        `inference.sampler`, or the class itself.
    fp : InferenceFile
        Open inference hdf file containing the samples to load for determing
        burn in.

    Returns
    -------
    burn_in_idx : array
        Array of indices giving the burn-in index for each chain.
    is_burned_in : array
        Array of booleans indicating whether each chain is burned in.
    """
    # get the posteriors
    # Note: multi-tempered samplers should just return the coldest chain by
    # default
    chain_stats = sampler.read_samples(fp, ['loglr', 'prior'],
        samples_group=fp.stats_group, thin_interval=1, thin_start=0,
        thin_end=None, flatten=False)
    chain_posteriors = chain_stats['loglr'] + chain_stats['prior']
    dim = float(len(fp.variable_args))
    # find the posterior to compare against
    max_p = chain_posteriors.max()
    criteria = max_p - dim/2
    nwalkers = chain_posteriors.shape[-2]
    niterations = chain_posteriors.shape[-1]
    burn_in_idx = numpy.repeat(niterations, nwalkers).astype(int)
    is_burned_in = numpy.zeros(nwalkers, dtype=bool)
    # find the first iteration in each chain where the logplr has exceeded
    # max_p - dim/2
    for ii in range(nwalkers):
        chain = chain_posteriors[...,ii,:]
        # numpy.where will return a tuple with multiple arrays if the chain is
        # more than 1D (which can happen for multi-tempered samplers). Always
        # taking the last array ensures we are looking at the indices that
        # count out iterations
        idx = numpy.where(chain >= criteria)[-1]
        if idx.size != 0:
            burn_in_idx[ii] = idx[0]
            is_burned_in[ii] = True
    return burn_in_idx, is_burned_in


def acl_or_max_posterior(sampler, fp):
    """Burn in that uses the minimum of `n_acl` and `max_posteior`.
    """
    acl_bidx, acl_ibi = n_acl(sampler, fp)
    maxp_bidx, maxp_ibi = max_posterior(sampler, fp)
    burn_in_idx = numpy.min([acl_bidx, maxp_bidx], axis=0)
    is_burned_in = acl_ibi | maxp_ibi
    return burn_in_idx, is_burned_in

def posterior_step(sampler, fp):
    """Burn in based on the last time a chain made a jump > dim/2.

    Parameters
    ----------
    sampler : pycbc.inference.sampler
        Sampler to determine burn in for. May be either an instance of a
        `inference.sampler`, or the class itself.
    fp : InferenceFile
        Open inference hdf file containing the samples to load for determing
        burn in.

    Returns
    -------
    burn_in_idx : array
        Array of indices giving the burn-in index for each chain.
    is_burned_in : array
        Array of booleans indicating whether each chain is burned in.
        By definition of this function, all values are set to True.
    """
    # get the posteriors
    # Note: multi-tempered samplers should just return the coldest chain by
    # default
    chain_stats = sampler.read_samples(fp, ['loglr', 'prior'],
        samples_group=fp.stats_group, thin_interval=1, thin_start=0,
        thin_end=None, flatten=False)
    chain_posteriors = chain_stats['loglr'] + chain_stats['prior']
    nwalkers = chain_posteriors.shape[-2]
    dim = float(len(fp.variable_args))
    burn_in_idx = numpy.zeros(nwalkers).astype(int)
    criteria = dim/2.
    # find the last iteration in each chain where the logplr has
    # jumped by more than dim/2
    for ii in range(nwalkers):
        chain = chain_posteriors[...,ii,:]
        dp = abs(numpy.diff(chain))
        idx = numpy.where(dp >= criteria)[-1]
        if idx.size != 0:
            burn_in_idx[ii] = idx[-1] + 1
    return burn_in_idx, numpy.ones(nwalkers, dtype=bool)


def half_chain(sampler, fp):
    """Takes the second half of the iterations as post-burn in.

    Parameters
    ----------
    sampler : pycbc.inference.sampler
        This option is not used; it is just here give consistent API as the
        other burn in functions.
    fp : InferenceFile
        Open inference hdf file containing the samples to load for determing
        burn in.

    Returns
    -------
    burn_in_idx : array
        Array of indices giving the burn-in index for each chain.
    is_burned_in : array
        Array of booleans indicating whether each chain is burned in.
        By definition of this function, all values are set to True.
    """
    nwalkers = fp.nwalkers
    niterations = fp.niterations
    return numpy.repeat(niterations/2, nwalkers).astype(int), \
           numpy.ones(nwalkers, dtype=bool)


def use_sampler(sampler, fp=None):
    """Uses the sampler's burn_in function.

    Parameters
    ----------
    sampler : pycbc.inference.sampler
        Sampler to determine burn in for. Must be an instance of an
        `inference.sampler` that has a `burn_in` function.
    fp : InferenceFile, optional
        This option is not used; it is just here give consistent API as the
        other burn in functions.

    Returns
    -------
    burn_in_idx : array
        Array of indices giving the burn-in index for each chain.
    is_burned_in : array
        Array of booleans indicating whether each chain is burned in.
        Since the sampler's burn in function will run until all chains
        are burned, all values are set to True.
    """
    sampler.burn_in()
    return sampler.burn_in_iterations, \
           numpy.ones(len(sampler.burn_in_iterations), dtype=bool)


burn_in_functions = {
    'ks_test': ks_test,
    'n_acl': n_acl,
    'max_posterior': max_posterior,
    'acl_or_max_posterior': acl_or_max_posterior,
    'posterior_step': posterior_step,
    'half_chain': half_chain,
    'use_sampler': use_sampler,
    }

class BurnIn(object):
    """Class to estimate the number of burn in iterations.

    Parameters
    ----------
    function_names : list, optional
        List of name of burn in functions to use. All names in the provided
        list muset be in the `burn_in_functions` dict. If none provided, will
        use no burn-in functions.
    min_iterations : int, optional
        Minimum number of burn in iterations to use. The burn in iterations
        returned by evaluate will be the maximum of this value
        and the values returned by the burn in functions provided in
        `function_names`. Default is 0.

    Examples
    --------
    Initialize a `BurnIn` instance that will use `max_posterior` and
    `posterior_step` as the burn in criteria:

    >>> from pycbc import inference
    >>> burn_in = inference.BurnIn(['max_posterior', 'posterior_step'])

    Use this `BurnIn` instance to find the burn-in iteration of each walker
    in an inference result file:

    >>> from pycbc.io import InferenceFile
    >>> fp = InferenceFile('inference.hdf', 'r')
    >>> burn_in.evaluate(inference.samplers[fp.sampler_name], fp)
    array([11486, 11983, 11894, ..., 11793, 11888, 11981])

    """

    def __init__(self, function_names, min_iterations=0):
        if function_names is None:
            function_names = []
        self.min_iterations = min_iterations
        self.burn_in_functions = {fname: burn_in_functions[fname]
                                  for fname in function_names}

    def evaluate(self, sampler, fp):
        """Evaluates sampler's chains to find burn in.

        Parameters
        ----------
        sampler : pycbc.inference.sampler
            Sampler to determine burn in for. May be either an instance of a
            `inference.sampler`, or the class itself.
        fp : InferenceFile
            Open inference hdf file containing the samples to load for
            determing burn in.

        Returns
        -------
        burnidx : array
            Array of indices giving the burn-in index for each chain.
        is_burned_in : array
            Array of booleans indicating whether each chain is burned in.
        """
        # if the number of iterations is < than the minimium desired,
        # just return the number of iterations and all False
        if fp.niterations < self.min_iterations:
            return numpy.repeat(self.min_iterations, fp.nwalkers), \
                   numpy.zeros(fp.nwalkers, dtype=bool)
        # if the file already has burn in iterations saved, use those as a
        # base
        try:
            burnidx = fp['burn_in_iterations'][:]
        except KeyError:
            # just use the minimum
            burnidx = numpy.repeat(self.min_iterations, fp.nwalkers)
        # start by assuming is burned in; the &= below will make this false
        # if any test yields false
        is_burned_in = numpy.ones(fp.nwalkers, dtype=bool)
        if self.burn_in_functions != {}:
            newidx = []
            for func in self.burn_in_functions.values():
                idx, state = func(sampler, fp)
                newidx.append(idx)
                is_burned_in &= state
            newidx = numpy.vstack(newidx).max(axis=0)
            # update the burn in idx if any test yields a larger iteration
            mask = burnidx < newidx
            burnidx[mask] = newidx[mask]
        # if any burn-in idx are less than the min iterations, set to the
        # min iterations
        burnidx[burnidx < self.min_iterations] = self.min_iterations
        return burnidx, is_burned_in

    def update(self, sampler, fp):
        """Evaluates burn in and saves the updated indices to the given file.

        Parameters
        ----------
        sampler : pycbc.inference.sampler
            Sampler to determine burn in for. May be either an instance of a
            `inference.sampler`, or the class itself.
        fp : InferenceFile
            Open inference hdf file containing the samples to load for
            determing burn in.

        Returns
        -------
        burnidx : array
            Array of indices giving the burn-in index for each chain.
        is_burned_in : array
            Array of booleans indicating whether each chain is burned in.
        """
        burnidx, is_burned_in = self.evaluate(sampler, fp)
        sampler.burn_in_iterations = burnidx
        sampler.write_burn_in_iterations(fp, burnidx, is_burned_in)
        return burnidx, is_burned_in
