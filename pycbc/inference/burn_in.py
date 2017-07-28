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
    array :
        Array of indices giving the burn-in index for each chain.
    """
    # get the posteriors
    # Note: multi-tempered samplers should just return the coldest chain by
    # default
    chain_stats = sampler.read_samples(fp, ['loglr', 'prior'],
        samples_group=fp.stats_group, thin_interval=1, thin_start=0,
        thin_end=None, flatten=False)
    chain_posteriors = chain_stats['loglr'] + chain_stats['prior']
    dim = len(fp.variable_args)
    # find the posterior to compare against
    max_p = chain_posteriors.max()
    criteria = max_p - dim/2
    nwalkers = chain_posteriors.shape[-2]
    niterations = chain_posteriors.shape[-1]
    burn_in_idx = numpy.repeat(niterations, nwalkers).astype(int)
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
    return burn_in_idx


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
    array :
        Array of indices giving the burn-in index for each chain.
    """
    # get the posteriors
    # Note: multi-tempered samplers should just return the coldest chain by
    # default
    chain_stats = sampler.read_samples(fp, ['loglr', 'prior'],
        samples_group=fp.stats_group, thin_interval=1, thin_start=0,
        thin_end=None, flatten=False)
    chain_posteriors = chain_stats['loglr'] + chain_stats['prior']
    nwalkers = chain_posteriors.shape[-2]
    dim = len(fp.variable_args)
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
    return burn_in_idx


def half_chain(sampler, fp):
    """Takes the second half of the iterations as post-burn in.

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
    array :
        Array of indices giving the burn-in index for each chain.
    """ 
    nwalkers = sampler.nwalkers
    niterations = fp.niterations
    return numpy.repeat(niterations/2, nwalkers).astype(int)


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
    array :
        Array of indices giving the burn-in index for each chain.
    """
    sampler.burn_in()
    return sampler.burn_in_iterations


burn_in_functions = {
    'max_posterior': max_posterior,
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
        array :
            Array of indices giving the burn-in index for each chain.
        """
        # if the file already has burn in iterations saved, use those as a
        # base
        try:
            burnidx = fp['burn_in_iterations'][:]
        except KeyError:
            # just use the minimum
            burnidx = numpy.repeat(self.min_iterations, fp.nwalkers)
        if self.burn_in_functions != {}:
            newidx = numpy.vstack([func(sampler, fp)
                for func in self.burn_in_functions.values()]).max(axis=0)
            mask = burnidx < newidx
            burnidx[mask] = newidx[mask]
        burnidx[burnidx < self.min_iterations] = self.min_iterations
        return burnidx

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
        """
        burnidx = self.evaluate(sampler, fp)
        sampler.burn_in_iterations = burnidx
        sampler.write_burn_in_iterations(fp, burnidx)

