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
        Sampler instance to determine burn in for.
    fp : InferenceFile
        Open inference hdf file containing the samples to load for determing
        burn in.

    Returns
    -------
    burnin_idx :
        Array of indices giving the burn-in index for each chain.
    """
    # get the posteriors
    chain_stats = sampler.read_samples(fp, ['loglr', 'prior'],
        samples_group=fp.stats_group, thin_interval=1, thin_start=0,
        thin_end=None, flatten=False)
    chain_posteriors = chain_stats['loglr'] + chain_stats['prior']
    dim = len(sampler.variable_args)
    # find the posterior to compare against
    maxP = chain_posteriors.max()
    criteria = maxP - dim/2
    nwalkers = chain_posteriors.shape[-2]
    niterations = chain_posteriors.shape[-1]
    burnin_idx = numpy.repeat(niterations, nwalkers).astype(int)
    for ii in range(nwalkers):
        chain = chain_posteriors[...,ii,:]
        idx = numpy.where(chain >= criteria)[-1]
        if idx.size != 0:
            burnin_idx[ii] = idx[0]
    return burnin_idx


def posterior_step(sampler, fp):
    """Burn in based on the last time a chain made a jump > dim/2.

    Parameters
    ----------
    sampler : pycbc.inference.sampler
        Sampler instance to determine burn in for.
    fp : InferenceFile
        Open inference hdf file containing the samples to load for determing
        burn in.

    Returns
    -------
    burnin_idx :
        Array of indices giving the burn-in index for each chain.
    """
    # get the posteriors
    chain_stats = sampler.read_samples(fp, ['loglr', 'prior'],
        samples_group=fp.stats_group, thin_interval=1, thin_start=0,
        thin_end=None, flatten=False)
    chain_posteriors = chain_stats['loglr'] + chain_stats['prior']
    dim = len(sampler.variable_args)
    burnin_idx = numpy.zeros(niterations, nwalkers).astype(int)
    for ii in range(nwalkers):
        chain = chain_posterior[...,ii,:]
        criteria = chain.max() - dim/2.
        dp = numpy.diff(chain)
        idx = numpy.where(dp >= criteria)[-1]
        if idx.size != 0:
            burnin_idx[ii] = idx[-1] + 1
    return burnin_idx


def half_chain(sampler, fp):
    """Takes the second half of the iterations as post-burn in.
    """ 
    nwalkers = sampler.nwalkers
    niterations = fp.niterations
    return numpy.repeat(niterations/2, nwalkers).astype(int)


def use_sampler(sampler, fp):
    """Uses the sampler's burn_in function.
    """
    sampler.burn_in()
    return numpy.repeat(sampler.burn_in_iterations,
                        sampler.nwalkers).astype(int)


burnin_functions = {
    'max_posterior': posterior_based,
    'posterior_step': posterior_step,
    'half_chain': half_chain,
    'use_sampler': use_sampler,
    }

class BurnIn(object):
    """Class to estimate the number of burn in iterations.
    """

    def __init__(self, burnin_functions, min_iterations=0):
        if burnin_functions is None:
            burnin_functions = []
        self.min_iterations = min_iterations
        self.burnin_functions = {fname: burnin_funcs[fname]
                                 for fname in burnin_functions}

    def evaluate(self, sampler, fp):
        """Evaluates sampler's chains to find burn in.
        """
        # if the file already has burn in iterations saved, use those as a
        # base
        try:
            burnidx = fp['burn_in_iterations'][:]
        except KeyError:
            # just use the minimum
            burnidx = numpy.repeat(self.min_iterations, sampler.nwalkers)
        if self.burnin_functions != {}:
            newidx = numpy.vstack([func(sampler, fp)
                for func in self.burnin_functions.values()]).max(axis=0)
            mask = burnidx < newidx
            burnidx[mask] = newidx[mask]
        burnidx[burnidx < self.min_iterations] = self.min_iterations
        return burnidx

    def update(self, sampler, fp):
        """Evaluates burn in saves updated indices to the given file.
        """
        burnidx = self.evaluate(sampler, fp)
        sampler.write_burnin_iterations(fp, burnidx)

