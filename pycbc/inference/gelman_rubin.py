# Copyright (C) 2017  Christopher M. Biwer
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
""" This modules provides functions for evaluating the Gelman-Rubin convergence
diagnostic statistic.
"""

import numpy

def walk(chains, start, end, step):
    """ Calculates Gelman-Rubin conervergence statistic along chains of data.
    This function will advance along the chains and calculate the
    statistic for each step.

    Parameters
    ----------
    chains : iterable
        An iterable of numpy.array instances that contain the samples
        for each chain. Each chain has shape (nparameters, niterations).
    start : float
        Start index of blocks to calculate all statistics.
    end : float
        Last index of blocks to calculate statistics.
    step : float
        Step size to take for next block.

    Returns
    -------
    starts : numpy.array
        1-D array of start indexes of calculations.
    ends : numpy.array
        1-D array of end indexes of caluclations.
    stats : numpy.array
        Array with convergence statistic. It has
        shape (nparameters, ncalculations).
    """

    # get number of chains, parameters, and iterations
    chains = numpy.array(chains)
    _, nparameters, _ = chains.shape

    # get end index of blocks
    ends = numpy.arange(start, end, step)
    stats = numpy.zeros((nparameters, len(ends)))

    # get start index of blocks
    starts = numpy.array(len(ends) * [start])

    # loop over end indexes and calculate statistic
    for i, e in enumerate(ends):
        tmp = chains[:,:,0:e]
        stats[:, i] = gelman_rubin(tmp)

    return starts, ends, stats

def gelman_rubin(chains, auto_burn_in=True):
    """ Calculates the univariate Gelman-Rubin convergence statistic
    which compares the evolution of multiple chains in a Markov-Chain Monte
    Carlo process and computes their difference to determine their convergence.
    The between-chain and within-chain variances are computed for each sampling
    parameter, and a weighted combination of the two is used to determine the
    convergence. As the chains converge, the point scale reduction factor
    should go to 1.

    Parameters
    ----------
    chains : iterable
        An iterable of numpy.array instances that contain the samples
        for each chain. Each chain has shape (nparameters, niterations).
    auto_burn_in : bool
        If True, then only use later half of samples provided.

    Returns
    -------
    psrf : numpy.array
        A numpy.array of shape (nparameters) that has the point estimates of
        the potential scale reduction factor.
    """

    # remove first half of samples
    # this will have shape (nchains, nparameters, niterations)
    if auto_burn_in:
        _, _, niterations = numpy.array(chains).shape
        chains = numpy.array([chain[:, niterations / 2 + 1:]
                              for chain in chains])

    # get number of chains, parameters, and iterations
    chains = numpy.array(chains)
    nchains, nparameters, niterations = chains.shape

    # calculate the covariance matrix for each chain
    # this will have shape (nchains, nparameters, nparameters)
    chains_covs = numpy.array([numpy.cov(chain) for chain in chains])
    if nparameters == 1:
        chains_covs = chains_covs.reshape((nchains, 1, 1))

    # calculate W the within-chain variance
    # this will have shape (nparameters, nparameters)
    w = numpy.zeros(chains_covs[0].shape)
    for i, row in enumerate(chains_covs[0]):
        for j, _ in enumerate(row):
            w[i, j] = numpy.mean(chains_covs[:, i, j])
    if nparameters == 1:
        w = w.reshape((1, 1))

    # calculate B the between-chain variance
    # this will have shape (nparameters, nparameters)
    means = numpy.zeros((nparameters, nchains))
    for i, chain in enumerate(chains):
        means[:, i] = numpy.mean(chain, axis=1).transpose()
    b = niterations * numpy.cov(means)
    if nparameters == 1:
        b = b.reshape((1, 1))

    # get diagonal elements of W and B
    # these will have shape (nparameters)
    w_diag = numpy.diag(w)
    b_diag = numpy.diag(b)

    # get variance for each chain
    # this will have shape (nparameters, nchains)
    var = numpy.zeros((nparameters, nchains))
    for i, chain_cov in enumerate(chains_covs):
        var[:, i] = numpy.diag(chain_cov)

    # get mean of means
    # this will have shape (nparameters)
    mu_hat = numpy.mean(means, axis=1)

    # get variance of variances
    # this will have shape (nparameters)
    s = numpy.var(var, axis=1)

    # get V the combined variance of all chains
    # this will have shape (nparameters)
    v = (niterations - 1.) * w_diag / niterations + (1. + 1. / nchains) * b_diag / niterations

    # get factors in variance of V calculation
    # this will have shape (nparameters)
    k = 2 * b_diag**2 / (nchains - 1)
    mid_term = numpy.cov(var, means**2)[nparameters:2 * nparameters, 0:nparameters].T
    end_term = numpy.cov(var, means)[nparameters:2 * nparameters, 0:nparameters].T
    wb = niterations / nchains * numpy.diag(mid_term - 2 * mu_hat * end_term)

    # get variance of V
    # this will have shape (nparameters)
    var_v = ((niterations - 1.)**2 * s + (1. + 1. / nchains)**2 * k + \
             2. * (niterations - 1.) * (1. + 1. / nchains) * wb) / niterations**2

    # get degrees of freedom
    # this will have shape (nparameters)
    dof = (2. * v**2) / var_v

    # more degrees of freedom factors
    # this will have shape (nparameters)
    df_adj = (dof + 3.) / (dof + 1.)

    # estimate R
    # this will have shape (nparameters)
    r2_fixed = (niterations - 1.) / niterations
    r2_random = (1. + 1. / nchains) * (1. / niterations) * (b_diag / w_diag)
    r2_estimate = r2_fixed + r2_random

    # calculate PSRF the potential scale reduction factor
    # this will have shape (nparameters)
    psrf = numpy.sqrt(r2_estimate * df_adj)

    return psrf

