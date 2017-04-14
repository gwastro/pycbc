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

def gelman_rubin_walk(x, seg_length, seg_stride, end_idx, seg_start=0):

    # lists to hold statistic and end index
    stats = []
    ends = []

    # get the beginning of all segments
    starts = numpy.arange(seg_start, end_idx, seg_stride)

    # loop over all segments
    for start in starts:

        # find the end of the first segment
        x_start_end = int(start + seg_length)

        # get first segment
        #x_start = [xx[start:x_start_end] for xx in x]
        x_start = [xx[0:x_start_end] for xx in x]

        # compute statistic
        stats.append((gelman_rubin(x_start)))

        # store end of first segment
        ends.append(x_start_end)

    return numpy.array(starts), numpy.array(ends), numpy.array(stats)

def gelman_rubin(chains):
    """ Calculates the univariate Gelman-Rubin convergence statistic.

    Parameters
    ----------
    chains : iterable
        An iterable of numpy.array instances that contain the samples
        for each chain. Each chain has shape (nparameters, niterations).
    """

def gelman_rubin(chains, auto_burn_in=True):
    """ chains is a list of numpy.array, ie. [chain1, chain2].
    each chain is shape (nparameters, niterations)
    """

    # remove first have of samples
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

    # calculate W
    # this will have shape (nparameters, nparameters)
    w = numpy.zeros(chains_covs[0].shape)
    for i, row in enumerate(chains_covs[0]):
        for j, _ in enumerate(row):
            w[i, j] = numpy.mean(chains_covs[:, i, j])
    if nparameters == 1:
        w = w.reshape((1, 1))

    # calculate B
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

    # get means
    # this will have shape (nparameters)
    mu_hat = numpy.mean(means, axis=1)

    # get variance of variance
    # this will have shape (nparameters)
    s = numpy.var(var, axis=1)

    # get factors in variance of V calculation
    # this will have shape (nparameters)
    k = 2 * b_diag**2 / (nchains - 1)
    mid_term = numpy.cov(var, means**2)[nparameters:2 * nparameters, 0:nparameters].T
    end_term = numpy.cov(var, means)[nparameters:2 * nparameters, 0:nparameters].T
    wb = niterations / nchains * numpy.diag(mid_term - 2 * mu_hat * end_term)

    # get V
    # this will have shape (nparameters)
    v = (niterations - 1.) * w_diag / niterations + (1. + 1. / nchains) * b_diag / niterations

    # get variance of V
    # this will have shape (nparameters)
    var_v = ((niterations - 1.)**2 * s + (1. + 1. / nchains)**2 * k + \
             2. * (niterations - 1.) * (1. + 1. / nchains) * wb) / niterations**2

    # get degrees of freedom
    # this will have shape (nparameters)
    dof = (2. * v**2) / var_v

    # more degrees of freedom
    # this will have shape (nparameters)
    df_adj = (dof + 3.) / (dof + 1.)
    b_dof = nchains - 1
    w_dof = (2. * w_diag**2) / s

    # estimate R
    # this will have shape (nparameters)
    r2_fixed = (niterations - 1.) / niterations
    r2_random = (1. + 1. / nchains) * (1. / niterations) * (b_diag / w_diag)
    r2_estimate = numpy.sqrt(r2_fixed + r2_random)

    # PSRF
    # this will have shape (nparameters)
    psrf = r2_estimate * df_adj

    return psrf

