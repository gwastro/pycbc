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

def gelman_rubin(chains):
    """ Calculates the Gelman-Rubin convergence statistic.

    Parameters
    ----------
    chains : iterable
        An iterable of numpy.array instances that contain the samples
        for each chain.
    """

    # calculate mean of each chain
    chains_means = numpy.array([chain.mean() for chain in chains])

    # calculate overall mean
    overall_mean = chains_means.mean()

    # calculate variance of each chain
    chains_vars = numpy.array([chain.var() for chain in chains])

    # calculate the between-chain variance
    n = len(chains[0])
    m = len(chains)
    b = n / (m - 1.0) * sum([(chain_mean - overall_mean)**2
                             for chain_mean in chains_means])

    # calculate the within-chain variance
    w = 1.0 / m * sum([chain_var for chain_var in chains_vars])

    # calculate the pooled variance
    v = (n - 1.0) / m * w + (m + 1.0) / (m * n) * b

    # calculate the potential scale reduction factor
    psrf = numpy.sqrt(v / w)

    return psrf
