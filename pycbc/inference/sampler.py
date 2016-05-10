# Copyright (C) 2016  Christopher M. Biwer
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
This modules provides classes and functions for using different sampler
packages for parameter estimation.
"""

class _BaseSampler(object):
    """ Base container class for running the MCMC sampler.

    Parameters
    ----------
    sampler : sampler class
        An instance of the sampler class from its package.
    """

    def __init__(self, sampler):
        self.sampler = sampler

    @property
    def acceptance_fraction(self):
        """ This function should return the fraction of walkers that accepted
        each step as an array.
        """
        return ValueError("acceptance_fraction function not set.")

    @property
    def chain(self):
        """ This function should return the past samples as a
        initerations x nwalker x ndim array.
        """
        return ValueError("chain function not set.")

    def burn_in(self, initial_values):
        """ This function should burn in the MCMC.
        """
        raise ValueError("burn_in function not set.")

    def run_mcmc(self, niterations):
        """ This function should run the MCMC for the number of samples.
        """
        raise ValueError("run_mcmc function not set.")

class KombineSampler(_BaseSampler):
    """ This class is used to construct the sampler from the kombine package.

    Parameters
    ----------
    likelihood_evaluator : likelihood class
        An instance of the likelihood class from the
        pycbc.inference.likelihood module.
    nwalkers : int
        Number of walkers to use in MCMC.
    ndim : int
        Number of dimensions in the parameter space. If transd is True this is
        the number of unique dimensions across the parameter spaces.
    transd : bool
        If True, the sampler will operate across parameter spaces using a
        kombine.clustered_kde.TransdimensionalKDE proposal distribution. In
        this mode a masked array with samples in each of the possible sets of
        dimensions must be given for the initial ensemble distribution.
    nprocesses : {None, int}
        Number of processes to use with multiprocessing. If None, all available
        cores are used.
    """

    def __init__(self, likelihood_evaluator, nwalkers=0, ndim=0,
                        transd=False, nprocess=None):

        try:
            import kombine
        except ImportError:
            raise ImportError("kombine is not installed.")

        # construct sampler 
        sampler = kombine.Sampler(nwalkers, ndim, self.likelihood_evaluator,
                                          transd=transd, processes=nprocesses)
        super(KombineSampler, self).__init__(self, sampler)

    @property
    def acceptance_fraction(self):
        """ Get the fraction of walkers that accepted each step as an arary.
        """
        return self.sampler.acceptance_fraction

    @property
    def chain(self):
        """ Get all past samples as an niterations x nwalker x ndim array.
        """
        return self.sampler.chain

    def burn_in(self, initial_values):
        """ Evolve an ensemble until the acceptance rate becomes roughly
        constant. This is done by splitting acceptances in half and checking
        for statistical consistency. This isn’t guaranteed to return a fully
        burned-in ensemble, but usually does.

        Parameters
        ----------
        initial_values : numpy.array
            An nwalkers x ndim array of initial values for walkers.
        """
        self.sampler.burn_in(initial_values)

    def run_mcmc(self, niterations):
        """ Advance the MCMC for a number of samples.

        Parameters
        ----------
        niterations : int
            Number of samples to get from MCMC.
        """
        self.sampler.run_mcmc(niterations)

samplers = {
    "kombine" : KombineSampler,
}
