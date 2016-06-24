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
    """ Base container class for running the inference sampler that will
    generate the posterior distributions.

    Parameters
    ----------
    sampler : sampler class
        An instance of the sampler class from its package.
    """

    def __init__(self, likelihood_evaluator):
        self.likelihood_evaluator = likelihood_evaluator
        self.burn_in_iterations = 0

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
        """ This function should burn in the sampler.
        """
        raise ValueError("burn_in function not set.")

    def run(self, niterations):
        """ This function should run the sampler.
        """
        raise ValueError("run function not set.")

class KombineSampler(_BaseSampler):
    """ This class is used to construct the MCMC sampler from the kombine
    package.

    Parameters
    ----------
    likelihood_evaluator : likelihood class
        An instance of the likelihood class from the
        pycbc.inference.likelihood module.
    nwalkers : int
        Number of walkers to use in sampler.
    ndim : int
        Number of dimensions in the parameter space. If kwargs["transd"] is
        True this is the number of unique dimensions across the parameter
        spaces.
    processes : {None, int}
        Number of processes to use with multiprocessing. If None, all available
        cores are used.
    **kwargs
        Any other keyword arguments are passed to the initialization of
        Kombine.Sampler.
    """

    def __init__(self, likelihood_evaluator, nwalkers=0, ndim=0,
                                        processes=None, **kwargs):

        try:
            import kombine
        except ImportError:
            raise ImportError("kombine is not installed.")

        # construct sampler for use in KombineSampler
        self._sampler = kombine.Sampler(nwalkers, ndim, likelihood_evaluator,
                                          processes=processes, **kwargs)

        # initialize
        super(KombineSampler, self).__init__(likelihood_evaluator)

    @property
    def acceptance_fraction(self):
        """ Get the fraction of walkers that accepted each step as an arary.
        """
        return self._sampler.acceptance_fraction

    @property
    def chain(self):
        """ Get all past samples as an niterations x nwalker x ndim array.
        """
        return self._sampler.chain

    def burn_in(self, initial_values):
        """ Evolve an ensemble until the acceptance rate becomes roughly
        constant. This is done by splitting acceptances in half and checking
        for statistical consistency. This isn't guaranteed to return a fully
        burned-in ensemble, but usually does.

        Parameters
        ----------
        initial_values : numpy.array
            An nwalkers x ndim array of initial values for walkers.

        Returns
        -------
        p : numpy.array
            An array of current walker positions with shape (nwalkers, ndim).
        lnpost : numpy.array
            The list of log posterior probabilities for the walkers at
            positions p, with shape (nwalkers, ndim).
        lnprop : numpy.array
            The list of log proposal densities for the walkers at positions p,
            with shape (nwalkers, ndim).
        """
        if self.burn_in_iterations == 0:
            p, post, q = self._sampler.burnin(initial_values)
            self.burn_in_iterations = self.chain.shape[0]
        else:
            raise ValueError("Burn in has already been performed")
        return p, post, q

    def run(self, niterations, **kwargs):
        """ Advance the sampler for a number of samples.

        Parameters
        ----------
        niterations : int
            Number of samples to get from sampler.

        Returns
        -------
        p : numpy.array
            An array of current walker positions with shape (nwalkers, ndim).
        lnpost : numpy.array
            The list of log posterior probabilities for the walkers at
            positions p, with shape (nwalkers, ndim).
        lnprop : numpy.array
            The list of log proposal densities for the walkers at positions p,
            with shape (nwalkers, ndim).
        """
        return self._sampler.run_mcmc(niterations, **kwargs)

samplers = {
    "kombine" : KombineSampler,
}
