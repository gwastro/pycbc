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
    def lnpost(self):
        """ This function should return the natural logarithm of the likelihood
        as an niterations x nwalker array.
        """
        return ValueError("lnpost function not set.")

    @property
    def chain(self):
        """ This function should return the past samples as a
        niterations x nwalker x ndim array.
        """
        return ValueError("chain function not set.")

    def burn_in(self, initial_values):
        """ This function should burn in the sampler.
        """
        raise ValueError("This sampler has no burn_in function.")

    def run(self, niterations):
        """ This function should run the sampler.
        """
        raise ValueError("run function not set.")

class _BaseMCMCSampler(_BaseSampler):
    """ This class is used to construct the MCMC sampler from the kombine-like
    packages.

    Parameters
    ----------
    likelihood_evaluator : likelihood class
        An instance of the likelihood class from the
        pycbc.inference.likelihood module.
    sampler : sampler instance
        An instance of an MCMC sampler similar to kombine or emcee.
    """
    def __init__(self, sampler, likelihood_evaluator):

        self._sampler = sampler 
        # initialize
        super(_BaseMCMCSampler, self).__init__(likelihood_evaluator)

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


class KombineSampler(_BaseMCMCSampler):
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
        Number of dimensions in the parameter space. If transd is True this is
        the number of unique dimensions across the parameter spaces.
    transd : bool
        If True, the sampler will operate across parameter spaces using a
        kombine.clustered_kde.TransdimensionalKDE proposal distribution. In
        this mode a masked array with samples in each of the possible sets of
        dimensions must be given for the initial ensemble distribution.
    processes : {None, int}
        Number of processes to use with multiprocessing. If None, all available
        cores are used.
    """

    def __init__(self, likelihood_evaluator, nwalkers=0, ndim=0,
                        transd=False, processes=None):

        try:
            import kombine
        except ImportError:
            raise ImportError("kombine is not installed.")

        # construct sampler for use in KombineSampler
        sampler = kombine.Sampler(nwalkers, ndim, likelihood_evaluator,
                                          transd=transd, processes=processes)

        # initialize
        super(KombineSampler, self).__init__(sampler, likelihood_evaluator)

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

    @property
    def lnpost(self):
        """ Get the natural logarithm of the likelihood as an 
        niterations x nwalkers array.
        """
        return self._sampler.lnpost

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
            res = self._sampler.burnin(initial_values)
            if len(res) == 4:
                p, post, q, _ = res
            else:
                p, post, q = res
            self.burn_in_iterations = self.chain.shape[0]
        else:
            raise ValueError("Burn in has already been performed")
        return p, post, q


class EmceeEnsembleSampler(_BaseMCMCSampler):
    """This class is used to construct an MCMC sampler from the emcee
    package's EnsembleSampler.

    Parameters
    ----------
    likelihood_evaluator : likelihood class
        An instance of the likelihood class from the
        pycbc.inference.likelihood module.
    nwalkers : int
        Number of walkers to use in sampler.
    ndim : int
        Number of dimensions in the parameter space. If transd is True this is
        the number of unique dimensions across the parameter spaces.
    processes : {None, int}
        Number of processes to use with multiprocessing. If None, all available
        cores are used.
    """

    def __init__(self, likelihood_evaluator, nwalkers=0, ndim=0,
                        processes=None):

        try:
            import emcee
        except ImportError:
            raise ImportError("emcee is not installed.")

        # initialize the pool to use
        if processes == 1:
            pool = None
        else:
            pool = emcee.interruptible_pool.InterruptiblePool(
                processes=processes)

        # construct the sampler
        sampler = emcee.EnsembleSampler(nwalkers, ndim, likelihood_evaluator,
            pool=pool)

        # initialize
        super(EmceeEnsembleSampler, self).__init__(sampler,
            likelihood_evaluator)

    @property
    def lnpost(self):
        """Get the natural logarithm of the likelihood as an 
        niterations x nwalkers array.
        """
        return self._sampler.lnprobability

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
        try:
            p0 = kwargs.pop('p0')
        except KeyError:
            raise ValueError("must provide a p0 in the keyword arguments")
        return self._sampler.run_mcmc(p0, niterations, **kwargs)


samplers = {
    "kombine" : KombineSampler,
    "emcee" : EmceeEnsembleSampler,
}
