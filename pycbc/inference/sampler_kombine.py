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
This modules provides classes and functions for using the kombine sampler
packages for parameter estimation.
"""

from pycbc.inference.sampler_base import BaseMCMCSampler

#
# =============================================================================
#
#                                   Samplers
#
# =============================================================================
#


class KombineSampler(BaseMCMCSampler):
    """This class is used to construct the MCMC sampler from the kombine
    package.

    Parameters
    ----------
    likelihood_evaluator : likelihood class
        An instance of the likelihood class from the
        pycbc.inference.likelihood module.
    nwalkers : int
        Number of walkers to use in sampler.
    transd : bool
        If True, the sampler will operate across parameter spaces using a
        kombine.clustered_kde.TransdimensionalKDE proposal distribution. In
        this mode a masked array with samples in each of the possible sets of
        dimensions must be given for the initial ensemble distribution.
    processes : {None, int}
        Number of processes to use with multiprocessing. If None, all available
        cores are used.
    min_burn_in : {None, int}
        Set the minimum number of burn in iterations to use. If None,
        `burn_in_iterations` will be initialized to `0`.
    """
    name = "kombine"

    def __init__(self, likelihood_evaluator, nwalkers, transd=False,
                 processes=None, min_burn_in=None):
        try:
            import kombine
        except ImportError:
            raise ImportError("kombine is not installed.")

        # construct sampler for use in KombineSampler
        ndim = len(likelihood_evaluator.waveform_generator.variable_args)
        sampler = kombine.Sampler(nwalkers, ndim, likelihood_evaluator,
                                  transd=transd, processes=processes)
        # initialize
        super(KombineSampler, self).__init__(sampler, likelihood_evaluator,
                                             min_burn_in=min_burn_in)
        self._nwalkers = nwalkers

    @classmethod
    def from_cli(cls, opts, likelihood_evaluator):
        """Create an instance of this sampler from the given command-line
        options.

        Parameters
        ----------
        opts : ArgumentParser options
            The options to parse.
        likelihood_evaluator : LikelihoodEvaluator
            The likelihood evaluator to use with the sampler.

        Returns
        -------
        KombineSampler
            A kombine sampler initialized based on the given arguments.
        """
        return cls(likelihood_evaluator, opts.nwalkers,
                   processes=opts.nprocesses, min_burn_in=opts.min_burn_in)

    def run(self, niterations, **kwargs):
        """Advance the sampler for a number of samples.

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
        blob0 = None
        if self.burn_in_iterations == 0:
            # no burn in, use the initial positions
            p0 = self.p0
            if self.likelihood_evaluator.return_meta:
                blob0 = [self.likelihood_evaluator(p0[wi, :])[1]
                         for wi in range(self.nwalkers)]
        else:
            p0 = None
            # kombine requires blob data to be specified
            if self.likelihood_evaluator.return_meta:
                blob0 = self._sampler.blobs[-1]
        kwargs['blob0'] = blob0
        res = self._sampler.run_mcmc(niterations, p0=p0, **kwargs)
        p, lnpost, lnprop = res[0], res[1], res[2]
        # update the positions
        self._pos = p
        return p, lnpost, lnprop

    @property
    def lnpost(self):
        """ Get the natural logarithm of the likelihood as an
        nwalkers x niterations array.
        """
        # kombine returns niterations x nwaklers
        return self._sampler.lnpost.transpose()

    @property
    def chain(self):
        """Get all past samples as an nwalker x niterations x ndim array."""
        # kombine returns niterations x nwalkers x ndim
        return self._sampler.chain.transpose((1, 0, 2))

    def burn_in(self):
        """Use kombine's `burnin` routine to advance the sampler.

        If a minimum number of burn-in iterations was specified, this will run
        the burn-in until it has advanced at least as many steps as desired.
        The initial positions (p0) must be set prior to running.

        For more details, see `kombine.sampler.burnin`.

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
        # check that we haven't already burned in
        if self.pos is not None:
            raise ValueError("burn-in already run")
        # run once
        p0 = self.p0
        if self.likelihood_evaluator.return_meta:
            blob0 = [self.likelihood_evaluator(p0[wi, :])[1]
                     for wi in range(self.nwalkers)]
        else:
            blob0 = None
        res = self._sampler.burnin(self.p0, blob0=blob0)
        p, post, q = res[0], res[1], res[2]
        # continue running until minimum burn in is satistfied
        while self.niterations < self.burn_in_iterations:
            p0 = p
            res = self._sampler.burnin(p0)
            p, post, q = res[0], res[1], res[2]
            # update position
            self._pos = p
        self.burn_in_iterations = self.niterations
        return p, post, q
