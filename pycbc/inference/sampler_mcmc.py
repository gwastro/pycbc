# Copyright (C) 2017  Vivien Raymond
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
This modules provides classes and functions for using a MCMC sampler
for parameter estimation.
"""

import numpy
import logging
from pycbc.inference.sampler_base import BaseMCMCSampler

#
# =============================================================================
#
#                                   Samplers
#
# =============================================================================
#

class MCMCSampler(BaseMCMCSampler):
    """This class is used to construct the MCMC sampler.

    Parameters
    ----------
    likelihood_evaluator : LikelihoodEvaluator
        An instance of a pycbc.inference.likelihood evaluator.
    """
    name = "mcmc"

    def __init__(self, likelihood_evaluator):
        self._chain = []
        self._blobs = []
        # Using p0 to store the last sample would require to store sparately the
        # last posterior value.
        self._lastsample = []
        self._lastblob = []
        sampler = self
        # initialize
        super(MCMCSampler, self).__init__(
              sampler, likelihood_evaluator)
        self.dtype=numpy.dtype([(name, None) for name in
                                            ('lnpost',)+self.sampling_args])
        # Harcoding the number of walkers to 1.
        # nwalkers should not be a BaseMCMCSampler property.
        self._nwalkers = 1

    @classmethod
    def from_cli(cls, opts, likelihood_evaluator, pool=None, likelihood_call=None):
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
        MCMCSampler
            A MCMC sampler initialized based on the given arguments.
        """
        return cls(likelihood_evaluator)

    @property
    def chain(self):
        """This function should return the past samples as a
        [additional dimensions x] niterations x ndim array, where ndim are the
        number of sampling args, niterations the number of iterations, and
        additional dimensions are any additional dimensions used by the
        sampler (e.g, walkers, temperatures).
        """
        #Adding the nwalkers dimention, and converting to an ndarray.
        return self._chain[list(self.sampling_args)].view(numpy.float).reshape(
                                            (1,) + self._chain.shape + (-1,))
        # Copy needed to avoid numpy 1.13 warning

    @property
    def blobs(self):
        """This function should return the blobs with a shape
        nwalkers x niteration as requested by the BaseMCMCSampler class.
        """
        #Adding the nwalkers dimention.
        return [self._blobs]

    def clear_chain(self):
        """This function should clear the current chain of samples from memory.
        """
        # store the iteration that the clear is occuring on
        self.lastclear = self.niterations
        self._chain = []
        self._blobs = []

    @property
    def niterations(self):
        """Get the current number of iterations."""
        return len(self._chain)+self.lastclear

    @property
    def lnpost(self):
        """This function should return the natural logarithm of the likelihood
        function used by the sampler as an
        [additional dimensions] x niterations array.
        """
        return self._chain['lnpost']

    @property
    def p0(self):
        """Since this is just a single chain, forces p0 to have shape (nparams,).
        """
        p0 = super(MCMCSampler, self).p0
        if p0.ndim == 2:
            p0 = p0.flatten()
        return p0

    def run(self, niterations):
        """This function should run the sampler.
        """

        if not self._lastsample:
            # first time running, use the initial positions
            # set_p0() was called in pycbc_inference, so self.p0 is set
            result = self.likelihood_evaluator(self.p0)
            try:
                logplr, blob = result
            except TypeError:
                # likelihood evaluator doesn't return blobs
                logplr = result
                blob = None

            logging.info("Starting logplr value %f", logplr)

            start_sample = numpy.insert(self.p0,0,logplr)
            start = 0

        else:
            start_sample = self._lastsample
            start = -1 # So restarts do not re-write the same sample.
            blob = self._lastblob

        self._chain = numpy.empty(niterations,dtype=self.dtype)
        # numpy >=1.14 only accepts tuples
        self._chain[start] = tuple(start_sample)
        self._blobs = [None]*niterations
        self._blobs[start] = blob

        for i in range(start, niterations-1):

            logplr_old = self._chain['lnpost'][i]
            # As _chain is a structured numpy array and self.sampling_args is a
            # tuple, a list() conversion is needed here.
            # This is not ideal being in the inner loop.
            samples = self._chain[list(self.sampling_args)][i]

            # Dummy proposal
            samples_prop = [sample + numpy.random.normal(loc=0.0, scale=1.0)
                            for sample in samples]

            result = self.likelihood_evaluator(samples_prop)
            try:
                logplr_prop, blob = result
            except TypeError:
                # likelihood evaluator doesn't return blobs
                logplr_prop = result
                blob = None

            acceptance_ratio=numpy.exp(logplr_prop - logplr_old)
            u=numpy.random.uniform()
            if acceptance_ratio >= u:
                self._chain[i+1]=numpy.insert(samples_prop,0,logplr_prop)
                self._blobs[i+1]=blob
                logging.info("Step %i, acceptance ratio %f >= %f, accepted",
                                        i+1, acceptance_ratio, u)
            else:
                self._chain[i+1]=self._chain[i]
                self._blobs[i+1]=self._blobs[i]
                logging.info("Step %i, acceptance ratio %f < %f, rejected",
                                        i+1, acceptance_ratio, u)

        self._lastsample = self._chain[-1]
        self._lastblob = self._blobs[-1]

        return

    def write_acceptance_fraction(self, fp, start_iteration=None,
                                  max_iterations=None):
        # Overwrite the BaseMCMCSampler function since we do not compute the
        # required acceptance fraction.
        pass
