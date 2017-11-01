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
from pycbc.inference.sampler_base import _BaseSampler

#
# =============================================================================
#
#                                   Samplers
#
# =============================================================================
#

class MCMCSampler(_BaseSampler):
    """This class is used to construct the MCMC sampler.

    Parameters
    ----------
    likelihood_evaluator : LikelihoodEvaluator
        An instance of a pycbc.inference.likelihood evaluator.
    niterations : int
        Number of iterations to use in sampler.
    """
    name = "mcmc"

    def __init__(self, likelihood_evaluator, niterations):
        self.likelihood_evaluator = likelihood_evaluator
        self._niterations = niterations
        sampling_args=likelihood_evaluator.sampling_args
        ndim = len(sampling_args)
        dtype=numpy.dtype([(name, None) for name in
                            ['lnpost','lnlike']+sampling_args])
        self.samples_chain = numpy.full([niterations,ndim+2],
                                        numpy.nan,dtype=dtype)

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
        return cls(likelihood_evaluator, opts.niterations)

    @property
    def ifos(self):
        """Returns the ifos that were sampled over."""
        return self.likelihood_evaluator.waveform_generator.detector_names

    @property
    def variable_args(self):
        """Returns the variable args used by the likelihood evaluator.
        """
        return self.likelihood_evaluator.variable_args

    @property
    def sampling_args(self):
        """Returns the sampling args used by the likelihood evaluator.
        """
        return self.likelihood_evaluator.sampling_args

    @property
    def chain(self):
        """This function should return the past samples as a
        [additional dimensions x] niterations x ndim array, where ndim are the
        number of variable args, niterations the number of iterations, and
        additional dimensions are any additional dimensions used by the
        sampler (e.g, walkers, temperatures).
        """
        return self.samples_chain[:self.niterations]

    @property
    def samples(self):
        """This function should return the past samples as a [additional
        dimensions x] niterations field array, where the fields are union
        of the sampling args and the variable args.
        """
        return self.samples_chain[:self.niterations]

    @property
    def clear_chain(self):
        """This function should clear the current chain of samples from memory.
        """
        del self.samples_chain
        return None

    @property
    def niterations(self):
        """Get the current number of iterations."""
        return len(numpy.where(self.samples_chain[0]!=numpy.nan)[0])

    @property
    def acceptance_fraction(self):
        """This function should return the fraction of walkers that accepted
        each step as an array.
        """
        return NotImplementedError("acceptance_fraction function not set.")

    @property
    def lnpost(self):
        """This function should return the natural logarithm of the likelihood
        function used by the sampler as an
        [additional dimensions] x niterations array.
        """
        return self.samples_chain[:self.niterations]['lnpost']

    @property
    def likelihood_stats(self):
        """This function should return the prior and likelihood ratio of
        samples as an [additional dimensions] x niterations
        array. If the likelihood evaluator did not return that info to the
        sampler, it should return None.
        """
        return NotImplementedError("likelihood stats not set")

    def burn_in(self, initial_values):
        """This function should burn in the sampler.
        """
        raise NotImplementedError("This sampler has no burn_in function.")

    def set_p0(self, samples=None, prior=None):
        """Sets the initial position of the MCMC.

        Parameters
        ----------
        samples : FieldArray, optional
            Use the given samples to set the initial positions. The samples
            will be transformed to the likelihood evaluator's `sampling_args`
            space.
        prior : PriorEvaluator, optional
            Use the given prior to set the initial positions rather than
            `likelihood_evaultor`'s prior.

        Returns
        -------
        p0 : array
            An ndim array of the initial positions that were set.
        """
        # create a (nwalker, ndim) array for initial positions

        ndim = len(self.variable_args)
        p0 = numpy.ones(ndim)
        # if samples are given then use those as initial positions
        if samples is not None:
            # transform to sampling parameter space
            samples = self.likelihood_evaluator.apply_sampling_transforms(
                samples)
        # draw random samples if samples are not provided
        else:
            samples = self.likelihood_evaluator.prior_rvs(size=1,prior=prior)
        # convert to 1D array
        for i, param in enumerate(self.sampling_args):
            p0[i] = samples[param]
        self._p0 = p0
        return p0

    @property
    def p0(self):
        if self._p0 is None:
            raise ValueError("initial positions not set; run set_p0")
        return self._p0

    def run(self, niterations):
        """This function should run the sampler.
        """
        if self.niterations == 0:
            # first time running, use the initial positions
            samples = self.p0
            loglr = self.likelihood_evaluator.loglr(p0)
            logplr = self.likelihood_evaluator.prior(p0) + loglr
            self.samples_chain[0]=numpy.insert(samples,0,[logplr,loglr])
            self.niterations=1

        for i in range(niterations):

            logplr_old,loglr_old = self.samples_chain[i][:2]
            samples = self.samples_chain[i][2:]
            samples_prop = samples + numpy.random.normal(loc=0.0, scale=0.1,
                                                        size=len(samples))
            loglr_prop = self.likelihood_evaluator.loglr(samples_prop)
            logplr_prop = self.likelihood_evaluator.prior(samples_prop) \
                        + loglr_prop
            if logplr_prop / logplr_old > numpy.random.uniform():
                self.samples_chain[i+1]=numpy.insert(samples_prop,0,
                                        [logplr_prop,loglr_prop])
            else:
                self.samples_chain[i+1]=self.samples_chain[i]

        return samples_prop, logplr_prop, loglr_prop

    @classmethod
    def calculate_logevidence(cls, fp):
        """This function should calculate the log evidence and its error using
        the results in the given file. If the sampler does not support evidence
        calculation, then this will raise a NotImplementedError.
        """
        raise NotImplementedError("this sampler does not support evidence "
                                  "calculation")

    # write and read functions
    def write_metadata(self, fp):
        """Writes metadata about this sampler to the given file. Metadata is
        written to the file's `attrs`.

        Parameters
        ----------
        fp : InferenceFile
            A file handler to an open inference file.
        """
        fp.attrs['sampler'] = self.name
        fp.attrs['likelihood_evaluator'] = self.likelihood_evaluator.name
        fp.attrs['ifos'] = self.ifos
        fp.attrs['variable_args'] = list(self.variable_args)
        fp.attrs['sampling_args'] = list(self.sampling_args)
        fp.attrs["niterations"] = self.niterations
        fp.attrs["lognl"] = self.likelihood_evaluator.lognl
        sargs = self.likelihood_evaluator.waveform_generator.static_args
        fp.attrs["static_args"] = sargs.keys()
        for arg, val in sargs.items():
            fp.attrs[arg] = val

    @staticmethod
    def write_logevidence(fp, lnz, dlnz):
        """Writes the given log evidence and its error to the given file.
        Results are saved to the file's 'log_evidence' and 'dlog_evidence'
        attributes.

        Parameters
        ----------
        fp : InferenceFile
            A file handler to an open inference file.
        lnz : float
            The log of the evidence.
        dlnz : float
            The error in the estimate of the log evidence.
        """
        fp.attrs['log_evidence'] = lnz
        fp.attrs['dlog_evidence'] = dlnz

    @staticmethod
    def write_burn_in_iterations(fp, burn_in_iterations):
        """Writes the burn in iterations to the given file.

        Parameters
        ----------
        fp : InferenceFile
            A file handler to an open inference file.
        burn_in_iterations : array
            Array of values giving the iteration of the burn in of each walker.
        """
        try:
            fp['burn_in_iterations'][:] = burn_in_iterations
        except KeyError:
            fp['burn_in_iterations'] = burn_in_iterations
        fp.attrs['burn_in_iterations'] = burn_in_iterations.max()
