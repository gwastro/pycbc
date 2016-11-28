# Copyright (C) 2016  Collin Capano
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
This modules provides classes and functions for using the emcee sampler
packages for parameter estimation.
"""

import numpy
from pycbc.inference.sampler_base import BaseMCMCSampler
from pycbc.io import WaveformArray, FieldArray
from pycbc.filter import autocorrelation

#
# =============================================================================
#
#                                   Samplers
#
# =============================================================================
#


class EmceeEnsembleSampler(BaseMCMCSampler):
    """This class is used to construct an MCMC sampler from the emcee
    package's EnsembleSampler.

    Parameters
    ----------
    likelihood_evaluator : likelihood class
        An instance of the likelihood class from the
        pycbc.inference.likelihood module.
    nwalkers : int
        Number of walkers to use in sampler.
    processes : {None, int}
        Number of processes to use with multiprocessing. If None, all available
        cores are used.
    burn_in_iterations : {None, int}
        Set the number of burn in iterations to use. If None,
        `burn_in_ieterations` will be initialized to 0.
    """
    name = "emcee"

    def __init__(self, likelihood_evaluator, nwalkers, processes=None,
                 burn_in_iterations=None):

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
        ndim = len(likelihood_evaluator.waveform_generator.variable_args)
        sampler = emcee.EnsembleSampler(nwalkers, ndim, likelihood_evaluator,
                                        pool=pool)
        # initialize
        super(EmceeEnsembleSampler, self).__init__(
              sampler, likelihood_evaluator, min_burn_in=burn_in_iterations)
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
        EmceeEnsembleSampler
            An emcee sampler initialized based on the given arguments.
        """
        # check that if not skipping burn in, more than one burn in iteration
        # has been specified
        if not opts.skip_burn_in and (
                opts.min_burn_in is None or opts.min_burn_in == 0):
            raise ValueError("{name} requires that you provide a ".format(
                name=cls.name) + " non-zero --min-burn-in if not skipping "
                "burn-in")
        return cls(likelihood_evaluator, opts.nwalkers,
                   processes=opts.nprocesses,
                   burn_in_iterations=opts.min_burn_in)

    @property
    def lnpost(self):
        """Get the natural logarithm of the likelihood as an
        nwalkers x niterations array.
        """
        # emcee returns nwalkers x niterations
        return self._sampler.lnprobability

    @property
    def chain(self):
        """Get all past samples as an nwalker x niterations x ndim array."""
        # emcee returns the chain as nwalker x niterations x ndim
        return self._sampler.chain

    def run(self, niterations, **kwargs):
        """Advance the ensemble for a number of samples.

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
        rstate :
            The current state of the random number generator.
        """
        pos = self._pos
        if pos is None:
            pos = self.p0
        res = self._sampler.run_mcmc(pos, niterations, **kwargs)
        p, lnpost, rstate = res[0], res[1], res[2]
        # update the positions
        self._pos = p
        return p, lnpost, rstate

    def burn_in(self, **kwargs):
        """Advance the ensemble by the number of the sampler's
        `burn_in_iterations`.

        Returns
        -------
        p : numpy.array
            An array of current walker positions with shape (nwalkers, ndim).
        lnpost : {None, numpy.array}
            The list of log posterior probabilities for the walkers at
            positions p, with shape (nwalkers, ndim).
        rstate :
            The current state of the random number generator.
        """
        if self.burn_in_iterations == 0:
            raise ValueError("must specify a non-zero number of iterations "
                             "to burn in")
        return self.run(self.burn_in_iterations, **kwargs)

    # Emcee defines acceptance fraction differently, so have to override
    # write functions
    def write_acceptance_fraction(self, fp):
        """Write acceptance_fraction data to file. Results are written to
        `fp[acceptance_fraction]`.

        Parameters
        -----------
        fp : InferenceFile
            A file handler to an open inference file.
        """
        dataset_name = "acceptance_fraction"
        try:
            fp[dataset_name][:] = self.acceptance_fraction
        except KeyError:
            # dataset doesn't exist yet, create it
            fp[dataset_name] = self.acceptance_fraction

    @staticmethod
    def read_acceptance_fraction(fp, walkers=None):
        """Reads the acceptance fraction from the given file.

        Parameters
        -----------
        fp : InferenceFile
            An open file handler to read the samples from.
        walkers : {None, (list of) int}
            The walker index (or a list of indices) to retrieve. If None,
            samples from all walkers will be obtained.

        Returns
        -------
        array
            Array of acceptance fractions with shape (requested walkers,).
        """
        group = 'acceptance_fraction'
        if walkers is None:
            wmask = numpy.ones(fp.nwalkers, dtype=bool)
        else:
            wmask = numpy.zeros(fp.nwalkers, dtype=bool)
            wmask[walkers] = True
        return fp[group][wmask]

    def write_results(self, fp, max_iterations=None):
        """Writes metadata, samples, likelihood stats, and acceptance fraction
        to the given file. See the write function for each of those for
        details.

        Parameters
        -----------
        fp : InferenceFile
            A file handler to an open inference file.
        max_iterations : {None, int}
            If results have not previously been written to the
            file, new datasets will be created. By default, the size of these
            datasets will be whatever the length of the sampler's chain is at
            this point. If you intend to run more iterations in the future,
            set this value to that size so that the array in the file will be
            large enough to accomodate future data.
        """
        self.write_metadata(fp)
        self.write_chain(fp, max_iterations=max_iterations)
        self.write_likelihood_stats(fp, max_iterations=max_iterations)
        self.write_acceptance_fraction(fp)


class EmceePTSampler(BaseMCMCSampler):
    """This class is used to construct a parallel-tempered MCMC sampler from
    the emcee package's PTSampler.

    Parameters
    ----------
    likelihood_evaluator : likelihood class
        An instance of the likelihood class from the
        pycbc.inference.likelihood module.
    ntemps : int
        Number of temeratures to use in the sampler.
    nwalkers : int
        Number of walkers to use in sampler.
    processes : {None, int}
        Number of processes to use with multiprocessing. If None, all available
        cores are used.
    burn_in_iterations : {None, int}
        Set the number of burn in iterations to use. If None,
        `burn_in_ieterations` will be initialized to 0.
    """
    name = "emcee_pt"

    def __init__(self, likelihood_evaluator, ntemps, nwalkers, processes=None,
                 burn_in_iterations=None):

        try:
            import emcee
        except ImportError:
            raise ImportError("emcee is not installed.")

        # construct the sampler: PTSampler needs the likelihood and prior
        # functions separately
        likelihood_evaluator.set_callfunc('loglikelihood')

        # initialize the pool to use
        if processes == 1:
            pool = None
        else:
            pool = emcee.interruptible_pool.InterruptiblePool(
                processes=processes)

        ndim = len(likelihood_evaluator.waveform_generator.variable_args)
        sampler = emcee.PTSampler(ntemps, nwalkers, ndim,
                                  likelihood_evaluator,
                                  likelihood_evaluator._prior,
                                  pool=pool)
        # initialize
        super(EmceePTSampler, self).__init__(
              sampler, likelihood_evaluator, min_burn_in=burn_in_iterations)
        self._nwalkers = nwalkers
        self._ntemps = ntemps

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
        EmceePTSampler
            An emcee sampler initialized based on the given arguments.
        """
        # check that if not skipping burn in, more than one burn in iteration
        # has been specified
        if not opts.skip_burn_in and (
                opts.min_burn_in is None or opts.min_burn_in == 0):
            raise ValueError("%s requires that you provide a non-zero " % (
                cls.name) + "--min-burn-in if not skipping burn-in")
        return cls(likelihood_evaluator, opts.ntemps, opts.nwalkers,
                   processes=opts.nprocesses,
                   burn_in_iterations=opts.min_burn_in)

    @property
    def ntemps(self):
        return self._ntemps

    @property
    def chain(self):
        """Get all past samples as an ntemps x nwalker x niterations x ndim
        array.
        """
        # emcee returns the chain as ntemps x nwalker x niterations x ndim
        return self._sampler.chain

    @property
    def likelihood_stats(self):
        """Returns the log likelihood ratio and log prior as a FieldArray.
        The returned array has shape ntemps x nwalkers x niterations.
        """
        # likelihood has shape ntemps x nwalkers x niterations
        logl = self._sampler.lnlikelihood
        # get prior from posterior
        logp = self._sampler.lnprobability - logl
        # compute the likelihood ratio
        loglr = logl - self.likelihood_evaluator.lognl
        return FieldArray.from_kwargs(loglr=loglr, prior=logp)

    @property
    def lnpost(self):
        """Get the natural logarithm of the likelihood + the prior as an
        ntemps x nwalkers x niterations array.
        """
        # emcee returns ntemps x nwalkers x niterations
        return self._sampler.lnprobability

    def set_p0(self, prior_distributions):
        """Sets the initial position of the walkers.

        Parameters
        ----------
        prior_distributions : list
            A list of priors to retrieve random values from (the sort of
            thing returned by `prior.read_distributions_from_config`).

        Returns
        -------
        p0 : array
            An ntemps x nwalkers x ndim array of the initial positions that
            were set.
        """
        # loop over all walkers and then parameters
        # find the distribution that has that parameter in it and draw a
        # random value from the distribution
        ntemps = self.ntemps
        nwalkers = self.nwalkers
        ndim = len(self.variable_args)
        pmap = dict([[param, k] for k, param in enumerate(self.variable_args)])
        p0 = numpy.ones((ntemps, nwalkers, ndim))
        for dist in prior_distributions:
            ps = dist.rvs(size=ntemps*nwalkers).reshape(ntemps, nwalkers)
            for param in dist.params:
                p0[:, :, pmap[param]] = ps[param]
        self._p0 = p0
        return p0

    def run(self, niterations, **kwargs):
        """Advance the ensemble for a number of samples.

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
        rstate :
            The current state of the random number generator.
        """
        pos = self._pos
        if pos is None:
            pos = self.p0
        res = self._sampler.run_mcmc(pos, niterations, **kwargs)
        p, lnpost, rstate = res[0], res[1], res[2]
        # update the positions
        self._pos = p
        return p, lnpost, rstate

    def burn_in(self, **kwargs):
        """Advance the ensemble by the number of the sampler's
        `burn_in_iterations`.

        Returns
        -------
        p : numpy.array
            An array of current walker positions with shape (nwalkers, ndim).
        lnpost : {None, numpy.array}
            The list of log posterior probabilities for the walkers at
            positions p, with shape (nwalkers, ndim).
        rstate :
            The current state of the random number generator.
        """
        if self.burn_in_iterations == 0:
            raise ValueError("must specify a non-zero number of iterations "
                             "to burn in")
        return self.run(self.burn_in_iterations, **kwargs)

    # read/write functions

    # add ntemps and betas to metadata
    def write_metadata(self, fp):
        """Writes metadata about this sampler to the given file. Metadata is
        written to the file's `attrs`.

        Parameters
        ----------
        fp : InferenceFile
            A file handler to an open inference file.
        """
        super(EmceePTSampler, self).write_metadata(fp)
        fp.attrs["ntemps"] = self.ntemps
        fp.attrs["betas"] = self._sampler.betas

    def write_acceptance_fraction(self, fp):
        """Write acceptance_fraction data to file. Results are written to
        `fp[acceptance_fraction/temp{k}]` where k is the temperature.

        Parameters
        -----------
        fp : InferenceFile
            A file handler to an open inference file.
        """
        group = "acceptance_fraction/temp{tk}"
        # acf has shape ntemps x nwalkers
        acf = self.acceptance_fraction
        for tk in range(fp.ntemps):
            try:
                fp[group.format(tk=tk)][:] = acf[tk, :]
            except KeyError:
                # dataset doesn't exist yet, create it
                fp[group.format(tk=tk)] = acf[tk, :]

    @staticmethod
    def read_acceptance_fraction(fp, temps=None, walkers=None):
        """Reads the acceptance fraction from the given file.

        Parameters
        -----------
        fp : InferenceFile
            An open file handler to read the samples from.
        temps : {None, (list of) int}
            The temperature index (or a list of indices) to retrieve. If None,
            acfs from all temperatures and all walkers will be retrieved.
        walkers : {None, (list of) int}
            The walker index (or a list of indices) to retrieve. If None,
            samples from all walkers will be obtained.

        Returns
        -------
        array
            Array of acceptance fractions with shape (requested temps,
            requested walkers).
        """
        group = 'acceptance_fraction/temp{tk}'
        if temps is None:
            temps = numpy.arange(fp.ntemps)
        if walkers is None:
            wmask = numpy.ones(fp.nwalkers, dtype=bool)
        else:
            wmask = numpy.zeros(fp.nwalkers, dtype=bool)
            wmask[walkers] = True
        arrays = []
        for tk in temps:
            arrays.append(fp[group.format(tk=tk)][wmask])
        return numpy.vstack(arrays)

    def write_chain(self, fp, max_iterations=None):
        """Writes the samples from the current chain to the given file. Results
        are written to: `fp[fp.samples_group/{vararg}/temp{k}/walker{i}]`,
        where `{vararg}` is the name of a variable arg, `{k}` is a temperature
        index (smaller = colder), and `{i}` is the index of a walker.

        Parameters
        -----------
        fp : InferenceFile
            A file handler to an open inference file.
        max_iterations : {None, int}
            If samples have not previously been written to the file, a new
            dataset will be created. By default, the size of this dataset will
            be whatever the length of the sampler's chain is at this point. If
            you intend to run more iterations, set this value to that size so
            that the array in the file will be large enough to accomodate
            future data.
        """
        # chain is ntemps x nwalkers x niterations x ndim
        samples = self.chain
        ntemps, nwalkers, niterations, _ = samples.shape

        # map sample values to the values that were actually passed to the
        # waveform generator and prior evaluator
        samples = numpy.array(
            self.likelihood_evaluator._prior.apply_boundary_conditions(
            samples.transpose(3,0,1,2))).transpose(1,2,3,0)

        group = fp.samples_group + '/{name}/temp{tk}/walker{wi}'

        # create an empty array if desired, in case this is the first time
        # writing
        if max_iterations is not None:
            if max_iterations < niterations:
                raise IndexError("The provided max size is less than the "
                                 "number of iterations")
            out = numpy.zeros(max_iterations, dtype=samples.dtype)

        # create indices for faster sub-looping
        widx = numpy.arange(nwalkers)
        tidx = numpy.arange(ntemps)

        # loop over number of dimensions
        for pi, param in enumerate(self.variable_args):
            # loop over number of temps
            for tk in tidx:
                # loop over number of walkers
                for wi in widx:
                    dataset_name = group.format(name=param, tk=tk, wi=wi)
                    try:
                        fp[dataset_name][:niterations] = samples[tk, wi, :, pi]
                    except KeyError:
                        # dataset doesn't exist yet, see if a larger array is
                        # desired
                        if max_iterations is not None:
                            out[:niterations] = samples[tk, wi, :, pi]
                            fp[dataset_name] = out
                        else:
                            fp[dataset_name] = samples[tk, wi, :, pi]

    def write_likelihood_stats(self, fp, max_iterations=None):
        """Writes the given likelihood array to the given file. Results are
        written to: `fp[fp.stats_group/{field}/temp{k}/walker{i}]`, where
        `{field}` is the name of stat (`loglr`, `prior`), `{k}` is a
        temperature index (smaller = colder) and `{i}` is the index of a
        walker.

        Parameters
        -----------
        fp : InferenceFile
            A file handler to an open inference file.
        max_iterations : {None, int}
            See `write_chain` for details.
        """
        lls = self.likelihood_stats
        ntemps, nwalkers, niterations = lls.shape

        group = fp.stats_group + '/{field}/temp{tk}/walker{wi}'

        if max_iterations is not None and max_iterations < niterations:
            raise IndexError("The provided max size is less than the "
                             "number of iterations")

        # create indices for faster sub-looping
        widx = numpy.arange(nwalkers)
        tidx = numpy.arange(ntemps)

        # loop over stats
        for stat in lls.fieldnames:
            arr = lls[stat]
            # create an empty array if desired, in case this is the first time
            # writing
            if max_iterations is not None:
                out = numpy.zeros(max_iterations, dtype=arr.dtype)
            # loop over temps
            for tk in tidx:
                # loop over number of walkers
                for wi in widx:
                    dataset_name = group.format(field=stat, tk=tk, wi=wi)
                    try:
                        fp[dataset_name][:niterations] = arr[tk, wi, :]
                    except KeyError:
                        # dataset doesn't exist yet, see if a larger array is
                        # desired
                        if max_iterations is not None:
                            out[:niterations] = arr[tk, wi, :]
                            fp[dataset_name] = out
                        else:
                            fp[dataset_name] = arr[tk, wi, :]

    def write_results(self, fp, max_iterations=None):
        """Writes metadata, samples, lnpost, lnprior,  and acceptance fraction
        to the given file. See the write function for each of those for
        details.

        Parameters
        -----------
        fp : InferenceFile
            A file handler to an open inference file.
        max_iterations : {None, int}
            If results have not previously been written to the
            file, new datasets will be created. By default, the size of these
            datasets will be whatever the length of the sampler's chain is at
            this point. If you intend to run more iterations in the future,
            set this value to that size so that the array in the file will be
            large enough to accomodate future data.
        """
        self.write_metadata(fp)
        self.write_chain(fp, max_iterations=max_iterations)
        self.write_likelihood_stats(fp, max_iterations=max_iterations)
        self.write_acceptance_fraction(fp)

    @staticmethod
    def _read_fields(fp, fields_group, fields, array_class,
                     thin_start=None, thin_interval=None, thin_end=None,
                     iteration=None, temps=None, walkers=None, flatten=True):
        """Base function for reading samples and likelihood stats. See
        `read_samples` and `read_likelihood_stats` for details.

        Parameters
        -----------
        fp : InferenceFile
            An open file handler to read the samples from.
        fields_group : str
            The name of the group to retrieve the desired fields.
        fields : list
            The list of field names to retrieve. Must be names of groups in
            `fp[fields_group/]`.
        array_class : FieldArray or similar
            The type of array to return. Must have a `from_kwargs` attribute.

        For other details on keyword arguments, see `read_samples` and
        `read_likelihood_stats`.

        Returns
        -------
        array_class
            An instance of the given array class populated with values
            retrieved from the fields.
        """
        # walkers to load
        if walkers is None:
            walkers = range(fp.nwalkers)
        if isinstance(walkers, int):
            walkers = [walkers]

        # temperatures to load
        if temps is None:
            temps = 0
        if temps == 'all':
            temps = range(fp.ntemps)
        if isinstance(temps, int):
            temps = [temps]

        # get the slice to use
        if iteration is not None:
            get_index = iteration
        else:
            if thin_end is None:
                # use the number of current iterations
                thin_end = fp.niterations
            get_index = fp.get_slice(thin_start=thin_start, thin_end=thin_end,
                                     thin_interval=thin_interval)

        # load
        arrays = {}
        group = fields_group + '/{name}/temp{tk}/walker{wi}'
        for name in fields:
            these_arrays = numpy.array(
                [[fp[group.format(name=name, wi=wi, tk=tk)][get_index]
                 for wi in walkers]
                 for tk in temps])
            if flatten:
                these_arrays = these_arrays.flatten()
            arrays[name] = these_arrays
        return array_class.from_kwargs(**arrays)

    @classmethod
    def read_samples(cls, fp, parameters,
                     thin_start=None, thin_interval=None, thin_end=None,
                     iteration=None, temps=0, walkers=None, flatten=True,
                     samples_group=None, array_class=None):
        """Reads samples for the given parameter(s).

        Parameters
        -----------
        fp : InferenceFile
            An open file handler to read the samples from.
        parameters : (list of) strings
            The parameter(s) to retrieve. A parameter can be the name of any
            field in `fp[fp.samples_group]`, a virtual field or method of
            `WaveformArray` (as long as the file contains the necessary fields
            to derive the virtual field or method), and/or a function of
            these.
        thin_start : int
            Index of the sample to begin returning samples. Default is to read
            samples after burn in. To start from the beginning set thin_start
            to 0.
        thin_interval : int
            Interval to accept every i-th sample. Default is to use the
            `fp.acl`. If `fp.acl` is not set, then use all samples
            (set thin_interval to 1).
        thin_end : int
            Index of the last sample to read. If not given then
            `fp.niterations` is used.
        iteration : int
            Get a single iteration. If provided, will override the
            `thin_{start/interval/end}` arguments.
        walkers : {None, (list of) int}
            The walker index (or a list of indices) to retrieve. If None,
            samples from all walkers will be obtained.
        temps : {None, (list of) int, 'all'}
            The temperature index (or list of indices) to retrieve. If None,
            only samples from the coldest (= 0) temperature chain will be
            retrieved. To retrieve all temperates pass 'all', or a list of
            all of the temperatures.
        flatten : {True, bool}
            The returned array will be one dimensional, with all desired
            samples from all desired walkers concatenated together. If False,
            the returned array will have dimension requested temps x requested
            walkers x requested iterations.
        samples_group : {None, str}
            The group in `fp` from which to retrieve the parameter fields. If
            None, searches in `fp.samples_group`.
        array_class : {None, array class}
            The type of array to return. The class must have a `from_kwargs`
            class method and a `parse_parameters` method. If None, will return
            a WaveformArray.

        Returns
        -------
        array_class
            Samples for the given parameters, as an instance of a the given
            `array_class` (`WaveformArray` if `array_class` is None).
        """
        # get the group to load from
        if samples_group is None:
            samples_group = fp.samples_group
        # get the type of array class to use
        if array_class is None:
            array_class = WaveformArray
        # get the names of fields needed for the given parameters
        possible_fields = dict([[str(name), float]
                               for name in fp[fp.samples_group].keys()])
        loadfields = array_class.parse_parameters(
            parameters, possible_fields=possible_fields)
        return cls._read_fields(
                fp, samples_group, loadfields, array_class,
                thin_start=thin_start, thin_interval=thin_interval,
                thin_end=thin_end, iteration=iteration, temps=temps,
                walkers=walkers, flatten=flatten)

    @classmethod
    def read_likelihood_stats(
            cls, fp,
            thin_start=None, thin_interval=None, thin_end=None, iteration=None,
            temps=None, walkers=None, flatten=True, stats_group=None,
            array_class=None):
        """Reads the likelihood stats from the given file.

        Parameters
        -----------
        fp : InferenceFile
            An open file handler to read the stats from.
        thin_start : int
            Index of the sample to begin returning stats. Default is to read
            stats after burn in. To start from the beginning set thin_start
            to 0.
        thin_interval : int
            Interval to accept every i-th sample. Default is to use the
            `fp.acl`. If `fp.acl` is not set, then use all stats
            (set thin_interval to 1).
        thin_end : int
            Index of the last sample to read. If not given then
            `fp.niterations` is used.
        iteration : int
            Get a single iteration. If provided, will override the
            `thin_{start/interval/end}` arguments.
        temps : {None, (list of) int, 'all'}
            The temperature index (or list of indices) to retrieve. If None,
            only samples from the coldest (= 0) temperature chain will be
            retrieved. To retrieve all temperates pass 'all', or a list of
            all of the temperatures.
        walkers : {None, (list of) int}
            The walker index (or a list of indices) to retrieve. If None,
            stats from all walkers will be obtained.
        flatten : {True, bool}
            The returned array will be one dimensional, with all desired
            stats from all desired walkers concatenated together. If False,
            the returned array will have dimension requested temps x requested
            walkers x requested iterations.
        stats_group : {None, str}
            The group in `fp` from which to retrieve the stats. If
            None, searches in `fp.stats_group`.
        array_class : {None, array class}
            The type of array to return. The class must have a `from_kwargs`
            class method. If None, will return a FieldArray.

        Returns
        -------
        array_class
            The likelihood stats, as an instance of the given
            `array_class` (`FieldArray` if `array_class` is None).
        """
        if stats_group is None:
            stats_group = fp.stats_group
        if array_class is None:
            array_class = FieldArray
        fields = fp[stats_group].keys()
        return cls._read_fields(fp, stats_group, fields, array_class,
                                thin_start=thin_start,
                                thin_interval=thin_interval,
                                thin_end=thin_end, iteration=iteration,
                                temps=temps, walkers=walkers, flatten=flatten)

    @classmethod
    def compute_acls(cls, fp, start_index=None, end_index=None):
        """Computes the autocorrleation length for all variable args for all
        walkers for all temps in the given file. If the returned acl is inf,
        will default to the number of requested iterations.

        Parameters
        -----------
        fp : InferenceFile
            An open file handler to read the samples from.
        start_index : {None, int}
            The start index to compute the acl from. If None, will try to use
            the number of burn-in iterations in the file; otherwise, will start
            at the first sample.
        end_index : {None, int}
            The end index to compute the acl to. If None, will go to the end
            of the current iteration.

        Returns
        -------
        WaveformArray
            An ntemps x nwalkers `WaveformArray` containing the acl for each
            walker and temp for each variable argument, with the variable
            arguments as fields.
        """
        acls = {}
        if end_index is None:
            end_index = fp.niterations
        tidx = numpy.arange(fp.ntemps)
        widx = numpy.arange(fp.nwalkers)
        for param in fp.variable_args:
            these_acls = numpy.zeros((fp.ntemps, fp.nwalkers), dtype=int)
            for tk in tidx:
                for wi in widx:
                    samples = cls.read_samples(
                            fp, param,
                            thin_start=start_index, thin_interval=1,
                            thin_end=end_index,
                            walkers=wi, temps=tk)[param]
                    acl = autocorrelation.calculate_acl(samples)
                    these_acls[tk, wi] = int(min(acl, samples.size))
            acls[param] = these_acls
        return WaveformArray.from_kwargs(**acls)

    @staticmethod
    def write_acls(fp, acls):
        """Writes the given autocorrelation lengths to the given file. The acl
        of each walker at each temperature and each parameter is saved to
        `fp[fp.samples_group/{param}/temp{k}/walker{i}].attrs['acl']`; the
        maximum over all the walkers for a given temperature and param is
        saved to `fp[fp.samples_group/{param}/temp{k}].attrs['acl']`; the
        maximum over all of the temperatures and walkers is saved to
        `fp[fp.samples_group/{param}].attrs['acl']`; the maximum over all the
        parameters, temperatures, and walkers is saved to the file's 'acl'
        attribute.

        Parameters
        ----------
        fp : InferenceFile
            An open file handler to write the samples to.
        acls : WaveformArray
            An array of autocorrelation lengths (the sort of thing returned by
            `compute_acls`).

        Returns
        -------
        acl
            The maximum of the acls that was written to the file.
        """
        # write the individual acls
        pgroup = fp.samples_group + '/{param}'
        tgroup = pgroup + '/temp{tk}'
        group = tgroup + '/walker{wi}'
        tidx = numpy.arange(fp.ntemps)
        overall_max = 0
        for param in acls.fieldnames:
            max_acls = []
            for tk in tidx:
                max_acl = 0
                for wi, acl in enumerate(acls[param][tk, :]):
                    fp[group.format(param=param, tk=tk,
                                    wi=wi)].attrs['acl'] = acl
                    max_acl = max(max_acl, acl)
                # write the maximum over the walkers
                fp[tgroup.format(param=param, tk=tk)].attrs['acl'] = max_acl
                max_acls.append(max_acl)
            # write the maximum over the temperatures
            this_max = max(max_acls)
            fp[pgroup.format(param=param)].attrs['acl'] = this_max
            overall_max = max(overall_max, this_max)
            # write the maximum over the params
            fp[pgroup.format(param=param)].attrs['acl'] = max_acl
        # write the maximum over all params
        fp.attrs['acl'] = overall_max
        return fp.attrs['acl']

    @staticmethod
    def read_acls(fp):
        """Reads the acls of all the walker chains saved in the given file.

        Parameters
        ----------
        fp : InferenceFile
            An open file handler to read the acls from.

        Returns
        -------
        WaveformArray
            An ntemps x nwalkers `WaveformArray` containing the acls for
            every temp and walker, with the variable arguments as fields.
        """
        group = fp.samples_group + '/{param}/temp{tk}/walker{wi}'
        widx = numpy.arange(fp.nwalkers)
        tidx = numpy.arange(fp.ntemps)
        arrays = {}
        for param in fp.variable_args:
            arrays[param] = numpy.array([
                [fp[group.format(param=param, tk=tk, wi=wi)].attrs['acl']
                    for wi in widx]
                for tk in tidx])
        return WaveformArray.from_kwargs(**arrays)

    @classmethod
    def calculate_logevidence(cls, fp, thin_start=None, thin_end=None,
                              thin_interval=None):
        """Calculates the log evidence from the given file using emcee's
        thermodynamic integration.

        Parameters
        ----------
        fp : InferenceFile
            An open file handler to read the stats from.
        thin_start : int
            Index of the sample to begin returning stats. Default is to read
            stats after burn in. To start from the beginning set thin_start
            to 0.
        thin_interval : int
            Interval to accept every i-th sample. Default is to use the
            `fp.acl`. If `fp.acl` is not set, then use all stats
            (set thin_interval to 1).
        thin_end : int
            Index of the last sample to read. If not given then
            `fp.niterations` is used.

        Returns
        -------
        lnZ : float
            The estimate of log of the evidence.
        dlnZ : float
            The error on the estimate.
        """
        try:
            import emcee
        except ImportError:
            raise ImportError("emcee is not installed.")

        logstats = cls.read_likelihood_stats(fp, thin_start=thin_start,
                                             thin_end=thin_end,
                                             thin_interval=thin_interval,
                                             temps='all', flatten=False)
        # get the likelihoods
        logls = logstats['loglr'] + fp.lognl
        # we need the betas that were used
        betas = fp.attrs['betas']
        # annoyingly, theromdynaimc integration in PTSampler is an instance
        # method, so we'll implement a dummy one
        ntemps = fp.ntemps
        nwalkers = fp.nwalkers
        ndim = len(fp.variable_args)
        dummy_sampler = emcee.PTSampler(ntemps, nwalkers, ndim, None,
                                        None, betas=betas)
        return dummy_sampler.thermodynamic_integration_log_evidence(
            logls=logls, fburnin=0.)
