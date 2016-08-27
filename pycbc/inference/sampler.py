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

import numpy
from pycbc.io import WaveformArray, FieldArray
from pycbc.filter import autocorrelation

#
# =============================================================================
#
#                                   Helper functions
#
# =============================================================================
#
def get_slice(fp, thin_start=None, thin_interval=None, thin_end=None):
    """Formats a slice using the given arguments that can be used to retrieve
    a thinned array from an InferenceFile.

    Parameters
    ----------
    fp : InferenceFile
        An open inference file. The `burn_in_iterations` and acl will try to
        be obtained from the file's attributes.
    thin_start : {None, int}
        The starting index to use. If None, will try to retrieve the
        `burn_in_iterations` from the given file. If no `burn_in_iterations`
        exists, will default to the start of the array.
    thin_interval : {None, int}
        The interval to use. If None, will try to retrieve the acl from the
        given file. If no acl attribute exists, will default to 1.
    thin_end : {None, int}
        The end index to use. If None, will retrieve to the end of the array.

    Returns
    -------
    slice :
        The slice needed.
    """
    # default is to skip burn in samples
    if thin_start is None:
        try:
            thin_start = fp.burn_in_iterations
        except KeyError:
            pass
    # default is to use stored ACL and accept every i-th sample
    if thin_interval is None:
        try:
            thin_interval = int(numpy.ceil(fp.acl))
        except KeyError:
            pass
    return slice(thin_start, thin_end, thin_interval)


#
# =============================================================================
#
#                                   Samplers
#
# =============================================================================
#

class _BaseSampler(object):
    """Base container class for running the inference sampler that will
    generate the posterior distributions.

    Parameters
    ----------
    likelihood_evaluator : LikelihoodEvaluator
        An instance of a pycbc.inference.likelihood evaluator.
    """
    name = None

    def __init__(self, likelihood_evaluator):
        self.likelihood_evaluator = likelihood_evaluator

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
        cls
            A sampler initialized based on the given arguments.
        """
        raise NotImplementedError("from_cli function not set")

    @property
    def ifos(self):
        """Returns the ifos that were sampled over."""
        return self.likelihood_evaluator.waveform_generator.detector_names

    @property
    def variable_args(self):
        """Returns the variable args used by the likelihood evaluator.
        """
        return self.likelihood_evaluator.waveform_generator.variable_args

    @property
    def chain(self):
        """This function should return the past samples as a
        [additional dimensions x] niterations x ndim array, where ndim are the
        number of variable args, niterations the number of iterations, and
        additional dimeionions are any additional dimensions used by the
        sampler (e.g, walkers, temperatures).
        """
        return NotImplementedError("chain function not set.")

    @property
    def niterations(self):
        """Get the current number of iterations."""
        return self.chain.shape[-2]

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
        return NotImplementedError("lnpost function not set.")

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

    def run(self, niterations):
        """This function should run the sampler.
        """
        raise NotImplementedError("run function not set.")

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
        fp.attrs['variable_args'] = self.variable_args
        fp.attrs["niterations"] = self.niterations
        fp.attrs["lognl"] = self.likelihood_evaluator.lognl
        sargs = self.likelihood_evaluator.waveform_generator.static_args
        fp.attrs["static_args"] = sargs.keys()
        for arg,val in sargs.items():
            fp.attrs[arg] = val

class _BaseMCMCSampler(_BaseSampler):
    """This class is used to construct the MCMC sampler from the kombine-like
    packages.

    Parameters
    ----------
    sampler : sampler instance
        An instance of an MCMC sampler similar to kombine or emcee.
    likelihood_evaluator : likelihood class
        An instance of the likelihood class from the
        pycbc.inference.likelihood module.
    min_burn_in : {None, int}
        Set the minimum number of burn in iterations to use. If None,
        `burn_in_iterations` will be initialized to `0`.

    Attributes
    ----------
    sampler :
        The MCMC sampler instance used.
    p0 : nwalkers x ndim array
        The initial position of the walkers. Set by using set_p0. If not set
        yet, a ValueError is raised when the attribute is accessed.
    pos : {None, array}
        An array of the current walker positions.
    """
    name = None
    def __init__(self, sampler, likelihood_evaluator, min_burn_in=None):
        self._sampler = sampler 
        self._pos = None
        self._p0 = None
        if min_burn_in is None:
            min_burn_in = 0
        self.burn_in_iterations = min_burn_in
        # initialize
        super(_BaseMCMCSampler, self).__init__(likelihood_evaluator)

    @property
    def sampler(self):
        return self._sampler

    @property
    def pos(self):
        return self._pos

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
            An nwalkers x ndim array of the initial positions that were set.
        """
        # loop over all walkers and then parameters
        # find the distribution that has that parameter in it and draw a
        # random value from the distribution
        nwalkers = self.nwalkers
        ndim = len(self.variable_args)
        pmap = dict([[param,k] for k,param in enumerate(self.variable_args)])
        p0 = numpy.ones((nwalkers, ndim))
        for dist in prior_distributions:
            ps = dist.rvs(size=nwalkers)
            for param in dist.params:
                p0[:,pmap[param]] = ps[param]
        self._p0 = p0
        return p0

    @property
    def p0(self):
        if self._p0 is None:
            raise ValueError("initial positions not set; run set_p0")
        return self._p0

    @property
    def nwalkers(self):
        """Get the number of walkers."""
        return self.chain.shape[-3]

    @property
    def acceptance_fraction(self):
        """Get the fraction of walkers that accepted each step as an array.
        """
        return self._sampler.acceptance_fraction

    @property
    def likelihood_stats(self):
        """Returns the likelihood stats as a FieldArray, with field names
        corresponding to the type of data returned by the likelihood evaluator.
        The returned array has shape nwalkers x niterations. If no additional
        stats were returned to the sampler by the likelihood evaluator, returns
        None.
        """
        stats = numpy.array(self._sampler.blobs)
        if stats.size == 0:
            return None
        arrays = dict([[field, stats[:,:,fi]]
                    for fi,field in
                        enumerate(self.likelihood_evaluator.metadata_fields)])
        return FieldArray.from_kwargs(**arrays).transpose()

    # write and read functions
    def write_metadata(self, fp):
        """Writes metadata about this sampler to the given file. Metadata is
        written to the file's `attrs`.

        Parameters
        ----------
        fp : InferenceFile
            A file handler to an open inference file.
        """
        super(_BaseMCMCSampler, self).write_metadata(fp)
        # add info about walkers, burn in
        fp.attrs["nwalkers"] = self.nwalkers
        fp.attrs['burn_in_iterations'] = self.burn_in_iterations

    def write_chain(self, fp, max_iterations=None):
        """Writes the samples from the current chain to the given file. Results
        are written to: `fp[fp.samples_group/{vararg}/walker{i}]`, where
        `{vararg}` is the name of a variable arg, and `{i}` is the index of
        a walker.

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
        # chain is nwalkers x niterations x ndim
        samples = self.chain
        nwalkers, niterations, _ = samples.shape

        group = fp.samples_group + '/{name}/walker{wi}'

        # create an empty array if desired, in case this is the first time
        # writing
        if max_iterations is not None:
            if max_iterations < niterations:
                raise IndexError("The provided max size is less than the "
                    "number of iterations")
            out = numpy.zeros(max_iterations, dtype=samples.dtype)

        # loop over number of dimensions
        widx = numpy.arange(nwalkers)
        for pi,param in enumerate(self.variable_args):
            # loop over number of walkers
            for wi in widx:
                dataset_name = group.format(name=param, wi=wi)
                try:
                    fp[dataset_name][:niterations] = samples[wi,:,pi]
                except KeyError:
                    # dataset doesn't exist yet, see if a larger array is
                    # desired
                    if max_iterations is not None:
                        out[:niterations] = samples[wi,:,pi]
                        fp[dataset_name] = out
                    else:
                        fp[dataset_name] = samples[wi,:,pi]

    def write_lnpost(self, fp, max_iterations=None):
        """Writes the `lnpost`s to the given file. Results are written to:
        `fp[fp.stats_group/lnpost/walker{i}]`, where `{i}` is the index of
        a walker.

        Parameters
        -----------
        fp : InferenceFile
            A file handler to an open inference file.
        max_iterations : {None, int}
            If `lnpost`s have not previously been written to the file, a new
            dataset will be created. By default, the size of this dataset will
            be whatever the length of the sampler's chain is at this point. If
            you intend to run more iterations, set this value to that size so
            that the array in the file will be large enough to accomodate
            future data.
        """
        # lnposts are an nwalkers x niterations array
        lnposts = self.lnpost
        nwalkers, niterations = lnposts.shape

        group = fp.stats_group + '/lnpost/walker{wi}'

        # create an empty array if desired, in case this is the first time
        # writing
        if max_iterations is not None:
            if max_iterations < niterations:
                raise IndexError("The provided max size is less than the "
                    "number of iterations")
            out = numpy.zeros(max_iterations, dtype=lnposts.dtype)

        # loop over number of walkers
        for wi in range(nwalkers):
            dataset_name = group.format(wi=wi)
            try:
                fp[dataset_name][:niterations] = lnposts[wi,:]
            except KeyError:
                # dataset doesn't exist yet, see if a larger array is
                # desired
                if max_iterations is not None:
                    out[:niterations] = lnposts[wi,:]
                    fp[dataset_name] = out
                else:
                    fp[dataset_name] = lnposts[wi,:]

    def write_likelihood_stats(self, fp, max_iterations=None):
        """Writes the `likelihood_stats` to the given file.  Results are
        written to: `fp[fp.stats_group/{field}/walker{i}]`, where `{i}` is
        the index of a walker and `{field}` is the name of each field returned
        by `likelihood_stats`. If nothing is returned by `likelihood_stats`,
        this does nothing.

        Parameters
        -----------
        fp : InferenceFile
            A file handler to an open inference file.
        max_iterations : {None, int}
            If the stats have not previously been written to the file, a new
            dataset will be created. By default, the size of this dataset will
            be whatever the length of the sampler's chain is at this point. If
            you intend to run more iterations, set this value to that size so
            that the array in the file will be large enough to accomodate
            future data.

        Returns
        -------
        stats : {FieldArray, None}
            The stats that were written, as a FieldArray. If there were no
            stats, returns None.
        """
        # stats is an nwalkers x niterations array
        stats = self.likelihood_stats
        if stats is None:
            return None
        nwalkers, niterations = stats.shape
        fields = stats.fieldnames

        group = fp.stats_group + '/{param}/walker{wi}'

        if max_iterations is not None and max_iterations < niterations:
            raise IndexError("The provided max size is less than the "
                "number of iterations")

        for param in fields:
            # create an empty array if desired, in case this is the first time
            # writing
            out = numpy.zeros(max_iterations, dtype=stats.dtype[param])
            # loop over number of walkers
            for wi in range(nwalkers):
                dataset_name = group.format(param=param, wi=wi)
                try:
                    fp[dataset_name][:niterations] = stats[param][wi,:]
                except KeyError:
                    # dataset doesn't exist yet, see if a larger array is
                    # desired
                    if max_iterations is not None:
                        out[:niterations] = stats[param][wi,:]
                        fp[dataset_name] = out
                    else:
                        fp[dataset_name] = stats[param][wi,:]
        return stats

    def write_acceptance_fraction(self, fp, max_iterations=None):
        """Write acceptance_fraction data to file. Results are written to
        `fp[acceptance_fraction]`.

        Parameters
        -----------
        fp : InferenceFile
            A file handler to an open inference file.
        max_iterations : {None, int}
            If acceptance fraction have not previously been written to the
            file, a new dataset will be created. By default, the size of this
            dataset will be whatever the length of the sampler's chain is at
            this point. If you intend to run more iterations, set this value
            to that size so that arrays in the file will be large enough to
            accomodate future data.
        """
        dataset_name = "acceptance_fraction"
        acf = self.acceptance_fraction
        try:
            fp[dataset_name][:acf.size] = acf
        except KeyError:
            # dataset doesn't exist yet, see if a larger array is
            # desired
            if max_iterations is not None:
                out = numpy.zeros(max_iterations, dtype=acf.dtype)
                out[:acf.size] = acf
            else:
                out = acf
            fp[dataset_name] = out

    def write_results(self, fp, max_iterations=None):
        """Writes metadata, samples, likelihood stats, and acceptance fraction
        to the given file. Also computes and writes the autocorrleation lengths
        of the chains. See the various write function for details.

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
        self.write_acceptance_fraction(fp, max_iterations=max_iterations)

    @classmethod
    def _read_fields(cls, fp, fields_group, fields, array_class,
            thin_start=None, thin_interval=None, thin_end=None, iteration=None,
            walkers=None, flatten=True):
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

        # get the slice to use
        if iteration is not None:
            get_index = iteration
        else:
            if thin_end is None:
                # use the number of current iterations
                thin_end = fp.niterations
            get_index = get_slice(fp, thin_start=thin_start, thin_end=thin_end,
                thin_interval=thin_interval)

        # load
        arrays = {}
        group = fields_group + '/{name}/walker{wi}'
        for name in fields:
            these_arrays = [
                    fp[group.format(name=name, wi=wi)][get_index]
                    for wi in walkers]
            if flatten:
                arrays[name] = numpy.hstack(these_arrays)
            else:
                arrays[name] = numpy.vstack(these_arrays)
        return array_class.from_kwargs(**arrays)

    @classmethod
    def read_samples(cls, fp, parameters,
            thin_start=None, thin_interval=None, thin_end=None, iteration=None,
            walkers=None, flatten=True, samples_group=None,
            array_class=None):
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
        flatten : {True, bool}
            The returned array will be one dimensional, with all desired
            samples from all desired walkers concatenated together. If False,
            the returned array will have dimension requested walkers
            x requested iterations.
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
        loadfields = array_class.parse_parameters(parameters,
            possible_fields=possible_fields)
        return cls._read_fields(fp, samples_group, loadfields, array_class,
                thin_start=thin_start, thin_interval=thin_interval,
                thin_end=thin_end, iteration=iteration, walkers=walkers,
                flatten=flatten)

    @classmethod
    def read_likelihood_stats(cls, fp,
            thin_start=None, thin_interval=None, thin_end=None, iteration=None,
            walkers=None, flatten=True, stats_group=None, array_class=None):
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
        walkers : {None, (list of) int}
            The walker index (or a list of indices) to retrieve. If None,
            stats from all walkers will be obtained.
        flatten : {True, bool}
            The returned array will be one dimensional, with all desired
            stats from all desired walkers concatenated together. If False,
            the returned array will have dimension requested walkers
            x requested iterations.
        stats_group : {None, str}
            The group in `fp` from which to retrieve the stats. If
            None, searches in `fp.stats_group`.
        array_class : {None, array class}
            The type of array to return. The class must have a `from_kwargs`
            class method. If None, will return a FieldArray.

        Returns
        -------
        array_class
            The likelihood stats, as an instance of a the given
            `array_class` (`FieldArray` if `array_class` is None).
        """
        if stats_group is None:
            stats_group = fp.stats_group
        if array_class is None:
            array_class = FieldArray
        fields = fp[stats_group].keys()
        return cls._read_fields(fp, stats_group, fields, array_class,
                thin_start=thin_start, thin_interval=thin_interval,
                thin_end=thin_end, iteration=iteration, walkers=walkers,
                flatten=flatten)

    @staticmethod
    def read_acceptance_fraction(fp,
            thin_start=None, thin_interval=None, thin_end=None,
            iteration=None):
        """Reads the acceptance fraction from the given file.

        Parameters
        -----------
        fp : InferenceFile
            An open file handler to read the samples from.
        walkers : {None, (list of) int}
            The walker index (or a list of indices) to retrieve. If None,
            samples from all walkers will be obtained.
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

        Returns
        -------
        array
            Array of acceptance fractions with shape (requested iterations,).
        """
        # get the slice to use
        if iteration is not None:
            get_index = iteration
        else:
            if thin_end is None:
                # use the number of current iterations
                thin_end = fp.niterations
            get_index = get_slice(fp, thin_start=thin_start, thin_end=thin_end,
                thin_interval=thin_interval)
        acfs = fp['acceptance_fraction'][get_index]
        if iteration is not None:
            acfs = numpy.array([acfs])
        return acfs

    @classmethod
    def compute_acls(cls, fp, start_index=None, end_index=None):
        """Computes the autocorrleation length for all variable args and all
        walkers in the given file. If the returned acl is inf, will default
        to the number of requested iterations.

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
            An nwalkers-long `WaveformArray` containing the acl for each walker
            and each variable argument, with the variable arguments as fields.
        """
        acls = {}
        if end_index is None:
            end_index = fp.niterations
        widx = numpy.arange(fp.nwalkers)
        for param in fp.variable_args:
            these_acls = []
            for wi in widx:
                samples = cls.read_samples(fp, param,
                        thin_start=start_index, thin_interval=1,
                        thin_end=end_index,
                        walkers=wi)[param]
                acl = autocorrelation.calculate_acl(samples)
                these_acls.append(min(acl, samples.size))
            acls[param] = numpy.array(these_acls, dtype=int)
        return WaveformArray.from_kwargs(**acls)

    @staticmethod
    def write_acls(fp, acls):
        """Writes the given autocorrelation lengths to the given file. The acl
        of each walker and each parameter is saved to
        `fp[fp.samples_group/{param}/walker{i}].attrs['acl']`; the maximum
        over all the walkers for a given param is saved to
        `fp[fp.samples_group/{param}].attrs['acl']`; the maximum over all the
        parameters and all of the walkers is saved to the file's 'acl'
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
        group = pgroup + '/walker{wi}'
        max_acls = []
        for param in acls.fieldnames:
            max_acl = 0
            for wi,acl in enumerate(acls[param]): 
                fp[group.format(param=param, wi=wi)].attrs['acl'] = acl
                max_acl = max(max_acl, acl)
            max_acls.append(max_acl)
            # write the maximum over the params
            fp[pgroup.format(param=param)].attrs['acl'] = max_acl
        # write the maximum over all params
        fp.attrs['acl'] = max(max_acls)
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
            An nwalkers-long `WaveformArray` containing the acl for each walker
            and each variable argument, with the variable arguments as fields.
        """
        group = fp.samples_group + '/{param}/walker{wi}'
        widx = numpy.arange(fp.nwalkers)
        arrays = {}
        for param in fp.variable_args:
            arrays[param] = numpy.array([
                fp[group.format(param=param, wi=wi)].attrs['acl']
                for wi in widx])
        return WaveformArray.from_kwargs(**arrays)


class KombineSampler(_BaseMCMCSampler):
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
        if self.burn_in_iterations == 0:
            # no burn in, use the initial positions
            p0 = self.p0
        else:
            p0 = None
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
        return self._sampler.chain.transpose((1,0,2))

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
        res = self._sampler.burnin(self.p0)
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
        super(EmceeEnsembleSampler, self).__init__(sampler,
            likelihood_evaluator, min_burn_in=burn_in_iterations)

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
            raise ValueError("%s requires that you provide a non-zero "%(
                cls.name) + "--min-burn-in if not skipping burn-in")
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

samplers = {
    KombineSampler.name : KombineSampler,
    EmceeEnsembleSampler.name : EmceeEnsembleSampler,
}
