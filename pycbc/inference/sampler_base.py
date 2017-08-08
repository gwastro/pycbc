# Copyright (C) 2016  Christopher M. Biwer, Collin Capano
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
from pycbc.io import FieldArray
from pycbc.filter import autocorrelation

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
        self._lastclear = 0

    @classmethod
    def from_cli(cls, opts, likelihood_evaluator, pool=None, likelihood_call=None):
        """This function create an instance of this sampler from the given
        command-line options.
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
        additional dimeionions are any additional dimensions used by the
        sampler (e.g, walkers, temperatures).
        """
        return NotImplementedError("chain function not set.")

    @property
    def samples(self):
        """This function should return the past samples as a [additional
        dimensions x] niterations field array, where the fields are union
        of the sampling args and the variable args.
        """
        return NotImplementedError("samples function not set.")
 
    @property
    def clear_chain(self):
        """This function should clear the current chain of samples from memory.
        """
        return NotImplementedError("clear chain function not set.")

    @property
    def niterations(self):
        """Get the current number of iterations."""
        return self.chain.shape[-2] + self._lastclear

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

class BaseMCMCSampler(_BaseSampler):
    """This class is used to construct the MCMC sampler from the kombine-like
    packages.

    Parameters
    ----------
    sampler : sampler instance
        An instance of an MCMC sampler similar to kombine or emcee.
    likelihood_evaluator : likelihood class
        An instance of the likelihood class from the
        pycbc.inference.likelihood module.

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

    def __init__(self, sampler, likelihood_evaluator):
        self._sampler = sampler
        self._pos = None
        self._p0 = None
        self._currentblob = None
        self._nwalkers = None
        self._lastclear = 0
        self.burn_in_iterations = None
        # initialize
        super(BaseMCMCSampler, self).__init__(likelihood_evaluator)

    @property
    def sampler(self):
        return self._sampler

    @property
    def pos(self):
        return self._pos

    def set_p0(self, samples=None, prior=None):
        """Sets the initial position of the walkers. Must be supplied a
        FieldArray of values or a list of Distributions or a PriorEvaluator.

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
            An nwalkers x ndim array of the initial positions that were set.
        """
        # create a (nwalker, ndim) array for initial positions
        nwalkers = self.nwalkers
        ndim = len(self.variable_args)
        p0 = numpy.ones((nwalkers, ndim))
        # if samples are given then use those as initial positions
        if samples is not None:
            # transform to sampling parameter space
            samples = self.likelihood_evaluator.apply_sampling_transforms(
                samples)
        # draw random samples if samples are not provided
        else:
            samples = self.likelihood_evaluator.prior_rvs(size=nwalkers,
                                                          prior=prior)
        # convert to 2D array
        for i, param in enumerate(self.sampling_args):
            p0[:, i] = samples[param]
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
        return self._nwalkers

    @property
    def acceptance_fraction(self):
        """Get the fraction of walkers that accepted each step as an array.
        """
        return self._sampler.acceptance_fraction

    @property
    def samples(self):
        """Returns the samples in the chain as a FieldArray.

        If the sampling args are not the same as the variable args, the
        returned samples will have both the sampling and the variable args.

        The returned FieldArray has dimension [additional dimensions x]
        nwalkers x niterations.
        """
        # chain is a [additional dimensions x] niterations x ndim array
        samples = self.chain
        sampling_args = self.sampling_args
        # convert to dictionary to apply boundary conditions
        samples = {param: samples[...,ii]
                   for ii,param in enumerate(sampling_args)}
        samples = self.likelihood_evaluator._prior.apply_boundary_conditions(
            **samples)
        # now convert to field array
        samples = FieldArray.from_arrays([samples[param]
                                          for param in sampling_args],
                                         names=sampling_args)
        # apply transforms to go to variable args space
        return self.likelihood_evaluator.apply_sampling_transforms(samples,
            inverse=True)

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
        # we'll force arrays to float; this way, if there are `None`s in the
        # blobs, they will be changed to `nan`s
        arrays = {field: stats[..., fi].astype(float)
                  for fi, field in
                  enumerate(self.likelihood_evaluator.metadata_fields)}
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
        super(BaseMCMCSampler, self).write_metadata(fp)
        # add info about walkers, burn in
        fp.attrs["nwalkers"] = self.nwalkers

    @staticmethod
    def write_samples_group(fp, samples_group, parameters, samples,
                            start_iteration=0, end_iteration=None,
                            index_offset=0, max_iterations=None):
        """Writes samples to the given file.

        Results are written to:
        
            `fp[samples_group/{vararg}/walker{i}]`,
            
        where `{vararg}` is the name of a variable arg, and `{i}` is the
        index of a walker.

        Parameters
        -----------
        fp : InferenceFile
            A file handler to an open inference file.
        samples_group : str
            Name of samples group to write.
        parameters : list
            The parameters to write to the file.
        samples : FieldArray
            The samples to write. Should be a FieldArray with fields containing
            the samples to write and shape nwalkers x niterations.
        start_iteration : {0, int}
            Write results starting from the given iteration.
        end_iteration : {None, int}
            Write results up to the given iteration.
        index_offset : int, optional
            Write the samples to the arrays on disk starting at
            `start_iteration` + `index_offset`. For example, if
            `start_iteration=0`, `end_iteration=1000` and `index_offset=500`,
            then `samples[0:1000]` will be written to indices `500:1500` in the
            arrays on disk. This is needed if you are adding new samples to
            a chain that was previously written to file, and you want to
            preserve the history (e.g., after a checkpoint). Default is 0.
        max_iterations : int, optional
            Set the maximum size that the arrays in the hdf file may be resized
            to. Only applies if the samples have not previously been written
            to file. The default (None) is to use the maximum size allowed by
            h5py.
        """
        # due to clearing memory, there can be a difference between indices in
        # memory and on disk
        nwalkers, niterations = samples.shape
        niterations += index_offset
        fa = start_iteration # file start index
        if end_iteration is None:
            end_iteration = niterations
        fb = end_iteration # file end index
        ma = fa - index_offset # memory start index
        mb = fb - index_offset # memory end index

        if max_iterations is not None and max_iterations < niterations:
            raise IndexError("The provided max size is less than the "
                             "number of iterations")

        group = samples_group + '/{name}/walker{wi}'

        # loop over number of dimensions
        widx = numpy.arange(nwalkers)
        for param in parameters:
            # loop over number of walkers
            for wi in widx:
                dataset_name = group.format(name=param, wi=wi)
                try:
                    if fb > fp[dataset_name].size:
                        # resize the dataset
                        fp[dataset_name].resize(fb, axis=0)
                    fp[dataset_name][fa:fb] = samples[param][wi, ma:mb]
                except KeyError:
                    # dataset doesn't exist yet
                    fp.create_dataset(dataset_name, (fb,),
                                      maxshape=(max_iterations,),
                                      dtype=samples[param].dtype)
                    fp[dataset_name][fa:fb] = samples[param][wi, ma:mb]

    def write_chain(self, fp, start_iteration=0, end_iteration=None,
                    max_iterations=None):
        """Writes the samples from the current chain to the given file.
        
        Results are written to:
        
            `fp[fp.samples_group/{field}/(temp{k}/)walker{i}]`,
        
        where `{i}` is the index of a walker, `{field}` is the name of each
        field returned by `likelihood_stats`, and, if the sampler is
        multitempered, `{k}` is the temperature.

        Parameters
        -----------
        fp : InferenceFile
            A file handler to an open inference file.
        start_iteration : {0, int}
            Write results starting from the given iteration.
        end_iteration : {None, int}
            Write results up to the given iteration.
        max_iterations : int, optional
            Set the maximum size that the arrays in the hdf file may be resized
            to. Only applies if the samples have not previously been written
            to file. The default (None) is to use the maximum size allowed by
            h5py.
        samples_group : str
            Name of samples group to write.
        """
        # samples is a nwalkers x niterations field array
        samples = self.samples
        parameters = self.variable_args
        samples_group = fp.samples_group
        # write data
        self.write_samples_group(
                         fp, samples_group, parameters, samples,
                         start_iteration=start_iteration,
                         end_iteration=end_iteration,
                         index_offset=self._lastclear,
                         max_iterations=max_iterations)

    def write_likelihood_stats(self, fp, start_iteration=0, end_iteration=None,
                               max_iterations=None):
        """Writes the `likelihood_stats` to the given file.
        
        Results are written to:
        
            `fp[fp.stats_group/{field}/(temp{k}/)walker{i}]`,
        
        where `{i}` is the index of a walker, `{field}` is the name of each
        field returned by `likelihood_stats`, and, if the sampler is
        multitempered, `{k}` is the temperature.  If nothing is returned by
        `likelihood_stats`, this does nothing.

        Parameters
        -----------
        fp : InferenceFile
            A file handler to an open inference file.
        start_iteration : {0, int}
            Write results starting from the given iteration.
        end_iteration : {None, int}
            Write results up to the given iteration.
        max_iterations : int, optional
            Set the maximum size that the arrays in the hdf file may be resized
            to. Only applies if the samples have not previously been written
            to file. The default (None) is to use the maximum size allowed by
            h5py.

        Returns
        -------
        stats : {FieldArray, None}
            The stats that were written, as a FieldArray. If there were no
            stats, returns None.
        """
        samples = self.likelihood_stats
        if samples is None:
            return None
        # ensure the prior is in the variable args parameter space
        if 'logjacobian' in samples.fieldnames:
            samples['prior'] -= samples['logjacobian']
        parameters = samples.fieldnames
        samples_group = fp.stats_group
        # write data
        self.write_samples_group(
                         fp, samples_group, parameters, samples,
                         start_iteration=start_iteration,
                         end_iteration=end_iteration,
                         index_offset=self._lastclear,
                         max_iterations=max_iterations)
        return samples

    def write_acceptance_fraction(self, fp, start_iteration=0,
                                  end_iteration=None, max_iterations=None):
        """Write acceptance_fraction data to file. Results are written to
        `fp[acceptance_fraction]`.

        Parameters
        -----------
        fp : InferenceFile
            A file handler to an open inference file.
        max_iterations : int, optional
            Set the maximum size that the arrays in the hdf file may be resized
            to. Only applies if the acceptance fraction has not previously been
            written to the file. The default (None) is to use the maximum size
            allowed by h5py.
        """
        dataset_name = "acceptance_fraction"
        acf = self.acceptance_fraction

        if end_iteration is None:
            end_iteration = acf.size

        if max_iterations is not None and max_iterations < acf.size:
            raise IndexError("The provided max size is less than the "
                             "number of iterations")

        try:
            if end_iteration > fp[dataset_name].size:
                # resize the dataset
                fp[dataset_name].resize(end_iteration, axis=0)
            fp[dataset_name][start_iteration:end_iteration] = \
                acf[start_iteration:end_iteration]
        except KeyError:
            # dataset doesn't exist yet
            fp.create_dataset(dataset_name, (end_iteration,),
                              maxshape=(max_iterations,),
                              dtype=acf.dtype)
            fp[dataset_name][start_iteration:end_iteration] = \
                acf[start_iteration:end_iteration]

    def write_results(self, fp, start_iteration=0, end_iteration=None,
                      max_iterations=None):
        """Writes metadata, samples, likelihood stats, and acceptance fraction
        to the given file. Also computes and writes the autocorrleation lengths
        of the chains. See the various write function for details.

        Parameters
        -----------
        fp : InferenceFile
            A file handler to an open inference file.
        start_iteration : {0, int}
            Write results starting from the given iteration.
        end_iteration : {None, int}
            Write results up to the given iteration.
        max_iterations : int, optional
            Set the maximum size that the arrays in the hdf file may be resized
            to. Only applies if the acceptance fraction has not previously been
            written to the file. The default (None) is to use the maximum size
            allowed by h5py.
        """
        self.write_metadata(fp)
        self.write_chain(fp, start_iteration=start_iteration,
                         end_iteration=end_iteration,
                         max_iterations=max_iterations)
        self.write_likelihood_stats(fp, start_iteration=start_iteration,
                                    end_iteration=end_iteration,
                                    max_iterations=max_iterations)
        self.write_acceptance_fraction(fp, start_iteration=start_iteration,
                                    end_iteration=end_iteration,
                                    max_iterations=max_iterations)


    @staticmethod
    def _read_fields(fp, fields_group, fields, array_class,
                     thin_start=None, thin_interval=None, thin_end=None,
                     iteration=None, walkers=None, flatten=True):
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
            get_index = fp.get_slice(thin_start=thin_start, thin_end=thin_end,
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
                     thin_start=None, thin_interval=None, thin_end=None,
                     iteration=None, walkers=None, flatten=True,
                     samples_group=None, array_class=None):
        """Reads samples for the given parameter(s).

        Parameters
        -----------
        fp : InferenceFile
            An open file handler to read the samples from.
        parameters : (list of) strings
            The parameter(s) to retrieve. A parameter can be the name of any
            field in `fp[fp.samples_group]`, a virtual field or method of
            `FieldArray` (as long as the file contains the necessary fields
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
            a FieldArray.

        Returns
        -------
        array_class
            Samples for the given parameters, as an instance of a the given
            `array_class` (`FieldArray` if `array_class` is None).
        """
        # get the group to load from
        if samples_group is None:
            samples_group = fp.samples_group
        # get the type of array class to use
        if array_class is None:
            array_class = FieldArray
        # get the names of fields needed for the given parameters
        possible_fields = fp[samples_group].keys()
        loadfields = array_class.parse_parameters(parameters, possible_fields)
        return cls._read_fields(fp, samples_group, loadfields, array_class,
                                thin_start=thin_start,
                                thin_interval=thin_interval, thin_end=thin_end,
                                iteration=iteration, walkers=walkers,
                                flatten=flatten)

    @classmethod
    def n_independent_samples(cls, fp):
        """Returns the number of independent samples stored in a file.

        Parameters
        -----------
        fp : InferenceFile
            An open file handler to read.

        Returns
        -------
        int
            The number of independent samples.
        """
        # we'll just read a single parameter from the file
        samples = cls.read_samples(fp, fp.variable_args[0])
        return samples.size

    @staticmethod
    def read_acceptance_fraction(fp, thin_start=None, thin_interval=None,
                                 thin_end=None, iteration=None):
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
            get_index = fp.get_slice(thin_start=thin_start, thin_end=thin_end,
                                     thin_interval=thin_interval)
        acfs = fp['acceptance_fraction'][get_index]
        if iteration is not None:
            acfs = numpy.array([acfs])
        return acfs

    @classmethod
    def compute_acfs(cls, fp, start_index=None, end_index=None,
                     per_walker=False, walkers=None, parameters=None):
        """Computes the autocorrleation function of the variable args in the
        given file.

        By default, parameter values are averaged over all walkers at each
        iteration. The ACF is then calculated over the averaged chain. An
        ACF per-walker will be returned instead if ``per_walker=True``.

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
        per_walker : optional, bool
            Return the ACF for each walker separately. Default is False.
        walkers : optional, int or array
            Calculate the ACF using only the given walkers. If None (the
            default) all walkers will be used.
        parameters : optional, str or array
            Calculate the ACF for only the given parameters. If None (the
            default) will calculate the ACF for all of the variable args.

        Returns
        -------
        FieldArray
            A ``FieldArray`` of the ACF vs iteration for each parameter. If
            `per-walker` is True, the FieldArray will have shape
            ``nwalkers x niterations``.
        """
        acfs = {}
        if parameters is None:
            parameters = fp.variable_args
        if isinstance(parameters, str) or isinstance(parameters, unicode):
            parameters = [parameters]
        for param in parameters:
            if per_walker:
                # just call myself with a single walker
                if walkers is None:
                    walkers = numpy.arange(fp.nwalkers)
                arrays = [cls.compute_acfs(fp, start_index=start_index,
                                           end_index=end_index,
                                           per_walker=False, walkers=ii,
                                           parameters=param)[param]
                          for ii in walkers]
                acfs[param] = numpy.vstack(arrays)
            else:
                samples = cls.read_samples(fp, param,
                                           thin_start=start_index,
                                           thin_interval=1, thin_end=end_index,
                                           walkers=walkers,
                                           flatten=False)[param]
                samples = samples.mean(axis=0)
                acfs[param] = autocorrelation.calculate_acf(samples).numpy()
        return FieldArray.from_kwargs(**acfs)

    @classmethod
    def compute_acls(cls, fp, start_index=None, end_index=None):
        """Computes the autocorrleation length for all variable args in the
        given file.

        Parameter values are averaged over all walkers at each iteration.
        The ACL is then calculated over the averaged chain. If the returned ACL
        is `inf`,  will default to the number of current iterations.

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
        dict
            A dictionary giving the ACL for each parameter.
        """
        acls = {}
        for param in fp.variable_args:
            samples = cls.read_samples(fp, param,
                                           thin_start=start_index,
                                           thin_interval=1, thin_end=end_index,
                                           flatten=False)[param]
            samples = samples.mean(axis=0)
            acl = autocorrelation.calculate_acl(samples)
            if numpy.isinf(acl):
                acl = samples.size
            acls[param] = acl
        return acls

    @staticmethod
    def write_acls(fp, acls):
        """Writes the given autocorrelation lengths to the given file.
        
        The ACL of each parameter is saved to
        `fp[fp.samples_group/{param}].attrs['acl']`; the maximum over all the
        parameters is saved to the file's 'acl' attribute.

        Parameters
        ----------
        fp : InferenceFile
            An open file handler to write the samples to.
        acls : dict
            A dictionary of ACLs keyed by the parameter.

        Returns
        -------
        ACL
            The maximum of the acls that was written to the file.
        """
        # write the individual acls
        pgroup = fp.samples_group + '/{param}'
        for param in acls:
            fp[pgroup.format(param=param)].attrs['acl'] = acls[param]
        # write the maximum over all params
        fp.attrs['acl'] = max(acls.values())
        return fp.attrs['acl']

    @staticmethod
    def read_acls(fp):
        """Reads the acls of all the parameters in the given file.

        Parameters
        ----------
        fp : InferenceFile
            An open file handler to read the acls from.

        Returns
        -------
        dict
            A dictionary of the ACLs, keyed by the parameter name.
        """
        group = fp.samples_group + '/{param}'
        return {param: fp[group.format(param=param)]
                for param in fp.variable_args}

    def write_state(self, fp):
        """ Saves the state of the sampler in a file.
        """
        raise NotImplementedError("Writing state to file not implemented.")

    def set_state_from_file(self, fp):
        """ Sets the state of the sampler back to the instance saved in a file.
        """
        raise NotImplementedError("Reading state from file not implemented.")
