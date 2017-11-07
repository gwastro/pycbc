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
from pycbc.inference.sampler_base import _BaseSampler
from pycbc.io import FieldArray

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
    """
    name = "mcmc"

    def __init__(self, likelihood_evaluator,verbose=False):
        self.likelihood_evaluator = likelihood_evaluator
        self.verbose = verbose
        self._lastclear = 0
        self.last_sample = []
        self.samples_chain = []

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
        return cls(likelihood_evaluator,verbose=opts.verbose)

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
        # Need to convert self.likelihood_evaluator.sampling_args to a list
        # as it is sometimes a tuple
        return [s for s in self.likelihood_evaluator.sampling_args]

    @property
    def chain(self):
        """This function should return the past samples as a
        [additional dimensions x] niterations x ndim array, where ndim are the
        number of variable args, niterations the number of iterations, and
        additional dimensions are any additional dimensions used by the
        sampler (e.g, walkers, temperatures).
        """
        return self.samples_chain

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
        samples = {param: samples[param] for param in sampling_args}
        samples = self.likelihood_evaluator._prior.apply_boundary_conditions(
            **samples)
        # now convert to field array
        samples = FieldArray.from_arrays([samples[param]
                                          for param in sampling_args],
                                         names=sampling_args)
        # apply transforms to go to variable args space
        return self.likelihood_evaluator.apply_sampling_transforms(samples,
            inverse=True)

    def clear_chain(self):
        """This function should clear the current chain of samples from memory.
        """
        # store the iteration that the clear is occuring on
        self._lastclear = self.niterations
        self.last_sample = self.samples_chain[-1]
        self.samples_chain = []

    @property
    def niterations(self):
        """Get the current number of iterations."""
        return len(self.samples_chain)+self._lastclear

    @property
    def lnpost(self):
        """This function should return the natural logarithm of the likelihood
        function used by the sampler as an
        [additional dimensions] x niterations array.
        """
        return self.samples_chain['lnpost']

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
        # if samples are given then use those as initial positions
        if samples is not None:
            # transform to sampling parameter space
            samples = self.likelihood_evaluator.apply_sampling_transforms(
                samples)
        # draw random samples if samples are not provided
        else:
            samples = self.likelihood_evaluator.prior_rvs(size=1,prior=prior)

        self._p0 = samples
        return samples

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

            # Need to convert the sampling parameters to a list
            samples_list = [samples[param] for param in self.sampling_args]

            # The Jacobian is unused, is there a way to not compute it?
            logplr, (prior, loglr, logjacobian) = \
                                        self.likelihood_evaluator(samples_list)

            logplr = logplr if isinstance(logplr, numpy.float64) else logplr[0]
            prior = prior if isinstance(prior, numpy.float64) else prior[0]
            loglr = loglr if isinstance(loglr, numpy.float64) else loglr[0]

            start_sample=numpy.insert(samples_list,0,[logplr,loglr])
        else:
            start_sample=self.last_sample

        dtype=numpy.dtype([(name, None) for name in
                                    ['lnpost','lnlike']+self.sampling_args])
        self.samples_chain = numpy.empty(niterations,dtype=dtype)
        self.samples_chain[0] = start_sample

        for i in range(niterations-1):

            logplr_old,loglr_old = self.samples_chain[['lnpost','lnlike']][i]
            samples = self.samples_chain[self.sampling_args][i]

            # Dummy proposal
            samples_prop = [sample + numpy.random.normal(loc=0.0, scale=0.1)
                            for sample in samples]

            # The Jacobian is unused, is there a way to not compute it?
            logplr_prop, (prior_prop, loglr_prop, logjacobian) = \
                                        self.likelihood_evaluator(samples_prop)

            # There has to be a better way than the multiple if statements below
            if isinstance(prior_prop,numpy.float64):
                prior_prop = prior_prop
            elif prior_prop: # Because sometimes the likelihood returns None
                prior_prop = prior_prop[0]
            else:
                prior_prop = -numpy.inf

            if isinstance(logplr_prop,numpy.float64):
                logplr_prop = logplr_prop
            elif logplr_prop: # Because sometimes the likelihood returns None
                logplr_prop = logplr_prop[0]
            else:
                logplr_prop = -numpy.inf

            if isinstance(loglr_prop,numpy.float64):
                loglr_prop = loglr_prop
            elif loglr_prop: # Because sometimes the likelihood returns None
                loglr_prop = loglr_prop[0]
            else:
                loglr_prop = -numpy.inf

            acceptance_ratio=numpy.exp(logplr_prop - logplr_old)
            u=numpy.random.uniform()
            if acceptance_ratio >= u:
                self.samples_chain[i+1]=numpy.insert(samples_prop,0,
                                                    [logplr_prop,loglr_prop])
                logging.info("Step %i, acceptance ratio %f >= %f, accepted",
                                        i+1, acceptance_ratio, u)
            else:
                self.samples_chain[i+1]=self.samples_chain[i]
                logging.info("Step %i, acceptance ratio %f < %f, rejected",
                                        i+1, acceptance_ratio, u)


        return samples_prop, logplr_prop, loglr_prop

    @classmethod
    def calculate_logevidence(cls, fp):
        """This function should calculate the log evidence and its error using
        the results in the given file. If the sampler does not support evidence
        calculation, then this will raise a NotImplementedError.
        """
        raise NotImplementedError("this sampler does not support evidence "
                                  "calculation")

    def write_results(self, fp, start_iteration=0, end_iteration=None,
                      max_iterations=None, **metadata):
        """Writes metadata and samples to the given file.
        See the various write function for details.

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
        \**metadata :
            All other keyword arguments are passed to ``write_metadata``.
        """
        self.write_metadata(fp, **metadata)
        self.write_chain(fp, start_iteration=start_iteration,
                         end_iteration=end_iteration,
                         max_iterations=max_iterations)

    @staticmethod
    def write_samples_group(fp, samples_group, parameters, samples,
                            start_iteration=0, end_iteration=None,
                            index_offset=0, max_iterations=None):
        """Writes samples to the given file.

        Results are written to:

            `fp[samples_group/{vararg}]`,

        where `{vararg}` is the name of a variable arg.

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
        niterations = len(samples)
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

        group = samples_group + '/{name}'

        # loop over number of dimensions
        for param in parameters:
            dataset_name = group.format(name=param)
            try:
                if fb > fp[dataset_name].size:
                    # resize the dataset
                    fp[dataset_name].resize(fb, axis=0)
                fp[dataset_name][fa:fb] = samples[param][ma:mb]
            except KeyError:
                # dataset doesn't exist yet
                fp.create_dataset(dataset_name, (fb,),
                                  maxshape=(max_iterations,),
                                  dtype=samples[param].dtype)
                fp[dataset_name][fa:fb] = samples[param][ma:mb]

    def write_chain(self, fp, start_iteration=0, end_iteration=None,
                    max_iterations=None):
        """Writes the samples from the current chain to the given file.

        Results are written to:

            `fp[fp.samples_group/{field}]`,

        where `{field}` is the name of each field returned by
        `likelihood_stats`.

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

    def write_metadata(self, fp, **kwargs):
        """Writes metadata about this sampler to the given file. Metadata is
        written to the file's `attrs`.

        Parameters
        ----------
        fp : InferenceFile
            A file handler to an open inference file.
        \**kwargs :
            All keyword args are written to the file's ``attrs``.
        """
        super(MCMCSampler, self).write_metadata(fp, **kwargs)
        # line 375 of pycbc_inference requires nwalkers:
        #             burn_in_eval.update(sampler, fp)
        fp.attrs["nwalkers"] = 1

    def write_state(self, fp):
        """ Saves the state of the sampler in a file.
        """
        raise NotImplementedError("Writing state to file not implemented.")

    @classmethod
    def read_samples(cls, fp, parameters,
                     thin_start=None, thin_interval=None, thin_end=None,
                     iteration=None,
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
                                iteration=iteration)

    @staticmethod
    def _read_fields(fp, fields_group, fields, array_class,
                     thin_start=None, thin_interval=None, thin_end=None,
                     iteration=None):
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
        group = fields_group + '/{name}'
        for name in fields:
            these_arrays = fp[group.format(name=name)][get_index]
            arrays[name] = numpy.vstack(these_arrays)
        return array_class.from_kwargs(**arrays)


    @staticmethod
    def write_acls(fp, acls):
        # Dummy method as required by line 389 of pycbc_inference
        #                sampler.write_acls(fp, sampler.compute_acls(fp))
        return None
    @classmethod
    def compute_acls(cls, fp, start_index=None, end_index=None):
        # Dummy method as required by line 389 of pycbc_inference
        #                sampler.write_acls(fp, sampler.compute_acls(fp))
        return None
