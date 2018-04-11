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

from six import string_types

import numpy
from pycbc.inference.sampler_base import BaseMCMCSampler, _check_fileformat
from pycbc.io import FieldArray
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
    pool : function with map, Optional
        A provider of a map function that allows a function call to be run
        over multiple sets of arguments and possibly maps them to
        cores/nodes/etc.
    """
    name = "emcee"

    def __init__(self, likelihood_evaluator, nwalkers, pool=None,
                 likelihood_call=None):
        try:
            import emcee
        except ImportError:
            raise ImportError("emcee is not installed.")

        if likelihood_call is None:
            likelihood_call = likelihood_evaluator

        ndim = len(likelihood_evaluator.variable_args)
        sampler = emcee.EnsembleSampler(nwalkers, ndim,
                                        likelihood_call,
                                        pool=pool)
        # emcee uses it's own internal random number generator; we'll set it
        # to have the same state as the numpy generator
        rstate = numpy.random.get_state()
        sampler.random_state = rstate
        # initialize
        super(EmceeEnsembleSampler, self).__init__(
              sampler, likelihood_evaluator)
        self._nwalkers = nwalkers

    @classmethod
    def from_cli(cls, opts, likelihood_evaluator, pool=None,
                 likelihood_call=None):
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
        return cls(likelihood_evaluator, opts.nwalkers,
                   pool=pool, likelihood_call=likelihood_call)

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

    def clear_chain(self):
        """Clears the chain and blobs from memory.
        """
        # store the iteration that the clear is occuring on
        self.lastclear = self.niterations
        # now clear the chain
        self._sampler.reset()
        self._sampler.clear_blobs()

    def set_p0(self, samples_file=None, prior=None):
        """Sets the initial position of the walkers.

        Parameters
        ----------
        samples_file : InferenceFile, optional
            If provided, use the last iteration in the given file for the
            starting positions.
        prior : JointDistribution, optional
            Use the given prior to set the initial positions rather than
            `likelihood_evaultor`'s prior.

        Returns
        -------
        p0 : array
            An nwalkers x ndim array of the initial positions that were set.
        """
        # we define set_p0 here to ensure that emcee's internal random number
        # generator is set to numpy's after the distributions' rvs functions
        # are called
        super(EmceeEnsembleSampler, self).set_p0(samples_file=samples_file,
            prior=prior)
        # update the random state
        self._sampler.random_state = numpy.random.get_state()

    def write_state(self, fp):
        """Saves the state of the sampler in a file.
        """
        fp.write_random_state(state=self._sampler.random_state)

    def set_state_from_file(self, fp):
        """Sets the state of the sampler back to the instance saved in a file.
        """
        rstate = fp.read_random_state()
        # set the numpy random state
        numpy.random.set_state(rstate)
        # set emcee's generator to the same state
        self._sampler.random_state = rstate

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

    def write_results(self, fp, start_iteration=None,
                      max_iterations=None, **metadata):
        """Writes metadata, samples, likelihood stats, and acceptance fraction
        to the given file. See the write function for each of those for
        details.

        Parameters
        -----------
        fp : InferenceFile
            A file handler to an open inference file.
        start_iteration : int, optional
            Write results to the file's datasets starting at the given
            iteration. Default is to append after the last iteration in the
            file.
        max_iterations : int, optional
            Set the maximum size that the arrays in the hdf file may be resized
            to. Only applies if the samples have not previously been written
            to file. The default (None) is to use the maximum size allowed by
            h5py.
        \**metadata :
            All other keyword arguments are passed to ``write_metadata``.
        """
        self.write_metadata(fp, **metadata)
        self.write_chain(fp, start_iteration=start_iteration,
                         max_iterations=max_iterations)
        self.write_likelihood_stats(fp, start_iteration=start_iteration,
                                    max_iterations=max_iterations)
        self.write_acceptance_fraction(fp)
        self.write_state(fp)

# This is needed for two reason
# 1) pools freeze state when created and so classes *cannot be updated*
# 2) methods cannot be pickled.
class _callprior(object):
    """Calls the likelihood function's prior function, and ensures that no
    metadata is returned."""
    def __init__(self, likelihood_call):
        self.callable = likelihood_call

    def __call__(self, args):
        prior = self.callable(args, callfunc='prior')
        return prior if isinstance(prior, numpy.float64) else prior[0]

class _callloglikelihood(object):
    """Calls the likelihood function's loglikelihood function.
    """
    def __init__(self, likelihood_call):
        self.callable = likelihood_call

    def __call__(self, args):
        return self.callable(args, callfunc='loglikelihood')


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
    pool : function with map, Optional
        A provider of a map function that allows a function call to be run
        over multiple sets of arguments and possibly maps them to
        cores/nodes/etc.
    """
    name = "emcee_pt"

    def __init__(self, likelihood_evaluator, ntemps, nwalkers, pool=None,
                 likelihood_call=None):

        try:
            import emcee
        except ImportError:
            raise ImportError("emcee is not installed.")

        if likelihood_call is None:
            likelihood_call = likelihood_evaluator

        # construct the sampler: PTSampler needs the likelihood and prior
        # functions separately
        ndim = len(likelihood_evaluator.variable_args)
        sampler = emcee.PTSampler(ntemps, nwalkers, ndim,
                                  _callloglikelihood(likelihood_call),
                                  _callprior(likelihood_call),
                                  pool=pool)
        # initialize
        super(EmceePTSampler, self).__init__(
              sampler, likelihood_evaluator)
        self._nwalkers = nwalkers
        self._ntemps = ntemps

    @classmethod
    def from_cli(cls, opts, likelihood_evaluator, pool=None,
                 likelihood_call=None):
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
        return cls(likelihood_evaluator, opts.ntemps, opts.nwalkers,
                   pool=pool, likelihood_call=likelihood_call)

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

    def clear_chain(self):
        """Clears the chain and blobs from memory.
        """
        # store the iteration that the clear is occuring on
        self.lastclear = self.niterations
        # now clear the chain
        self._sampler.reset()

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
        kwargs = {'loglr': loglr, 'prior': logp}
        # if different coordinates were used for sampling, get the jacobian
        if self.likelihood_evaluator.sampling_transforms is not None:
            samples = self.samples
            # convert to dict
            d = {param: samples[param] for param in samples.fieldnames}
            logj = self.likelihood_evaluator.logjacobian(**d)
            kwargs['logjacobian'] = logj
        return FieldArray.from_kwargs(**kwargs)

    @property
    def lnpost(self):
        """Get the natural logarithm of the likelihood + the prior as an
        ntemps x nwalkers x niterations array.
        """
        # emcee returns ntemps x nwalkers x niterations
        return self._sampler.lnprobability

    def set_p0(self, samples_file=None, prior=None):
        """Sets the initial position of the walkers.

        Parameters
        ----------
        samples_file : InferenceFile, optional
            If provided, use the last iteration in the given file for the
            starting positions.
        prior : JointDistribution, optional
            Use the given prior to set the initial positions rather than
            `likelihood_evaultor`'s prior.

        Returns
        -------
        p0 : array
            An ntemps x nwalkers x ndim array of the initial positions that
            were set.
        """
        # create a (nwalker, ndim) array for initial positions
        ntemps = self.ntemps
        nwalkers = self.nwalkers
        ndim = len(self.variable_args)
        p0 = numpy.ones((ntemps, nwalkers, ndim))
        # if samples are given then use those as initial positions
        if samples_file is not None:
            samples = self.read_samples(samples_file, self.variable_args,
                iteration=-1, temps='all', flatten=False)[..., 0]
            # transform to sampling parameter space
            samples = self.likelihood_evaluator.apply_sampling_transforms(
                samples)
        # draw random samples if samples are not provided
        else:
            samples = self.likelihood_evaluator.prior_rvs(
                size=nwalkers*ntemps, prior=prior).reshape((ntemps, nwalkers))
        # convert to array
        for i, param in enumerate(self.sampling_args):
            p0[..., i] = samples[param]
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

    # read/write functions

    # add ntemps and betas to metadata
    def write_metadata(self, fp, **kwargs):
        """Writes metadata about this sampler to the given file. Metadata is
        written to the file's `attrs`.

        Parameters
        ----------
        fp : InferenceFile
            A file handler to an open inference file.
        \**kwargs :
            All keyword arguments are saved as separate arguments in the
            file attrs. If any keyword argument is a dictionary, the keyword
            will point to the list of keys in the the file's ``attrs``. Each
            key is then stored as a separate attr with its corresponding value.
        """
        super(EmceePTSampler, self).write_metadata(fp, **kwargs)
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
            arrays.extend(fp[group.format(tk=tk)][wmask])
        return arrays

    @staticmethod
    def write_samples_group(fp, samples_group, parameters, samples,
                             start_iteration=None, max_iterations=None):
        """Writes samples to the given file.

        Results are written to:

            ``fp[samples_group/{vararg}]``,

        where ``{vararg}`` is the name of a variable arg. The samples are
        written as an ``ntemps x nwalkers x niterations`` array.

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
        start_iteration : int, optional
            Write results to the file's datasets starting at the given
            iteration. Default is to append after the last iteration in the
            file.
        max_iterations : int, optional
            Set the maximum size that the arrays in the hdf file may be resized
            to. Only applies if the samples have not previously been written
            to file. The default (None) is to use the maximum size allowed by
            h5py.
        """
        ntemps, nwalkers, niterations = samples.shape
        if max_iterations is not None and max_iterations < niterations:
            raise IndexError("The provided max size is less than the "
                             "number of iterations")
        group = samples_group + '/{name}'
        # loop over number of dimensions
        for param in parameters:
            dataset_name = group.format(name=param)
            istart = start_iteration
            try:
                fp_niterations = fp[dataset_name].shape[-1]
                if istart is None:
                    istart = fp_niterations
                istop = istart + niterations
                if istop > fp_niterations:
                    # resize the dataset
                    fp[dataset_name].resize(istop, axis=2)
            except KeyError:
                # dataset doesn't exist yet
                if istart is not None and istart != 0:
                    raise ValueError("non-zero start_iteration provided, but "
                                     "dataset doesn't exist yet")
                istart = 0
                istop = istart + niterations
                fp.create_dataset(dataset_name, (ntemps, nwalkers, istop),
                                  maxshape=(ntemps, nwalkers, max_iterations),
                                  dtype=float, fletcher32=True)
            fp[dataset_name][:,:,istart:istop] = samples[param]

    def write_results(self, fp, start_iteration=None, max_iterations=None,
                      **metadata):
        """Writes metadata, samples, likelihood stats, and acceptance fraction
        to the given file. See the write function for each of those for
        details.

        Parameters
        -----------
        fp : InferenceFile
            A file handler to an open inference file.
        start_iteration : int, optional
            Write results to the file's datasets starting at the given
            iteration. Default is to append after the last iteration in the
            file.
        max_iterations : int, optional
            Set the maximum size that the arrays in the hdf file may be resized
            to. Only applies if the samples have not previously been written
            to file. The default (None) is to use the maximum size allowed by
            h5py.
        \**metadata :
            All other keyword arguments are passed to ``write_metadata``.
        """
        self.write_metadata(fp, **metadata)
        self.write_chain(fp, start_iteration=start_iteration,
                         max_iterations=max_iterations)
        self.write_likelihood_stats(fp, start_iteration=start_iteration,
                                    max_iterations=max_iterations)
        self.write_acceptance_fraction(fp)
        self.write_state(fp)

    @staticmethod
    def _read_oldstyle_fields(fp, fields_group, fields, array_class,
                     thin_start=None, thin_interval=None, thin_end=None,
                     iteration=None, temps=None, walkers=None, flatten=True):
        """Base function for reading samples and likelihood stats. See
        `read_samples` and `read_likelihood_stats` for details.

        This function is to provide backward compatability with older files.
        This will be removed in a future update.

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
            get_index = [iteration]
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
        if walkers is not None:
            widx = numpy.zeros(fp.nwalkers, dtype=bool)
            widx[walkers] = True
            nwalkers = widx.sum()
        else:
            widx = slice(None, None)
            nwalkers = fp.nwalkers
        # temperatures to load
        selecttemps = False
        if temps is None:
            tidx = 0
            ntemps = 1
        elif isinstance(temps, int):
            tidx = temps
            ntemps = 1
        else:
            # temps is either 'all' or a list of temperatures;
            # in either case, we'll get all of the temperatures from the file;
            # if not 'all', then we'll pull out the ones we want
            tidx = slice(None, None)
            selecttemps = temps != 'all'
            if selecttemps:
                ntemps = len(temps)
            else:
                ntemps = fp.ntemps
        # get the slice to use
        if iteration is not None:
            get_index = iteration
            niterations = 1
        else:
            if thin_end is None:
                # use the number of current iterations
                thin_end = fp.niterations
            get_index = fp.get_slice(thin_start=thin_start, thin_end=thin_end,
                                     thin_interval=thin_interval)
            # we'll just get the number of iterations from the returned shape
            niterations = None
        # load
        arrays = {}
        group = fields_group + '/{name}'
        for name in fields:
            arr = fp[group.format(name=name)][tidx, widx, get_index]
            if niterations is None:
                niterations = arr.shape[-1]
            # pull out the temperatures we need
            if selecttemps:
                arr = arr[temps, ...]
            if flatten:
                arr = arr.flatten()
            else:
                # ensure that the returned array is 3D
                arr = arr.reshape((ntemps, nwalkers, niterations))
            arrays[name] = arr
        return array_class.from_kwargs(**arrays)

    @classmethod
    @_check_fileformat
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
        return cls._read_fields(
                fp, samples_group, loadfields, array_class,
                thin_start=thin_start, thin_interval=thin_interval,
                thin_end=thin_end, iteration=iteration, temps=temps,
                walkers=walkers, flatten=flatten)

    @classmethod
    def compute_acfs(cls, fp, start_index=None, end_index=None,
                     per_walker=False, walkers=None, parameters=None,
                     temps=None):
        """Computes the autocorrleation function of the variable args in the
        given file.

        By default, parameter values are averaged over all walkers at each
        iteration. The ACF is then calculated over the averaged chain for each
        temperature. An ACF per-walker will be returned instead if
        ``per_walker=True``.

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
        temps : optional, (list of) int or 'all'
            The temperature index (or list of indices) to retrieve. If None
            (the default), the ACF will only be computed for the coldest (= 0)
            temperature chain. To compute an ACF for all temperates pass 'all',
            or a list of all of the temperatures.

        Returns
        -------
        FieldArray
            A ``FieldArray`` of the ACF vs iteration for each parameter. If
            `per-walker` is True, the FieldArray will have shape
            ``ntemps x nwalkers x niterations``. Otherwise, the returned
            array will have shape ``ntemps x niterations``.
        """
        acfs = {}
        if parameters is None:
            parameters = fp.variable_args
        if isinstance(parameters, string_types):
            parameters = [parameters]
        if isinstance(temps, int):
            temps = [temps]
        elif temps == 'all':
            temps = numpy.arange(fp.ntemps)
        elif temps is None:
            temps = [0]
        for param in parameters:
            subacfs = []
            for tk in temps:
                if per_walker:
                    # just call myself with a single walker
                    if walkers is None:
                        walkers = numpy.arange(fp.nwalkers)
                    arrays = [cls.compute_acfs(fp, start_index=start_index,
                                               end_index=end_index,
                                               per_walker=False, walkers=ii,
                                               parameters=param,
                                               temps=tk)[param][0,:]
                              for ii in walkers]
                    # we'll stack all of the walker arrays to make a single
                    # nwalkers x niterations array; when these are stacked
                    # below, we'll get a ntemps x nwalkers x niterations array
                    subacfs.append(numpy.vstack(arrays))
                else:
                    samples = cls.read_samples(fp, param,
                                               thin_start=start_index,
                                               thin_interval=1,
                                               thin_end=end_index,
                                               walkers=walkers, temps=tk,
                                               flatten=False)[param]
                    # contract the walker dimension using the mean, and flatten
                    # the (length 1) temp dimension
                    samples = samples.mean(axis=1)[0,:]
                    thisacf = autocorrelation.calculate_acf(samples).numpy()
                    subacfs.append(thisacf)
            # stack the temperatures
            # FIXME: the following if/else can be condensed to a single line
            # using numpy.stack, once the version requirements are bumped to
            # numpy >= 1.10
            if per_walker:
                nw, ni = subacfs[0].shape
                acfs[param] = numpy.zeros((len(temps), nw, ni), dtype=float)
                for tk in range(len(temps)):
                    acfs[param][tk,...] = subacfs[tk]
            else:
                acfs[param] = numpy.vstack(subacfs)
        return FieldArray.from_kwargs(**acfs)

    @classmethod
    def compute_acls(cls, fp, start_index=None, end_index=None):
        """Computes the autocorrleation length for all variable args and
        temperatures in the given file.

        Parameter values are averaged over all walkers at each iteration and
        temperature.  The ACL is then calculated over the averaged chain. If
        the returned ACL is `inf`,  will default to the number of current
        iterations.

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
            A dictionary of ntemps-long arrays of the ACLs of each parameter.
        """
        acls = {}
        if end_index is None:
            end_index = fp.niterations
        tidx = numpy.arange(fp.ntemps)
        for param in fp.variable_args:
            these_acls = numpy.zeros(fp.ntemps, dtype=int)
            for tk in tidx:
                samples = cls.read_samples(fp, param, thin_start=start_index,
                                           thin_interval=1, thin_end=end_index,
                                           temps=tk, flatten=False)[param]
                # contract the walker dimension using the mean, and flatten
                # the (length 1) temp dimension
                samples = samples.mean(axis=1)[0,:]
                acl = autocorrelation.calculate_acl(samples)
                if numpy.isinf(acl):
                    acl = samples.size
                these_acls[tk] = acl
            acls[param] = these_acls
        return acls

    @staticmethod
    def _oldstyle_read_acls(fp):
        """Deprecated: reads acls from older style files.

        Parameters
        ----------
        fp : InferenceFile
            An open file handler to read the acls from.

        Returns
        -------
        FieldArray
            An ntemps-long ``FieldArray`` containing the acls for every
            temperature, with the variable arguments as fields.
        """
        group = fp.samples_group + '/{param}/temp{tk}'
        tidx = numpy.arange(fp.ntemps)
        arrays = {}
        for param in fp.variable_args:
            arrays[param] = numpy.array([
                fp[group.format(param=param, tk=tk)].attrs['acl']
                for tk in tidx])
        return FieldArray.from_kwargs(**arrays)

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

        stats_group = fp.stats_group
        parameters = fp[stats_group].keys()
        logstats = cls.read_samples(fp, parameters, samples_group=stats_group,
                                    thin_start=thin_start,  thin_end=thin_end,
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
