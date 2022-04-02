# Copyright (C) 2018 Collin Capano
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# self.option) any later version.
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
"""Provides I/O support for multi-tempered sampler.
"""

import argparse
import numpy
from .base_mcmc import (CommonMCMCMetadataIO, thin_samples_for_writing,
                        _ensemble_get_index, _ensemble_get_walker_index,
                        _get_index)

class ParseTempsArg(argparse.Action):
    """Argparse action that will parse temps argument.

    If the provided argument is 'all', sets 'all' in the namespace dest. If a
    a sequence of numbers are provided, converts those numbers to ints before
    saving to the namespace.
    """
    def __init__(self, type=str, **kwargs): # pylint: disable=redefined-builtin
        # check that type is string
        if type != str:
            raise ValueError("the type for this action must be a string")
        super(ParseTempsArg, self).__init__(type=type, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        singlearg = isinstance(values, str)
        if singlearg:
            values = [values]
        if values[0] == 'all':
            # check that only a single value was provided
            if len(values) > 1:
                raise ValueError("if provide 'all', should not specify any "
                                 "other temps")
            temps = 'all'
        else:
            temps = []
            for val in values:
                try:
                    val = int(val)
                except ValueError:
                    pass
                temps.append(val)
            if singlearg:
                temps = temps[0]
        setattr(namespace, self.dest, temps)


class CommonMultiTemperedMetadataIO(CommonMCMCMetadataIO):
    """Adds support for reading/writing multi-tempered metadata to
    :py:class:`~pycbc.inference.io.base_mcmc.CommonMCMCMetadatIO`.
    """
    @property
    def ntemps(self):
        """Returns the number of temperatures used by the sampler."""
        return self[self.sampler_group].attrs['ntemps']

    def write_sampler_metadata(self, sampler):
        """Adds writing ntemps to file.
        """
        super(CommonMultiTemperedMetadataIO, self).write_sampler_metadata(
            sampler)
        self[self.sampler_group].attrs["ntemps"] = sampler.ntemps

    @staticmethod
    def extra_args_parser(parser=None, skip_args=None, **kwargs):
        """Adds --temps to MCMCIO parser.
        """
        if skip_args is None:
            skip_args = []
        parser, actions = CommonMCMCMetadataIO.extra_args_parser(
            parser=parser, skip_args=skip_args, **kwargs)
        if 'temps' not in skip_args:
            act = parser.add_argument(
                "--temps", nargs="+", default=0, action=ParseTempsArg,
                help="Get the given temperatures. May provide either a "
                     "sequence of integers specifying the temperatures to "
                     "plot, or 'all' for all temperatures. Default is to only "
                     "plot the coldest (= 0) temperature chain.")
            actions.append(act)
        return parser, actions


def write_samples(fp, samples, parameters=None, last_iteration=None,
                  samples_group=None, thin_by=None):
    """Writes samples to the given file.

    This works both for standard MCMC and ensemble MCMC samplers with
    parallel tempering.

    Results are written to ``samples_group/{vararg}``, where ``{vararg}``
    is the name of a model params. The samples are written as an
    ``ntemps x nwalkers x niterations`` array.

    Parameters
    -----------
    fp : BaseInferenceFile
        Open file handler to write files to. Must be an instance of
        BaseInferenceFile with CommonMultiTemperedMetadataIO methods added.
    samples : dict
        The samples to write. Each array in the dictionary should have
        shape ntemps x nwalkers x niterations.
    parameters : list, optional
        Only write the specified parameters to the file. If None, will
        write all of the keys in the ``samples`` dict.
    last_iteration : int, optional
        The iteration of the last sample. If the file's ``thinned_by``
        attribute is > 1, this is needed to determine where to start
        thinning the samples to match what has already been stored on disk.
    samples_group : str, optional
        Which group to write the samples to. Default (None) will result
        in writing to "samples".
    thin_by : int, optional
        Override the ``thinned_by`` attribute in the file with the given
        value. **Only set this if you are using this function to write
        something other than inference samples!**
    """
    ntemps, nwalkers, niterations = tuple(samples.values())[0].shape
    assert all(p.shape == (ntemps, nwalkers, niterations)
               for p in samples.values()), (
           "all samples must have the same shape")
    if samples_group is None:
        samples_group = fp.samples_group
    if parameters is None:
        parameters = list(samples.keys())
    # thin the samples
    samples = thin_samples_for_writing(fp, samples, parameters,
                                       last_iteration, samples_group,
                                       thin_by=thin_by)
    # loop over number of dimensions
    group = samples_group + '/{name}'
    for param in parameters:
        dataset_name = group.format(name=param)
        data = samples[param]
        # check that there's something to write after thinning
        if data.shape[2] == 0:
            # nothing to write, move along
            continue
        try:
            fp_niterations = fp[dataset_name].shape[-1]
            istart = fp_niterations
            istop = istart + data.shape[2]
            if istop > fp_niterations:
                # resize the dataset
                fp[dataset_name].resize(istop, axis=2)
        except KeyError:
            # dataset doesn't exist yet
            istart = 0
            istop = istart + data.shape[2]
            fp.create_dataset(dataset_name, (ntemps, nwalkers, istop),
                              maxshape=(ntemps, nwalkers, None),
                              dtype=data.dtype,
                              fletcher32=True)
        fp[dataset_name][:, :, istart:istop] = data


def read_raw_samples(fp, fields,
                     thin_start=None, thin_interval=None, thin_end=None,
                     iteration=None, temps='all', chains=None,
                     flatten=True, group=None):
    """Base function for reading samples from a collection of independent
    MCMC chains file with parallel tempering.

    This may collect differing numbering of samples from each chains,
    depending on the thinning settings for each chain. If not flattened the
    returned array will have dimensions requested temps x requested chains x
    max samples, where max samples is the largest number of samples retrieved
    from a single chain. Chains that retrieve fewer samples will be padded with
    ``numpy.nan``. If flattened, the NaNs are removed prior to returning.

    Parameters
    -----------
    fp : BaseInferenceFile
        Open file handler to read samples from. Must be an instance of
        BaseInferenceFile with CommonMultiTemperedMetadataIO methods added.
    fields : list
        The list of field names to retrieve.
    thin_start : array or int, optional
        Start reading from the given sample. May either provide an array
        indicating the start index for each chain, or an integer. If the
        former, the array must have the same length as the number of chains
        that will be retrieved. If the latter, the given value will be used
        for all chains. Default (None) is to use the file's ``thin_start``
        attribute.
    thin_interval : array or int, optional
        Only read every ``thin_interval``-th sample. May either provide an
        array indicating the interval to use for each chain, or an integer. If
        the former, the array must have the same length as the number of chains
        that will be retrieved. If the latter, the given value will be used for
        all chains. Default (None) is to use the file's ``thin_interval``
        attribute.
    thin_end : array or int, optional
        Stop reading at the given sample index. May either provide an
        array indicating the end index to use for each chain, or an integer. If
        the former, the array must have the same length as the number of chains
        that will be retrieved. If the latter, the given value will be used for
        all chains. Default (None) is to use the the file's ``thin_end``
        attribute.
    iteration : int, optional
        Only read the given iteration from all chains. If provided, it
        overrides the ``thin_(start|interval|end)`` options.
    temps : 'all' or (list of) int, optional
        The temperature index (or list of indices) to retrieve. To retrieve
        all temperates pass 'all', or a list of all of the temperatures.
        Default is 'all'.
    chains : (list of) int, optional
        Only read from the given chains. Default is to read all.
    flatten : bool, optional
        Remove NaNs and flatten the samples to 1D arrays before returning.
        Otherwise, the returned arrays will have shape (requested temps x
        requested chains x max requested iteration(s)), with chains that return
        fewer samples padded with NaNs. Default is True.
    group : str, optional
        The name of the group to read sample datasets from. Default is
        the file's ``samples_group``.

    Returns
    -------
    dict
        A dictionary of field name -> numpy array pairs.
    """
    if isinstance(fields, str):
        fields = [fields]
    if group is None:
        group = fp.samples_group
    group = group + '/{name}'
    # chains to load
    if chains is None:
        chains = numpy.arange(fp.nchains)
    elif not isinstance(chains, (list, numpy.ndarray)):
        chains = numpy.array([chains]).astype(int)
    get_index = _get_index(fp, chains, thin_start, thin_interval, thin_end,
                           iteration)
    # load the samples
    arrays = {}
    for name in fields:
        dset = group.format(name=name)
        # get the temperatures to load
        tidx, selecttemps, ntemps = _get_temps_index(temps, fp, dset)
        alist = []
        maxiters = 0
        for ii, cidx in enumerate(chains):
            idx = get_index[ii]
            # load the data
            thisarr = fp[dset][tidx, cidx, idx]
            if thisarr.size == 0:
                # no samples were loaded; skip this chain
                alist.append(None)
                continue
            if isinstance(idx, (int, numpy.int_)):
                # make sure the last dimension corresponds to iteration
                thisarr = thisarr.reshape(list(thisarr.shape)+[1])
            # pull out the temperatures we need
            if selecttemps:
                thisarr = thisarr[temps, ...]
            # make sure its 2D
            thisarr = thisarr.reshape(ntemps, thisarr.shape[-1])
            alist.append(thisarr)
            maxiters = max(maxiters, thisarr.shape[-1])
        # stack into a single array
        arr = numpy.full((ntemps, len(chains), maxiters), numpy.nan,
                         dtype=fp[dset].dtype)
        for ii, thisarr in enumerate(alist):
            if thisarr is not None:
                arr[:, ii, :thisarr.shape[-1]] = thisarr
        if flatten:
            # flatten and remove nans
            arr = arr.flatten()
            arr = arr[~numpy.isnan(arr)]
        arrays[name] = arr
    return arrays


def ensemble_read_raw_samples(fp, fields, thin_start=None,
                              thin_interval=None, thin_end=None,
                              iteration=None, temps='all', walkers=None,
                              flatten=True, group=None):
    """Base function for reading samples from ensemble MCMC file with
    parallel tempering.

    Parameters
    -----------
    fp : BaseInferenceFile
        Open file handler to write files to. Must be an instance of
        BaseInferenceFile with CommonMultiTemperedMetadataIO methods added.
    fields : list
        The list of field names to retrieve.
    thin_start : int, optional
        Start reading from the given iteration. Default is to start from
        the first iteration.
    thin_interval : int, optional
        Only read every ``thin_interval`` -th sample. Default is 1.
    thin_end : int, optional
        Stop reading at the given iteration. Default is to end at the last
        iteration.
    iteration : int, optional
        Only read the given iteration. If this provided, it overrides
        the ``thin_(start|interval|end)`` options.
    temps : 'all' or (list of) int, optional
        The temperature index (or list of indices) to retrieve. To retrieve
        all temperates pass 'all', or a list of all of the temperatures.
        Default is 'all'.
    walkers : (list of) int, optional
        Only read from the given walkers. Default (``None``) is to read all.
    flatten : bool, optional
        Flatten the samples to 1D arrays before returning. Otherwise, the
        returned arrays will have shape (requested temps x
        requested walkers x requested iteration(s)). Default is True.
    group : str, optional
        The name of the group to read sample datasets from. Default is
        the file's ``samples_group``.

    Returns
    -------
    dict
        A dictionary of field name -> numpy array pairs.
    """
    if isinstance(fields, str):
        fields = [fields]
    # walkers to load
    widx, nwalkers = _ensemble_get_walker_index(fp, walkers)
    # get the slice to use
    get_index = _ensemble_get_index(fp, thin_start, thin_interval, thin_end,
                                    iteration)
    # load
    if group is None:
        group = fp.samples_group
    group = group + '/{name}'
    arrays = {}
    for name in fields:
        dset = group.format(name=name)
        tidx, selecttemps, ntemps = _get_temps_index(temps, fp, dset)
        arr = fp[dset][tidx, widx, get_index]
        niterations = arr.shape[-1] if iteration is None else 1
        if selecttemps:
            # pull out the temperatures we need
            arr = arr[temps, ...]
        if flatten:
            arr = arr.flatten()
        else:
            # ensure that the returned array is 3D
            arr = arr.reshape((ntemps, nwalkers, niterations))
        arrays[name] = arr
    return arrays


def _get_temps_index(temps, fp, dataset):
    """Convenience function to determine which temperatures to load.

    Parameters
    -----------
    temps : 'all' or (list of) int
        The temperature index (or list of indices) to retrieve. To retrieve
        all temperates pass 'all', or a list of all of the temperatures.
    fp : BaseInferenceFile
        Open file handler to read samples from. Must be an instance of
        BaseInferenceFile with CommonMultiTemperedMetadataIO methods added.
    dataset : str
        The name of the dataset that samples will be loaded from.

    Returns
    -------
    tidx : slice or list of int
        The temperature indices to load from the file.
    selecttemps : bool
        Whether specific temperatures need to be pulled out of the samples
        array after it is loaded from the file.
    ntemps : int
        The number of temperatures that will be loaded.
    """
    if temps == 'all':
        # all temperatures were requested; just need to know how many
        ntemps = fp[dataset].shape[0]
        tidx = slice(None, None)
        selecttemps = False
    elif isinstance(temps, (int, numpy.int_)):
        # only a single temperature is requested
        ntemps = 1
        tidx = temps
        selecttemps = False
    else:
        # a select set of temperatures are requested
        tidx = slice(None, None)
        ntemps = len(temps)
        selecttemps = True
    return tidx, selecttemps, ntemps
