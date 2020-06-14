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

from __future__ import absolute_import
import argparse
from six import string_types
import numpy
from .base_mcmc import (CommonMCMCMetadataIO, thin_samples_for_writing,
                        nsamples_in_chain)

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
        singlearg = isinstance(values, string_types)
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
        super(MultiTemperedMetadataIO, self).write_sampler_metadata(sampler)
        self[self.sampler_group].attrs["ntemps"] = sampler.ntemps

    @staticmethod
    def extra_args_parser(parser=None, skip_args=None, **kwargs):
        """Adds --temps to MCMCIO parser.
        """
        if skip_args is None:
            skip_args = []
        parser, actions = MCMCMetadataIO.extra_args_parser(
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
                                maxshape=(ntemps, nwalkers,
                                          None),
                                dtype=data.dtype,
                                fletcher32=True)
        fp[dataset_name][:, :, istart:istop] = data


def read_raw_samples(self, fields,
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
    if isinstance(fields, string_types):
        fields = [fields]
    if group is None:
        group = self.samples_group
    group = group + '/{name}'
    # chains to load
    if chains is None:
        chains = numpy.arange(self.nchains)
    elif not isinstance(chains, (list, numpy.ndarray)):
        chains = numpy.array([chains]).astype(int)
    # temperatures to load
    selecttemps = False
    if isinstance(temps, (int, numpy.int32, numpy.int64)):
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
            ntemps = self.ntemps
    # iterations to load
    if iteration is not None:
        maxiters = 1
        get_index = [int(iteration)]*len(chains)
    else:
        if thin_start is None:
            thin_start = self.thin_start
        if not isinstance(thin_start, (numpy.ndarray, list)):
            thin_start = numpy.repeat(thin_start, len(chains), dtype=int)
        if thin_interval is None:
            thin_interval = self.thin_interval
        if not isinstance(thin_interval, (numpy.ndarray, list)):
            thin_interval = numpy.repeat(thin_interval, len(chains),
                                         dtype=int)
        if thin_end is None:
            thin_end = self.thin_end
        if not isinstance(thin_end, (numpy.ndarray, list)):
            thin_end = numpy.repeat(thin_end, len(chains))
        # figure out the maximum number of samples we will get from all chains
        maxiters = nsamples_in_chain(thin_start, thin_interval, thin_end).max()
        # the slices to use for each chain
        get_index = [self.get_slice(thin_start=thin_start[ci],
                                    thin_interval=thin_interval[ci],
                                    thin_end=thin_end[ci])
                     for ci in chains]
    # load the samples
    for name in fields:
        arr = numpy.full((ntemps, nchains, maxiter), numpy.nan)
        for ci in cidx:
            thisarr = self[group.format(name=name)][tidx, ci, get_index]
            # pull out the temperatures we need
            if selecttemps:
                thisarr = thisarr[temps, ...]
            # make sure its 2D
            thisarr = thisarr.reshape(ntemps, thisarr.shape[-1])
            arr[:, ci, :thisarr.shape[-1]] = thisarr
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
        Only read from the given walkers. Default is to read all.
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
    if isinstance(fields, string_types):
        fields = [fields]
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
    if isinstance(temps, (int, numpy.int32, numpy.int64)):
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
        get_index = int(iteration)
        niterations = 1
    else:
        if thin_start is None:
            thin_start = fp.thin_start
        if thin_interval is None:
            thin_interval = fp.thin_interval
        if thin_end is None:
            thin_end = fp.thin_end
        get_index = fp.get_slice(thin_start=thin_start,
                                 thin_interval=thin_interval,
                                 thin_end=thin_end)
        # we'll just get the number of iterations from the returned shape
        niterations = None
    # load
    if group is None:
        group = fp.samples_group
    group = group + '/{name}'
    arrays = {}
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
    return arrays
