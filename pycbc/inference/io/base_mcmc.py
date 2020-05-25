# Copyright (C) 2016 Christopher M. Biwer, Collin Capano
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
"""Provides I/O that is specific to MCMC samplers.
"""

from __future__ import (absolute_import, division)

from six import string_types

import numpy
import argparse

class MCMCMetadataIO(object):
    """Provides functions for reading/writing MCMC metadata to file.
    """
    def write_resume_point(self):
        """Keeps a list of the number of iterations that were in a file when a
        run was resumed from a checkpoint."""
        try:
            resume_pts = self.attrs["resume_points"].tolist()
        except KeyError:
            resume_pts = []
        try:
            niterations = self.niterations
        except KeyError:
            niterations = 0
        resume_pts.append(niterations)
        self.attrs["resume_points"] = resume_pts

    def write_niterations(self, niterations):
        """Writes the given number of iterations to the sampler group."""
        self[self.sampler_group].attrs['niterations'] = niterations

    @property
    def niterations(self):
        """Returns the number of iterations the sampler was run for."""
        return self[self.sampler_group].attrs['niterations']

    @property
    def nwalkers(self):
        """Returns the number of walkers used by the sampler."""
        return self[self.sampler_group].attrs['nwalkers']

    def _thin_data(self, group, params, thin_interval):
        """Thins data on disk by the given interval.

        This makes no effort to record the thinning interval that is applied.

        Parameters
        ----------
        group : str
            The group where the datasets to thin live.
        params : list
            The list of dataset names to thin.
        thin_interval : int
            The interval to thin the samples on disk by.
        """
        samples = self.read_raw_samples(params, thin_start=0,
                                        thin_interval=thin_interval,
                                        thin_end=None, flatten=False,
                                        group=group)
        # now resize and write the data back to disk
        fpgroup = self[group]
        for param in params:
            data = samples[param]
            # resize the arrays on disk
            fpgroup[param].resize(data.shape)
            # and write
            fpgroup[param][:] = data

    def thin(self, thin_interval):
        """Thins the samples on disk to the given thinning interval.

        The interval must be a multiple of the file's current ``thinned_by``.

        Parameters
        ----------
        thin_interval : int
            The interval the samples on disk should be thinned by.
        """
        # get the new interval to thin by
        new_interval = thin_interval / self.thinned_by
        if new_interval % 1:
            raise ValueError("thin interval ({}) must be a multiple of the "
                             "current thinned_by ({})"
                             .format(thin_interval, self.thinned_by))
        new_interval = int(new_interval)
        # now thin the data on disk
        params = list(self[self.samples_group].keys())
        self._thin_data(self.samples_group, params, new_interval)
        # store the interval that samples were thinned by
        self.thinned_by = thin_interval
        # If a default thin interval and thin start exist, reduce them by the
        # thinned interval. If the thin interval is not an integer multiple
        # of the original, we'll round up, to avoid getting samples from
        # before the burn in / at an interval less than the ACL.
        self.thin_start = int(numpy.ceil(self.thin_start/new_interval))
        self.thin_interval = int(numpy.ceil(self.thin_interval/new_interval))

    @property
    def thinned_by(self):
        """Returns interval samples have been thinned by on disk.

        This looks for ``thinned_by`` in the samples group attrs. If none is
        found, will just return 1.
        """
        try:
            thinned_by = self.attrs['thinned_by']
        except KeyError:
            thinned_by = 1
        return thinned_by

    @thinned_by.setter
    def thinned_by(self, thinned_by):
        """Sets the thinned_by attribute.

        This is the interval that samples have been thinned by on disk. The
        given value is written to
        ``self[self.samples_group].attrs['thinned_by']``.
        """
        self.attrs['thinned_by'] = int(thinned_by)

    def last_iteration(self, parameter=None, group=None):
        """Returns the iteration of the last sample of the given parameter.

        Parameters
        ----------
        parameter : str, optional
            The name of the parameter to get the last iteration for. If
            None provided, will just use the first parameter in ``group``.
        group : str, optional
            The name of the group to get the last iteration from. Default is
            the ``samples_group``.
        """
        if group is None:
            group = self.samples_group
        if parameter is None:
            try:
                parameter = list(self[group].keys())[0]
            except (IndexError, KeyError):
                # nothing has been written yet, just return 0
                return 0
        try:
            lastiter = self[group][parameter].shape[-1]
        except KeyError:
            # no samples have been written, just return 0
            lastiter = 0
        # account for thinning
        return lastiter * self.thinned_by

    def iterations(self, parameter):
        """Returns the iteration each sample occurred at."""
        return numpy.arange(0, self.last_iteration(parameter), self.thinned_by)

    def write_sampler_metadata(self, sampler):
        """Writes the sampler's metadata."""
        self.attrs['sampler'] = sampler.name
        self[self.sampler_group].attrs['nwalkers'] = sampler.nwalkers
        # write the model's metadata
        sampler.model.write_metadata(self)

    def write_acls(self, acls):
        """Writes the given autocorrelation lengths.

        The ACL of each parameter is saved to
        ``[sampler_group]/acls/{param}']``.  The maximum over all the
        parameters is saved to the file's 'acl' attribute.

        Parameters
        ----------
        acls : dict
            A dictionary of ACLs keyed by the parameter.

        Returns
        -------
        ACL
            The maximum of the acls that was written to the file.
        """
        group = self.sampler_group + '/acls/{}'
        # write the individual acls
        for param in acls:
            try:
                # we need to use the write_direct function because it's
                # apparently the only way to update scalars in h5py
                self[group.format(param)].write_direct(
                    numpy.array(acls[param]))
            except KeyError:
                # dataset doesn't exist yet
                self[group.format(param)] = acls[param]
        # write the maximum over all params
        acl = numpy.array(list(acls.values())).max()
        self[self.sampler_group].attrs['acl'] = acl
        # set the default thin interval to be the acl (if it is finite)
        if numpy.isfinite(acl):
            self.thin_interval = int(numpy.ceil(acl))

    def read_acls(self):
        """Reads the acls of all the parameters.

        Returns
        -------
        dict
            A dictionary of the ACLs, keyed by the parameter name.
        """
        group = self[self.sampler_group]['acls']
        return {param: group[param].value for param in group.keys()}

    def write_burn_in(self, burn_in):
        """Write the given burn-in data to the given filename."""
        group = self[self.sampler_group]
        group.attrs['burn_in_test'] = burn_in.burn_in_test
        group.attrs['is_burned_in'] = burn_in.is_burned_in
        group.attrs['burn_in_iteration'] = burn_in.burn_in_iteration
        # set the defaut thin_start to be the burn_in_index
        self.thin_start = burn_in.burn_in_index
        # write individual test data
        for tst in burn_in.burn_in_data:
            key = 'burn_in_tests/{}'.format(tst)
            try:
                attrs = group[key].attrs
            except KeyError:
                group.create_group(key)
                attrs = group[key].attrs
            self.write_kwargs_to_attrs(attrs, **burn_in.burn_in_data[tst])

    @staticmethod
    def extra_args_parser(parser=None, skip_args=None, **kwargs):
        """Create a parser to parse sampler-specific arguments for loading
        samples.

        Parameters
        ----------
        parser : argparse.ArgumentParser, optional
            Instead of creating a parser, add arguments to the given one. If
            none provided, will create one.
        skip_args : list, optional
            Don't parse the given options. Options should be given as the
            option string, minus the '--'. For example,
            ``skip_args=['iteration']`` would cause the ``--iteration``
            argument not to be included.
        \**kwargs :
            All other keyword arguments are passed to the parser that is
            created.

        Returns
        -------
        parser : argparse.ArgumentParser
            An argument parser with th extra arguments added.
        actions : list of argparse.Action
            A list of the actions that were added.
        """
        if parser is None:
            parser = argparse.ArgumentParser(**kwargs)
        elif kwargs:
            raise ValueError("No other keyword arguments should be provded if "
                             "a parser is provided.")
        if skip_args is None:
            skip_args = []
        actions = []
        if 'thin-start' not in skip_args:
            act = parser.add_argument(
                "--thin-start", type=int, default=None,
                help="Sample number to start collecting samples to plot. If "
                     "none provided, will use the input file's `thin_start` "
                     "attribute.")
            actions.append(act)
        if 'thin-interval' not in skip_args:
            act = parser.add_argument(
                "--thin-interval", type=int, default=None,
                help="Interval to use for thinning samples. If none provided, "
                     "will use the input file's `thin_interval` attribute.")
            actions.append(act)
        if 'thin-end' not in skip_args:
            act = parser.add_argument(
                "--thin-end", type=int, default=None,
                help="Sample number to stop collecting samples to plot. If "
                     "none provided, will use the input file's `thin_end` "
                     "attribute.")
            actions.append(act)
        if 'iteration' not in skip_args:
            act = parser.add_argument(
                "--iteration", type=int, default=None,
                help="Only retrieve the given iteration. To load "
                     "the last n-th sampe use -n, e.g., -1 will "
                     "load the last iteration. This overrides "
                     "the thin-start/interval/end options.")
            actions.append(act)
        if 'walkers' not in skip_args:
            act = parser.add_argument(
                "--walkers", type=int, nargs="+", default=None,
                help="Only retrieve samples from the listed "
                     "walkers. Default is to retrieve from all "
                     "walkers.")
            actions.append(act)
        return parser, actions


class SingleTempMCMCIO(object):
    """Provides functions for reading/writing samples from an MCMC sampler.

    These functions will work for samplers that have 1 or more walkers, with
    only a single temperature.
    """

    def write_samples(self, samples, parameters=None, last_iteration=None,
                      samples_group=None, thin_by=None):
        """Writes samples to the given file.

        Results are written to ``samples_group/{vararg}``, where ``{vararg}``
        is the name of a model params. The samples are written as an
        ``nwalkers x niterations`` array. If samples already exist, the new
        samples are appended to the current.

        If the current samples on disk have been thinned (determined by the
        ``thinned_by`` attribute in the samples group), then the samples will
        be thinned by the same amount before being written. The thinning is
        started at the sample in ``samples`` that occured at the iteration
        equal to the last iteration on disk plus the ``thinned_by`` interval.
        If this iteration is larger than the iteration of the last given
        sample, then none of the samples will be written.

        Parameters
        -----------
        samples : dict
            The samples to write. Each array in the dictionary should have
            shape nwalkers x niterations.
        parameters : list, optional
            Only write the specified parameters to the file. If None, will
            write all of the keys in the ``samples`` dict.
        last_iteration : int, optional
            The iteration of the last sample. If the file's ``thinned_by``
            attribute is > 1, this is needed to determine where to start
            thinning the samples such that the interval between the last sample
            currently on disk and the first new sample is the same as all of
            the other samples.
        samples_group : str, optional
            Which group to write the samples to. Default (None) will result
            in writing to "samples".
        thin_by : int, optional
            Override the ``thinned_by`` attribute in the file with the given
            value. **Only set this if you are using this function to write
            something other than inference samples!**
        """
        nwalkers, nsamples = list(samples.values())[0].shape
        assert all(p.shape == (nwalkers, nsamples)
                   for p in samples.values()), (
               "all samples must have the same shape")
        if samples_group is None:
            samples_group = self.samples_group
        if parameters is None:
            parameters = samples.keys()
        # thin the samples
        samples = thin_samples_for_writing(self, samples, parameters,
                                           last_iteration, samples_group,
                                           thin_by=thin_by)
        # loop over number of dimensions
        group = samples_group + '/{name}'
        for param in parameters:
            dataset_name = group.format(name=param)
            data = samples[param]
            # check that there's something to write after thinning
            if data.shape[1] == 0:
                # nothing to write, move along
                continue
            try:
                fp_nsamples = self[dataset_name].shape[-1]
                istart = fp_nsamples
                istop = istart + data.shape[1]
                if istop > fp_nsamples:
                    # resize the dataset
                    self[dataset_name].resize(istop, axis=1)
            except KeyError:
                # dataset doesn't exist yet
                istart = 0
                istop = istart + data.shape[1]
                self.create_dataset(dataset_name, (nwalkers, istop),
                                    maxshape=(nwalkers, None),
                                    dtype=data.dtype,
                                    fletcher32=True)
            self[dataset_name][:, istart:istop] = data

    def read_raw_samples(self, fields,
                         thin_start=None, thin_interval=None, thin_end=None,
                         iteration=None, walkers=None, flatten=True,
                         group=None):
        """Base function for reading samples.

        Parameters
        -----------
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
        walkers : int, optional
            Only read from the given walkers. Default is to read all.
        flatten : bool, optional
            Flatten the samples to 1D arrays before returning. Otherwise, the
            returned arrays will have shape (requested walkers x
            requested iteration(s)). Default is True.
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
            widx = numpy.zeros(self.nwalkers, dtype=bool)
            widx[walkers] = True
            nwalkers = widx.sum()
        else:
            widx = slice(0, None)
            nwalkers = self.nwalkers
        # get the slice to use
        if iteration is not None:
            get_index = int(iteration)
            niterations = 1
        else:
            get_index = self.get_slice(thin_start=thin_start,
                                       thin_end=thin_end,
                                       thin_interval=thin_interval)
            # we'll just get the number of iterations from the returned shape
            niterations = None
        # load
        if group is None:
            group = self.samples_group
        group = group + '/{name}'
        arrays = {}
        for name in fields:
            arr = self[group.format(name=name)][widx, get_index]
            if niterations is None:
                niterations = arr.shape[-1]
            if flatten:
                arr = arr.flatten()
            else:
                # ensure that the returned array is 2D
                arr = arr.reshape((nwalkers, niterations))
            arrays[name] = arr
        return arrays


def thin_samples_for_writing(fp, samples, parameters, last_iteration,
                             group, thin_by=None):
    """Thins samples for writing to disk.

    The thinning interval to use is determined by the given file handler's
    ``thinned_by`` attribute. If that attribute is 1, just returns the samples.

    Parameters
    ----------
    fp : MCMCMetadataIO instance
        The file the sampels will be written to. Needed to determine the
        thin interval used on disk.
    samples : dict
        Dictionary mapping parameter names to arrays of (unthinned) samples.
        The arrays are thinned along their last dimension.
    parameters : list of str
        The parameters to thin in ``samples`` before writing. All listed
        parameters must be in ``samples``.
    last_iteration : int
        The iteration that the last sample in ``samples`` occurred at. This is
        needed to figure out where to start the thinning in ``samples``, such
        that the interval between the last sample on disk and the first new
        sample is the same as all of the other samples.
    group : str
        The name of the group that the samples will be written to. This is
        needed to determine what the last iteration saved on disk was.
    thin_by : int, optional
        Override the ``thinned_by`` attribute in the file for with the given
        value. **Only do this if you are thinning something other than
        inference samples!**

    Returns
    -------
    dict :
        Dictionary of the thinned samples to write.
    """
    if thin_by is None:
        thin_by = fp.thinned_by
    if thin_by > 1:
        if last_iteration is None:
            raise ValueError("File's thinned_by attribute is > 1 ({}), "
                             "but last_iteration not provided."
                             .format(thin_by))
        thinned_samples = {}
        for param in parameters:
            data = samples[param]
            nsamples = data.shape[-1]
            # To figure out where to start:
            # the last iteration in the file + the file's thinning interval
            # gives the iteration of the next sample that should be written;
            # last_iteration - nsamples gives the iteration of the first
            # sample in samples. Subtracting the latter from the former - 1
            # (-1 to convert from iteration to index) therefore gives the index
            # in the samples data to start using samples.
            thin_start = fp.last_iteration(param, group) + thin_by \
                - (last_iteration - nsamples) - 1
            thinned_samples[param] = data[..., thin_start::thin_by]
    else:
        thinned_samples = samples
    return thinned_samples
