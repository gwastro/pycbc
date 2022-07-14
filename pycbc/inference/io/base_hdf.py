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
"""This modules defines functions for reading and writing samples that the
inference samplers generate.
"""


import sys
import logging
from io import StringIO

from abc import (ABCMeta, abstractmethod)

import numpy
import h5py

from pycbc.io import FieldArray
from pycbc.inject import InjectionSet
from pycbc.io import (dump_state, load_state)
from pycbc.workflow import WorkflowConfigParser
from pycbc.types import FrequencySeries


def format_attr(val):
    """Formats an attr so that it can be read in either python 2 or 3.

    In python 2, strings that are saved as an attribute in an hdf file default
    to unicode. Since unicode was removed in python 3, if you load that file
    in a python 3 environment, the strings will be read as bytes instead, which
    causes a number of issues. This attempts to fix that. If the value is
    a bytes string, then it will be decoded into a string. If the value is
    a numpy array of byte strings, it will convert the array to a list of
    strings.

    Parameters
    ----------
    val : obj
        The value to format. This will try to apply decoding to the value

    Returns
    -------
    obj
        If ``val`` was a byte string, the value as a ``str``. If the value
        was a numpy array of ``bytes_``, the value as a list of ``str``.
        Otherwise, just returns the value.
    """
    try:
        val = str(val.decode())
    except AttributeError:
        pass
    if isinstance(val, numpy.ndarray) and val.dtype.type == numpy.bytes_:
        val = val.astype(numpy.unicode_).tolist()
    return val


class BaseInferenceFile(h5py.File, metaclass=ABCMeta):
    """Base class for all inference hdf files.

    This is a subclass of the h5py.File object. It adds functions for
    handling reading and writing the samples from the samplers.

    Parameters
    -----------
    path : str
        The path to the HDF file.
    mode : {None, str}
        The mode to open the file, eg. "w" for write and "r" for read.
    """

    name = None
    samples_group = 'samples'
    sampler_group = 'sampler_info'
    data_group = 'data'
    injections_group = 'injections'
    config_group = 'config_file'

    def __init__(self, path, mode=None, **kwargs):
        super(BaseInferenceFile, self).__init__(path, mode, **kwargs)
        # check that file type matches self
        try:
            filetype = self.attrs['filetype']
        except KeyError:
            if mode == 'w':
                # first time creating the file, add this class's name
                filetype = self.name
                self.attrs['filetype'] = filetype
            else:
                filetype = None
        try:
            filetype = str(filetype.decode())
        except AttributeError:
            pass
        if filetype != self.name:
            raise ValueError("This file has filetype {}, whereas this class "
                             "is named {}. This indicates that the file was "
                             "not written by this class, and so cannot be "
                             "read by this class.".format(filetype, self.name))

    def __getattr__(self, attr):
        """Things stored in ``.attrs`` are promoted to instance attributes.

        Note that properties will be called before this, so if there are any
        properties that share the same name as something in ``.attrs``, that
        property will get returned.
        """
        return self.attrs[attr]

    def getattrs(self, group=None, create_missing=True):
        """Convenience function for getting the `attrs` from the file or group.

        Parameters
        ----------
        group : str, optional
            Get the attrs of the specified group. If None or ``/``, will
            retrieve the file's ``attrs``.
        create_missing: bool, optional
            If ``group`` is provided, but doesn't yet exist in the file, create
            the group. Otherwise, a KeyError will be raised. Default is True.

        Returns
        -------
        h5py.File.attrs
            An attrs instance of the file or requested group.
        """
        if group is None or group == "/":
            attrs = self.attrs
        else:
            try:
                attrs = self[group].attrs
            except KeyError as e:
                if create_missing:
                    self.create_group(group)
                    attrs = self[group].attrs
                else:
                    raise e
        return attrs

    @abstractmethod
    def write_samples(self, samples, **kwargs):
        """This should write all of the provided samples.

        This function should be used to write both samples and model stats.

        Parameters
        ----------
        samples : dict
            Samples should be provided as a dictionary of numpy arrays.
        \**kwargs :
            Any other keyword args the sampler needs to write data.
        """
        pass

    def parse_parameters(self, parameters, array_class=None):
        """Parses a parameters arg to figure out what fields need to be loaded.

        Parameters
        ----------
        parameters : (list of) strings
            The parameter(s) to retrieve. A parameter can be the name of any
            field in ``samples_group``, a virtual field or method of
            ``FieldArray`` (as long as the file contains the necessary fields
            to derive the virtual field or method), and/or a function of
            these.
        array_class : array class, optional
            The type of array to use to parse the parameters. The class must
            have a ``parse_parameters`` method. Default is to use a
            ``FieldArray``.

        Returns
        -------
        list :
            A list of strings giving the fields to load from the file.
        """
        # get the type of array class to use
        if array_class is None:
            array_class = FieldArray
        # get the names of fields needed for the given parameters
        possible_fields = self[self.samples_group].keys()
        return array_class.parse_parameters(parameters, possible_fields)

    def read_samples(self, parameters, array_class=None, **kwargs):
        """Reads samples for the given parameter(s).

        The ``parameters`` can be the name of any dataset in ``samples_group``,
        a virtual field or method of ``FieldArray`` (as long as the file
        contains the necessary fields to derive the virtual field or method),
        and/or any numpy function of these.

        The ``parameters`` are parsed to figure out what datasets are needed.
        Only those datasets will be loaded, and will be the base-level fields
        of the returned ``FieldArray``.

        The ``static_params`` are also added as attributes of the returned
        ``FieldArray``.

        Parameters
        -----------
        parameters : (list of) strings
            The parameter(s) to retrieve.
        array_class : FieldArray-like class, optional
            The type of array to return. The class must have ``from_kwargs``
            and ``parse_parameters`` methods. If None, will return a
            ``FieldArray``.
        \**kwargs :
            All other keyword arguments are passed to ``read_raw_samples``.

        Returns
        -------
        FieldArray :
            The samples as a ``FieldArray``.
        """
        # get the type of array class to use
        if array_class is None:
            array_class = FieldArray
        # get the names of fields needed for the given parameters
        possible_fields = self[self.samples_group].keys()
        loadfields = array_class.parse_parameters(parameters, possible_fields)
        samples = self.read_raw_samples(loadfields, **kwargs)
        # convert to FieldArray
        samples = array_class.from_kwargs(**samples)
        # add the static params and attributes
        addatrs = (list(self.static_params.items()) +
                   list(self[self.samples_group].attrs.items()))
        for (p, val) in addatrs:
            if p in loadfields:
                continue
            setattr(samples, format_attr(p), format_attr(val))
        return samples

    @abstractmethod
    def read_raw_samples(self, fields, **kwargs):
        """Low level function for reading datasets in the samples group.

        This should return a dictionary of numpy arrays.
        """
        pass

    @staticmethod
    def extra_args_parser(parser=None, skip_args=None, **kwargs):
        """Provides a parser that can be used to parse sampler-specific command
        line options for loading samples.

        This is optional. Inheriting classes may override this if they want to
        implement their own options.

        Parameters
        ----------
        parser : argparse.ArgumentParser, optional
            Instead of creating a parser, add arguments to the given one. If
            none provided, will create one.
        skip_args : list, optional
            Don't include the given options. Options should be given as the
            option string, minus the '--'. For example,
            ``skip_args=['iteration']`` would cause the ``--iteration``
            argument not to be included.
        \**kwargs :
            All other keyword arguments are passed to the parser that is
            created.

        Returns
        -------
        parser : argparse.ArgumentParser or None
            If this class adds extra arguments, an argument parser with the
            extra arguments. Otherwise, will just return whatever was passed
            for the ``parser`` argument (default is None).
        actions : list of argparse.Action
            List of the actions that were added.
        """
        return parser, []

    @staticmethod
    def _get_optional_args(args, opts, err_on_missing=False, **kwargs):
        """Convenience function to retrieve arguments from an argparse
        namespace.

        Parameters
        ----------
        args : list of str
            List of arguments to retreive.
        opts : argparse.namespace
            Namespace to retreive arguments for.
        err_on_missing : bool, optional
            If an argument is not found in the namespace, raise an
            AttributeError. Otherwise, just pass. Default is False.
        \**kwargs :
            All other keyword arguments are added to the return dictionary.
            Any keyword argument that is the same as an argument in ``args``
            will override what was retrieved from ``opts``.

        Returns
        -------
        dict :
            Dictionary mapping arguments to values retrieved from ``opts``. If
            keyword arguments were provided, these will also be included in the
            dictionary.
        """
        parsed = {}
        for arg in args:
            try:
                parsed[arg] = getattr(opts, arg)
            except AttributeError as e:
                if err_on_missing:
                    raise AttributeError(e)
                else:
                    continue
        parsed.update(kwargs)
        return parsed

    def samples_from_cli(self, opts, parameters=None, **kwargs):
        """Reads samples from the given command-line options.

        Parameters
        ----------
        opts : argparse Namespace
            The options with the settings to use for loading samples (the sort
            of thing returned by ``ArgumentParser().parse_args``).
        parameters : (list of) str, optional
            A list of the parameters to load. If none provided, will try to
            get the parameters to load from ``opts.parameters``.
        \**kwargs :
            All other keyword arguments are passed to ``read_samples``. These
            will override any options with the same name.

        Returns
        -------
        FieldArray :
            Array of the loaded samples.
        """
        if parameters is None and opts.parameters is None:
            parameters = self.variable_params
        elif parameters is None:
            parameters = opts.parameters
        # parse optional arguments
        _, extra_actions = self.extra_args_parser()
        extra_args = [act.dest for act in extra_actions]
        kwargs = self._get_optional_args(extra_args, opts, **kwargs)
        return self.read_samples(parameters, **kwargs)

    @property
    def static_params(self):
        """Returns a dictionary of the static_params. The keys are the argument
        names, values are the value they were set to.
        """
        return {arg: self.attrs[arg] for arg in self.attrs["static_params"]}

    @property
    def effective_nsamples(self):
        """Returns the effective number of samples stored in the file.
        """
        try:
            return self.attrs['effective_nsamples']
        except KeyError:
            return 0

    def write_effective_nsamples(self, effective_nsamples):
        """Writes the effective number of samples stored in the file."""
        self.attrs['effective_nsamples'] = effective_nsamples

    @property
    def thin_start(self):
        """The default start index to use when reading samples.

        Unless overridden by sub-class attribute, just returns 0.
        """
        return 0

    @property
    def thin_interval(self):
        """The default interval to use when reading samples.

        Unless overridden by sub-class attribute, just returns 1.
        """
        return 1

    @property
    def thin_end(self):
        """The defaut end index to use when reading samples.

        Unless overriden by sub-class attribute, just return ``None``.
        """
        return None

    @property
    def cmd(self):
        """Returns the (last) saved command line.

        If the file was created from a run that resumed from a checkpoint, only
        the last command line used is returned.

        Returns
        -------
        cmd : string
            The command line that created this InferenceFile.
        """
        cmd = self.attrs["cmd"]
        if isinstance(cmd, numpy.ndarray):
            cmd = cmd[-1]
        return cmd

    def write_logevidence(self, lnz, dlnz):
        """Writes the given log evidence and its error.

        Results are saved to file's 'log_evidence' and 'dlog_evidence'
        attributes.

        Parameters
        ----------
        lnz : float
            The log of the evidence.
        dlnz : float
            The error in the estimate of the log evidence.
        """
        self.attrs['log_evidence'] = lnz
        self.attrs['dlog_evidence'] = dlnz

    @property
    def log_evidence(self):
        """Returns the log of the evidence and its error, if they exist in the
        file. Raises a KeyError otherwise.
        """
        return self.attrs["log_evidence"], self.attrs["dlog_evidence"]

    def write_random_state(self, group=None, state=None):
        """Writes the state of the random number generator from the file.

        The random state is written to ``sampler_group``/random_state.

        Parameters
        ----------
        group : str
            Name of group to write random state to.
        state : tuple, optional
            Specify the random state to write. If None, will use
            ``numpy.random.get_state()``.
        """
        # Write out the default numpy random state
        group = self.sampler_group if group is None else group
        dataset_name = "/".join([group, "random_state"])
        if state is None:
            state = numpy.random.get_state()
        s, arr, pos, has_gauss, cached_gauss = state
        if dataset_name in self:
            self[dataset_name][:] = arr
        else:
            self.create_dataset(dataset_name, arr.shape, fletcher32=True,
                                dtype=arr.dtype)
            self[dataset_name][:] = arr
        self[dataset_name].attrs["s"] = s
        self[dataset_name].attrs["pos"] = pos
        self[dataset_name].attrs["has_gauss"] = has_gauss
        self[dataset_name].attrs["cached_gauss"] = cached_gauss

    def read_random_state(self, group=None):
        """Reads the state of the random number generator from the file.

        Parameters
        ----------
        group : str
            Name of group to read random state from.

        Returns
        -------
        tuple
            A tuple with 5 elements that can be passed to numpy.set_state.
        """
        # Read numpy randomstate
        group = self.sampler_group if group is None else group
        dataset_name = "/".join([group, "random_state"])
        arr = self[dataset_name][:]
        s = self[dataset_name].attrs["s"]
        pos = self[dataset_name].attrs["pos"]
        has_gauss = self[dataset_name].attrs["has_gauss"]
        cached_gauss = self[dataset_name].attrs["cached_gauss"]
        state = s, arr, pos, has_gauss, cached_gauss
        return state

    def write_strain(self, strain_dict, group=None):
        """Writes strain for each IFO to file.

        Parameters
        -----------
        strain : {dict, FrequencySeries}
            A dict of FrequencySeries where the key is the IFO.
        group : {None, str}
            The group to write the strain to. If None, will write to the top
            level.
        """
        subgroup = self.data_group + "/{ifo}/strain"
        if group is None:
            group = subgroup
        else:
            group = '/'.join([group, subgroup])
        for ifo, strain in strain_dict.items():
            self[group.format(ifo=ifo)] = strain
            self[group.format(ifo=ifo)].attrs['delta_t'] = strain.delta_t
            self[group.format(ifo=ifo)].attrs['start_time'] = \
                float(strain.start_time)

    def write_stilde(self, stilde_dict, group=None):
        """Writes stilde for each IFO to file.

        Parameters
        -----------
        stilde : {dict, FrequencySeries}
            A dict of FrequencySeries where the key is the IFO.
        group : {None, str}
            The group to write the strain to. If None, will write to the top
            level.
        """
        subgroup = self.data_group + "/{ifo}/stilde"
        if group is None:
            group = subgroup
        else:
            group = '/'.join([group, subgroup])
        for ifo, stilde in stilde_dict.items():
            self[group.format(ifo=ifo)] = stilde
            self[group.format(ifo=ifo)].attrs['delta_f'] = stilde.delta_f
            self[group.format(ifo=ifo)].attrs['epoch'] = float(stilde.epoch)

    def write_psd(self, psds, group=None):
        """Writes PSD for each IFO to file.

        PSDs are written to ``[{group}/]data/{detector}/psds/0``, where {group}
        is the optional keyword argument.

        Parameters
        -----------
        psds : dict
            A dict of detector name -> FrequencySeries.
        group : str, optional
            Specify a top-level group to write the data to. If ``None`` (the
            default), data will be written to the file's top level.
        """
        subgroup = self.data_group + "/{ifo}/psds/0"
        if group is None:
            group = subgroup
        else:
            group = '/'.join([group, subgroup])
        for ifo in psds:
            self[group.format(ifo=ifo)] = psds[ifo]
            self[group.format(ifo=ifo)].attrs['delta_f'] = psds[ifo].delta_f

    def write_injections(self, injection_file, group=None):
        """Writes injection parameters from the given injection file.

        Everything in the injection file is copied to
        ``[{group}/]injections_group``, where ``{group}`` is the optional
        keyword argument.

        Parameters
        ----------
        injection_file : str
            Path to HDF injection file.
        group : str, optional
            Specify a top-level group to write the injections group to. If
            ``None`` (the default), injections group will be written to the
            file's top level.
        """
        logging.info("Writing injection file to output")
        if group is None or group == '/':
            group = self.injections_group
        else:
            group = '/'.join([group, self.injections_group])
        try:
            with h5py.File(injection_file, "r") as fp:
                super(BaseInferenceFile, self).copy(fp, group)
        except IOError:
            logging.warn("Could not read %s as an HDF file", injection_file)

    def read_injections(self, group=None):
        """Gets injection parameters.

        Injections are retrieved from ``[{group}/]injections``.

        Parameters
        ----------
        group : str, optional
            Group that the injections group is in. Default (None) is to look
            in the top-level.

        Returns
        -------
        FieldArray
            Array of the injection parameters.
        """
        if group is None or group == '/':
            group = self.injections_group
        else:
            group = '/'.join([group, self.injections_group])
        injset = InjectionSet(self.filename, hdf_group=group)
        injections = injset.table.view(FieldArray)
        # close the new open filehandler to self
        injset._injhandler.filehandler.close()
        return injections

    def write_command_line(self):
        """Writes command line to attributes.

        The command line is written to the file's ``attrs['cmd']``. If this
        attribute already exists in the file (this can happen when resuming
        from a checkpoint), ``attrs['cmd']`` will be a list storing the current
        command line and all previous command lines.
        """
        cmd = [" ".join(sys.argv)]
        try:
            previous = self.attrs["cmd"]
            if isinstance(previous, str):
                # convert to list
                previous = [previous]
            elif isinstance(previous, numpy.ndarray):
                previous = previous.tolist()
        except KeyError:
            previous = []
        self.attrs["cmd"] = cmd + previous

    @staticmethod
    def get_slice(thin_start=None, thin_interval=None, thin_end=None):
        """Formats a slice to retrieve a thinned array from an HDF file.

        Parameters
        ----------
        thin_start : float or int, optional
            The starting index to use. If provided, the ``int`` will be taken.
        thin_interval : float or int, optional
            The interval to use. If provided the ceiling of it will be taken.
        thin_end : float or int, optional
            The end index to use. If provided, the ``int`` will be taken.

        Returns
        -------
        slice :
            The slice needed.
        """
        if thin_start is not None:
            thin_start = int(thin_start)
        if thin_interval is not None:
            thin_interval = int(numpy.ceil(thin_interval))
        if thin_end is not None:
            thin_end = int(thin_end)
        return slice(thin_start, thin_end, thin_interval)

    def copy_metadata(self, other):
        """Copies all metadata from this file to the other file.

        Metadata is defined as everything in the top-level ``.attrs``.

        Parameters
        ----------
        other : InferenceFile
            An open inference file to write the data to.
        """
        logging.info("Copying metadata")
        # copy attributes
        for key in self.attrs.keys():
            other.attrs[key] = self.attrs[key]

    def copy_info(self, other, ignore=None):
        """Copies "info" from this file to the other.

        "Info" is defined all groups that are not the samples group.

        Parameters
        ----------
        other : output file
            The output file. Must be an hdf file.
        ignore : (list of) str
            Don't copy the given groups.
        """
        logging.info("Copying info")
        # copy non-samples/stats data
        if ignore is None:
            ignore = []
        if isinstance(ignore, str):
            ignore = [ignore]
        ignore = set(ignore + [self.samples_group])
        copy_groups = set(self.keys()) - ignore
        for key in copy_groups:
            super(BaseInferenceFile, self).copy(key, other)

    def copy_samples(self, other, parameters=None, parameter_names=None,
                     read_args=None, write_args=None):
        """Should copy samples to the other files.

        Parameters
        ----------
        other : InferenceFile
            An open inference file to write to.
        parameters : list of str, optional
            List of parameters to copy. If None, will copy all parameters.
        parameter_names : dict, optional
            Rename one or more parameters to the given name. The dictionary
            should map parameter -> parameter name. If None, will just use the
            original parameter names.
        read_args : dict, optional
            Arguments to pass to ``read_samples``.
        write_args : dict, optional
            Arguments to pass to ``write_samples``.
        """
        # select the samples to copy
        logging.info("Reading samples to copy")
        if parameters is None:
            parameters = self.variable_params
        # if list of desired parameters is different, rename
        if set(parameters) != set(self.variable_params):
            other.attrs['variable_params'] = parameters
        if read_args is None:
            read_args = {}
        samples = self.read_samples(parameters, **read_args)
        logging.info("Copying {} samples".format(samples.size))
        # if different parameter names are desired, get them from the samples
        if parameter_names:
            arrs = {pname: samples[p] for p, pname in parameter_names.items()}
            arrs.update({p: samples[p] for p in parameters if
                         p not in parameter_names})
            samples = FieldArray.from_kwargs(**arrs)
            other.attrs['variable_params'] = samples.fieldnames
        logging.info("Writing samples")
        if write_args is None:
            write_args = {}
        other.write_samples({p: samples[p] for p in samples.fieldnames},
                            **write_args)

    def copy(self, other, ignore=None, parameters=None, parameter_names=None,
             read_args=None, write_args=None):
        """Copies metadata, info, and samples in this file to another file.

        Parameters
        ----------
        other : str or InferenceFile
            The file to write to. May be either a string giving a filename,
            or an open hdf file. If the former, the file will be opened with
            the write attribute (note that if a file already exists with that
            name, it will be deleted).
        ignore : (list of) strings
            Don't copy the given groups. If the samples group is included, no
            samples will be copied.
        parameters : list of str, optional
            List of parameters in the samples group to copy. If None, will copy
            all parameters.
        parameter_names : dict, optional
            Rename one or more parameters to the given name. The dictionary
            should map parameter -> parameter name. If None, will just use the
            original parameter names.
        read_args : dict, optional
            Arguments to pass to ``read_samples``.
        write_args : dict, optional
            Arguments to pass to ``write_samples``.

        Returns
        -------
        InferenceFile
            The open file handler to other.
        """
        if not isinstance(other, h5py.File):
            # check that we're not trying to overwrite this file
            if other == self.name:
                raise IOError("destination is the same as this file")
            other = self.__class__(other, 'w')
        # metadata
        self.copy_metadata(other)
        # info
        if ignore is None:
            ignore = []
        if isinstance(ignore, str):
            ignore = [ignore]
        self.copy_info(other, ignore=ignore)
        # samples
        if self.samples_group not in ignore:
            self.copy_samples(other, parameters=parameters,
                              parameter_names=parameter_names,
                              read_args=read_args,
                              write_args=write_args)
            # if any down selection was done, re-set the default
            # thin-start/interval/end
            p = tuple(self[self.samples_group].keys())[0]
            my_shape = self[self.samples_group][p].shape
            p = tuple(other[other.samples_group].keys())[0]
            other_shape = other[other.samples_group][p].shape
            if my_shape != other_shape:
                other.attrs['thin_start'] = 0
                other.attrs['thin_interval'] = 1
                other.attrs['thin_end'] = None
        return other

    @classmethod
    def write_kwargs_to_attrs(cls, attrs, **kwargs):
        """Writes the given keywords to the given ``attrs``.

        If any keyword argument points to a dict, the keyword will point to a
        list of the dict's keys. Each key is then written to the attrs with its
        corresponding value.

        Parameters
        ----------
        attrs : an HDF attrs
            The ``attrs`` of an hdf file or a group in an hdf file.
        \**kwargs :
            The keywords to write.
        """
        for arg, val in kwargs.items():
            if val is None:
                val = str(None)
            if isinstance(val, dict):
                attrs[str(arg)] = list(map(str, val.keys()))
                # just call self again with the dict as kwargs
                cls.write_kwargs_to_attrs(attrs, **val)
            else:
                attrs[str(arg)] = val

    def write_data(self, name, data, path=None, append=False):
        """Convenience function to write data.

        Given ``data`` is written as a dataset with ``name`` in ``path``.
        If the dataset or path do not exist yet, the dataset and path will
        be created.

        Parameters
        ----------
        name : str
            The name to associate with the data. This will be the dataset
            name (if data is array-like) or the key in the attrs.
        data : array, dict, or atomic
            The data to write. If a dictionary, a subgroup will be created
            for each key, and the values written there. This will be done
            recursively until an array or atomic (e.g., float, int, str), is
            found. Otherwise, the data is written to the given name.
        path : str, optional
            Write to the given path. Default (None) will write to the top
            level. If the path does not exist in the file, it will be
            created.
        append : bool, optional
            Append the data to what is currently in the file if ``path/name``
            already exists in the file, and if it does not, create the dataset
            so that its last dimension can be resized. The data can only
            be appended along the last dimension, and if it already exists in
            the data, it must be resizable along this dimension. If ``False``
            (the default) what is in the file will be overwritten, and the
            given data must have the same shape.
        """
        if path is None:
            path = '/'
        try:
            group = self[path]
        except KeyError:
            # create the group
            self.create_group(path)
            group = self[path]
        if isinstance(data, dict):
            # call myself for each key, value pair in the dictionary
            for key, val in data.items():
                self.write_data(key, val, path='/'.join([path, name]),
                                append=append)
        # if appending, we need to resize the data on disk, or, if it doesn't
        # exist yet, create a dataset that is resizable along the last
        # dimension
        elif append:
            # cast the data to an array if it isn't already one
            if isinstance(data, (list, tuple)):
                data = numpy.array(data)
            if not isinstance(data, numpy.ndarray):
                data = numpy.array([data])
            dshape = data.shape
            ndata = dshape[-1]
            try:
                startidx = group[name].shape[-1]
                group[name].resize(dshape[-1]+group[name].shape[-1],
                                   axis=len(group[name].shape)-1)
            except KeyError:
                # dataset doesn't exist yet
                group.create_dataset(name, dshape,
                                     maxshape=tuple(list(dshape)[:-1]+[None]),
                                     dtype=data.dtype, fletcher32=True)
                startidx = 0
            group[name][..., startidx:startidx+ndata] = data[..., :]
        else:
            try:
                group[name][()] = data
            except KeyError:
                # dataset doesn't exist yet
                group[name] = data

    def write_config_file(self, cp):
        """Writes the given config file parser.

        File is stored as a pickled buffer array to ``config_parser/{index}``,
        where ``{index}`` is an integer corresponding to the number of config
        files that have been saved. The first time a save is called, it is
        stored to ``0``, and incremented from there.

        Parameters
        ----------
        cp : ConfigParser
            Config parser to save.
        """
        # get the index of the last saved file
        try:
            index = list(map(int, self[self.config_group].keys()))
        except KeyError:
            index = []
        if index == []:
            # hasn't been written yet
            index = 0
        else:
            index = max(index) + 1
        # we'll store the config file as a text file that is pickled
        out = StringIO()
        cp.write(out)
        # now pickle it
        dump_state(out, self, path=self.config_group, dsetname=str(index))

    def read_config_file(self, return_cp=True, index=-1):
        """Reads the config file that was used.

        A ``ValueError`` is raised if no config files have been saved, or if
        the requested index larger than the number of stored config files.

        Parameters
        ----------
        return_cp : bool, optional
            If true, returns the loaded config file as
            :py:class:`pycbc.workflow.configuration.WorkflowConfigParser`
            type. Otherwise will return as string buffer. Default is True.
        index : int, optional
            The config file to load. If ``write_config_file`` has been called
            multiple times (as would happen if restarting from a checkpoint),
            there will be config files stored. Default (-1) is to load the
            last saved file.

        Returns
        -------
        WorkflowConfigParser or StringIO :
            The parsed config file.
        """
        # get the stored indices
        try:
            indices = sorted(map(int, self[self.config_group].keys()))
            index = indices[index]
        except KeyError:
            raise ValueError("no config files saved in hdf")
        except IndexError:
            raise ValueError("no config file matches requested index")
        cf = load_state(self, path=self.config_group, dsetname=str(index))
        cf.seek(0)
        if return_cp:
            cp = WorkflowConfigParser()
            cp.read_file(cf)
            return cp
        return cf

    def read_data(self, group=None):
        """Loads the data stored in the file as a FrequencySeries.

        Only works for models that store data as a frequency series in
        ``data/DET/stilde``. A ``KeyError`` will be raised if the model used
        did not store data in that path.

        Parameters
        ----------
        group : str, optional
            Group that the data group is in. Default (None) is to look in the
            top-level.

        Returns
        -------
        dict :
            Dictionary of detector name -> FrequencySeries.
        """
        fmt = '{}/{}/stilde'
        if group is None or group == '/':
            path = self.data_group
        else:
            path = '/'.join([group, self.data_group])
        data = {}
        for det in self[path].keys():
            group = self[fmt.format(path, det)]
            data[det] = FrequencySeries(
                group[()], delta_f=group.attrs['delta_f'],
                epoch=group.attrs['epoch'])
        return data

    def read_psds(self, group=None):
        """Loads the PSDs stored in the file as a FrequencySeries.

        Only works for models that store PSDs in
        ``data/DET/psds/0``. A ``KeyError`` will be raised if the model used
        did not store PSDs in that path.

        Parameters
        ----------
        group : str, optional
            Group that the data group is in. Default (None) is to look in the
            top-level.

        Returns
        -------
        dict :
            Dictionary of detector name -> FrequencySeries.
        """
        fmt = '{}/{}/psds/0'
        if group is None or group == '/':
            path = self.data_group
        else:
            path = '/'.join([group, self.data_group])
        psds = {}
        for det in self[path].keys():
            group = self[fmt.format(path, det)]
            psds[det] = FrequencySeries(
                group[()], delta_f=group.attrs['delta_f'])
        return psds
