# Copyright (C) 2016 Christopher M. Biwer
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

import os
import sys
import h5py
import numpy
import logging
from pycbc import DYN_RANGE_FAC
from pycbc.types import FrequencySeries
from pycbc.waveform import parameters as wfparams
import pycbc.inference.sampler
import pycbc.inference.likelihood
from pycbc.io import FieldArray

class _PosteriorOnlyParser(object):
    """Provides interface for reading/writing samples from/to an InferenceFile
    that contains flattened posterior samples.
    """
    @staticmethod
    def _read_fields(fp, fields_group, fields, array_class,
                     thin_start=None, thin_interval=None, thin_end=None,
                     iteration=None):
        """Reads fields from the given file.
        """
        if iteration is not None:
            get_index = iteration
        else:
            get_index = fp.get_slice(thin_start=thin_start, thin_end=thin_end,
                                     thin_interval=thin_interval)
        # load
        arrays = {}
        group = fields_group + '/{}'
        arrays = {field: fp[group.format(field)][get_index]
                  for field in fields}
        return array_class.from_kwargs(**arrays)

    @classmethod
    def read_samples(cls, fp, parameters, samples_group=None,
                     thin_start=0, thin_end=None, thin_interval=1,
                     iteration=None, array_class=None):
        """Reads posterior samples from a posterior-only file.
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
    def write_samples_group(fp, samples_group, fields, samples):
        """Writes the given samples to the given samples group.
        """
        for field in samples.fieldnames:
            grp = '{}/{}'.format(samples_group, field)
            fp[grp] = samples[field]

    @classmethod
    def n_independent_samples(cls, fp):
        """Returns the number of independent samples stored in the file.
        """
        return cls.read_samples(fp, fp.variable_args[0]).size


class InferenceFile(h5py.File):
    """ A subclass of the h5py.File object that has extra functions for
    handling reading and writing the samples from the samplers.

    Parameters
    -----------
    path : str
        The path to the HDF file.
    mode : {None, str}
        The mode to open the file, eg. "w" for write and "r" for read.
    """
    name = "hdf"
    samples_group = 'samples'
    stats_group = 'likelihood_stats'
    sampler_group = 'sampler_states'

    def __init__(self, path, mode=None, **kwargs):
        super(InferenceFile, self).__init__(path, mode, **kwargs)

    @property
    def posterior_only(self):
        """Whether the file only contains flattened posterior samples.
        """
        try:
            return self.attrs['posterior_only']
        except KeyError:
            return False

    @property
    def sampler_name(self):
        """Returns the name of the sampler that was used."""
        return self.attrs["sampler"]

    @property
    def sampler_class(self):
        """Returns the sampler class that was used."""
        try:
            sampler = self.sampler_name
        except KeyError:
            return None
        return pycbc.inference.sampler.samplers[sampler]

    @property
    def samples_parser(self):
        """Returns the class to use to read/write samples from/to the file."""
        if self.posterior_only:
            return _PosteriorOnlyParser
        else:
            return self.sampler_class

    @property
    def likelihood_eval_name(self):
        """Returns the name of the likelihood evaluator that was used."""
        return self.attrs["likelihood_evaluator"]

    @property
    def variable_args(self):
        """Returns list of variable_args.

        Returns
        -------
        variable_args : {list, str}
            List of str that contain variable_args keys.
        """
        return self.attrs["variable_args"]

    @property
    def static_args(self):
        """Returns a dictionary of the static_args. The keys are the argument
        names, values are the value they were set to.
        """
        return dict([[arg, self.attrs[arg]]
            for arg in self.attrs["static_args"]])

    @property
    def sampling_args(self):
        """Returns the parameters that were used to sample.

        Returns
        -------
        sampling_args : {list, str}
            List of the sampling args.
        """
        return self.attrs["sampling_args"]

    @property
    def lognl(self):
        """Returns the log noise likelihood."""
        return self.attrs["lognl"]

    @property
    def niterations(self):
        """Returns number of iterations performed.

        Returns
        -------
        niterations : int
            Number of iterations performed.
        """
        return self.attrs["niterations"]

    @property
    def n_independent_samples(self):
        """Returns the number of independent samples stored in the file.
        """
        return self.samples_parser.n_independent_samples(self)

    @property
    def burn_in_iterations(self):
        """Returns number of iterations in the burn in.
        """
        return self.attrs["burn_in_iterations"]

    @property
    def is_burned_in(self):
        """Returns whether or not the sampler is burned in.
        """
        return self.attrs["is_burned_in"]

    @property
    def nwalkers(self):
        """Returns number of walkers used.

        Returns
        -------
        nwalkesr : int
            Number of walkers used.
        """
        return self.attrs["nwalkers"]

    @property
    def ntemps(self):
        """Returns number of temperatures used."""
        return self.attrs["ntemps"]

    @property
    def acl(self):
        """ Returns the saved autocorelation length (ACL).

        Returns
        -------
        acl : {int, float}
            The ACL.
        """
        return self.attrs["acl"]

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

    @property
    def resume_points(self):
        """The iterations at which a run was resumed from checkpoint.

        Returns
        -------
        resume_points : array or None
            An array of integers giving the points at which the run resumed.

        Raises
        ------
        KeyError
            If the run never resumed from a checkpoint.
        """
        return self.attrs['resume_points']

    @property
    def log_evidence(self):
        """Returns the log of the evidence and its error, if they exist in the
        file. Raises a KeyError otherwise.
        """
        return self.attrs["log_evidence"], self.attrs["dlog_evidence"]

    def read_samples(self, parameters, samples_group=None, **kwargs):
        """Reads samples from the file.

        Parameters
        -----------
        parameters : (list of) strings
            The parameter(s) to retrieve. A parameter can be the name of any
            field in `samples_group`, a virtual field or method of
            `FieldArray` (as long as the file contains the necessary fields
            to derive the virtual field or method), and/or a function of
            these.
        samples_group : str
            Group in HDF InferenceFile that parameters belong to.
        \**kwargs :
            The rest of the keyword args are passed to the sampler's
            `read_samples` method.

        Returns
        -------
        FieldArray
            Samples for the given parameters, as an instance of a
            FieldArray.
        """
        # get the appropriate sampler class
        samples_group = samples_group if samples_group else self.samples_group
        return self.samples_parser.read_samples(self, parameters,
                                                samples_group=samples_group,
                                                **kwargs)

    def read_likelihood_stats(self, **kwargs):
        """Reads likelihood stats from self.

        Parameters
        -----------
        \**kwargs :
            The keyword args are passed to the sampler's `read_likelihood_stats`
            method.

        Returns
        -------
        stats : {FieldArray, None}
            Likelihood stats in the file, as a FieldArray. The fields of the
            array are the names of the stats that are in the `likelihood_stats`
            group.
        """
        parameters = self[self.stats_group].keys()
        return self.read_samples(parameters, samples_group=self.stats_group,
                                 **kwargs)

    def read_acceptance_fraction(self, **kwargs):
        """Returns the acceptance fraction that was written to the file.

        Parameters
        ----------
        \**kwargs :
            All keyword arguments are passed to the sampler's
            `read_acceptance_fraction` function.
        Returns
        -------
        numpy.array
            The acceptance fraction.
        """
        return self.sampler_class.read_acceptance_fraction(self, **kwargs)

    def read_acls(self):
        """Returns all of the individual chains' acls. See the `read_acls`
        function of this file's sampler for more details.
        """
        return self.sampler_class.read_acls(self)

    def read_label(self, parameter, error_on_none=False):
        """Returns the label for the parameter.

        Parameters
        -----------
        parameter : str
            Name of parameter to get a label for. Will first try to retrieve
            a label from this file's "label" attributes. If the parameter
            is not found there, will look for a label from
            pycbc.waveform.parameters.
        error_on_none : {False, bool}
            If True, will raise a ValueError if a label cannot be found, or if
            the label is None. Otherwise, the parameter will just be returned
            if no label can be found.

        Returns
        -------
        label : str
            A formatted string for the name of the paramter.
        """
        # get label
        try:
            label = self[parameter].attrs["label"]
        except KeyError:
            # try looking in pycbc.waveform.parameters
            try:
                label = getattr(wfparams, parameter).label
            except AttributeError:
                label = None
        if label is None:
            if error_on_none:
                raise ValueError("Cannot find a label for paramter %s" %(
                    parameter))
            else:
                return parameter
        return label

    def read_random_state(self, group=None):
        """ Reads the state of the random number generator from the file.

        Parameters
        ----------
        group : str
            Name of group to read random state from.

        Returns
        -------
        tuple
            A tuple with 5 elements that can be passed to numpy.set_state.
        """
        group = self.sampler_group if group is None else group
        dataset_name = "/".join([group, "random_state"])
        arr = self[dataset_name][:]
        s = self[dataset_name].attrs["s"]
        pos = self[dataset_name].attrs["pos"]
        has_gauss = self[dataset_name].attrs["has_gauss"]
        cached_gauss = self[dataset_name].attrs["cached_gauss"]
        return s, arr, pos, has_gauss, cached_gauss

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
        subgroup = "{ifo}/strain"
        if group is None:
            group = subgroup
        else:
            group = '/'.join([group, subgroup])
        for ifo,strain in strain_dict.items():
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
        subgroup = "{ifo}/stilde"
        if group is None:
            group = subgroup
        else:
            group = '/'.join([group, subgroup])
        for ifo,stilde in stilde_dict.items():
            self[group.format(ifo=ifo)] = stilde
            self[group.format(ifo=ifo)].attrs['delta_f'] = stilde.delta_f
            self[group.format(ifo=ifo)].attrs['epoch'] = float(stilde.epoch)

    def write_psd(self, psds, low_frequency_cutoff, group=None):
        """Writes PSD for each IFO to file.

        Parameters
        -----------
        psds : {dict, FrequencySeries}
            A dict of FrequencySeries where the key is the IFO.
        low_frequency_cutoff : {dict, float}
            A dict of the low-frequency cutoff where the key is the IFO. The
            minimum value will be stored as an attr in the File.
        group : {None, str}
            The group to write the strain to. If None, will write to the top
            level.
        """
        subgroup = "{ifo}/psds/0"
        if group is None:
            group = subgroup
        else:
            group = '/'.join([group, subgroup])
        self.attrs["low_frequency_cutoff"] = min(low_frequency_cutoff.values())
        for ifo in psds:
            self[group.format(ifo=ifo)] = psds[ifo]
            self[group.format(ifo=ifo)].attrs['delta_f'] = psds[ifo].delta_f

    def write_data(self, strain_dict=None, stilde_dict=None,
                   psd_dict=None, low_frequency_cutoff_dict=None,
                   group=None):
        """Writes the strain/stilde/psd.

        Parameters
        ----------
        strain_dict : {None, dict}
            A dictionary of strains. If None, no strain will be written.
        stilde_dict : {None, dict}
            A dictionary of stilde. If None, no stilde will be written.
        psd_dict : {None, dict}
            A dictionary of psds. If None, no psds will be written.
        low_freuency_cutoff_dict : {None, dict}
            A dictionary of low frequency cutoffs used for each detector in
            `psd_dict`; must be provided if `psd_dict` is not None.
        group : {None, str}
            The group to write the strain to. If None, will write to the top
            level.
        """
        # save PSD
        if psd_dict is not None:
            if low_frequency_cutoff_dict is None:
                raise ValueError("must provide low_frequency_cutoff_dict if "
                                 "saving psds to output")
            # apply dynamic range factor for saving PSDs since
            # plotting code expects it
            psd_dyn_dict = {}
            for key,val in psd_dict.iteritems():
                psd_dyn_dict[key] = FrequencySeries(val*DYN_RANGE_FAC**2,
                                                    delta_f=val.delta_f)
            self.write_psd(psds=psd_dyn_dict,
                           low_frequency_cutoff=low_frequency_cutoff_dict,
                           group=group)

        # save stilde
        if stilde_dict is not None:
            self.write_stilde(stilde_dict, group=group)

        # save strain if desired
        if strain_dict is not None:
            self.write_strain(strain_dict, group=group)

    def write_injections(self, injection_file, ifo):
        """ Writes injection parameters for an IFO to file.

        Parameters
        ----------
        injection_file : str
            Path to HDF injection file.
        ifo : str
            IFO name.
        """
        subgroup = "{ifo}/injections"
        self.create_group(subgroup.format(ifo=ifo))
        try:
            with h5py.File(injection_file, "r") as fp:
                for param in fp.keys():
                    self[subgroup.format(ifo=ifo)][param] = fp[param][:]
                for key in fp.attrs.keys():
                    self[subgroup.format(ifo=ifo)].attrs[key] = fp.attrs[key]
        except IOError:
            logging.warn("Could not read %s as an HDF file", injection_file)

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

    def write_random_state(self, group=None, state=None):
        """ Writes the state of the random number generator from the file.

        Parameters
        ----------
        group : str
            Name of group to read random state to.
        state : tuple, optional
            Specify the random state to write. If None, will use
            ``numpy.random.get_state()``.
        """
        group = self.sampler_group if group is None else group
        dataset_name = "/".join([group, "random_state"])
        if state is None:
            state = numpy.random.get_state()
        s, arr, pos, has_gauss, cached_gauss = state
        if group in self:
            self[dataset_name][:] = arr
        else:
            self.create_dataset(dataset_name, arr.shape, fletcher32=True,
                                dtype=arr.dtype)
            self[dataset_name][:] = arr
        self[dataset_name].attrs["s"] = s
        self[dataset_name].attrs["pos"] = pos
        self[dataset_name].attrs["has_gauss"] = has_gauss
        self[dataset_name].attrs["cached_gauss"] = cached_gauss

    def get_slice(self, thin_start=None, thin_interval=None, thin_end=None):
        """Formats a slice using the given arguments that can be used to
        retrieve a thinned array from an InferenceFile.

        Parameters
        ----------
        thin_start : {None, int}
            The starting index to use. If None, will try to retrieve the
            `burn_in_iterations` from the given file. If no
            `burn_in_iterations` exists, will default to the start of the
            array.
        thin_interval : {None, int}
            The interval to use. If None, will try to retrieve the acl from the
            given file. If no acl attribute exists, will default to 1.
        thin_end : {None, int}
            The end index to use. If None, will retrieve to the end of the
            array.

        Returns
        -------
        slice :
            The slice needed.
        """

        # default is to skip burn in samples
        if thin_start is None:
            try:
                thin_start = self.burn_in_iterations
                # if the sampler hasn't burned in, the burn_in_iterations will
                # be the same as the number of iterations, which would result
                # in 0 samples. In that case, just use the last one
                if thin_start == self.niterations:
                    thin_start = thin_start - 1
            except KeyError:
                pass

        # default is to use stored ACL and accept every i-th sample
        if thin_interval is None:
            try:
                thin_interval = int(numpy.ceil(self.acl))
            except KeyError:
                pass
        return slice(thin_start, thin_end, thin_interval)

    def copy_metadata(self, other):
        """Copies all metadata from this file to the other file.

        Metadata is defined as all data that is not in either the samples or
        stats group.

        Parameters
        ----------
        other : InferenceFile
            An open inference file to write the data to.
        """
        logging.info("Copying metadata")
        # copy non-samples/stats data
        for key in self.keys():
            if key not in [self.samples_group, self.stats_group]:
                super(InferenceFile, self).copy(key, other)
        # copy attributes
        for key in self.attrs.keys():
            other.attrs[key] = self.attrs[key]


    def copy(self, other, parameters=None, parameter_names=None,
             posterior_only=False, **kwargs):
        """Copies data in this file to another file.

        The samples and stats to copy may be down selected using the given
        kwargs. All other data (the "metadata") are copied exactly.

        Parameters
        ----------
        other : str or InferenceFile
            The file to write to. May be either a string giving a filename,
            or an open hdf file. If the former, the file will be opened with
            the write attribute (note that if a file already exists with that
            name, it will be deleted).
        parameters : list of str, optional
            List of parameters to copy. If None, will copy all parameters.
        parameter_names : dict, optional
            Rename one or more parameters to the given name. The dictionary
            should map parameter -> parameter name. If None, will just use the
            original parameter names.
        posterior_only : bool, optional
            Write the samples and likelihood stats as flattened arrays, and
            set other's posterior_only attribute. For example, if this file
            has a parameter's samples written to
            `{samples_group}/{param}/walker{x}`, then other will have all of
            the selected samples from all walkers written to
            `{samples_group}/{param}/`.
        \**kwargs :
            All other keyword arguments are passed to `read_samples`.

        Returns
        -------
        InferenceFile
            The open file handler to other.
        """
        if not isinstance(other, h5py.File):
            # check that we're not trying to overwrite this file
            if other == self.name:
                raise IOError("destination is the same as this file")
            other = InferenceFile(other, 'w')
        # copy metadata over
        self.copy_metadata(other)
        # update other's posterior attribute
        if posterior_only:
            other.attrs['posterior_only'] = posterior_only
        # select the samples to copy
        logging.info("Reading samples to copy")
        if parameters is None:
            parameters = self.variable_args
        # if list of desired parameters is different, rename variable args
        if set(parameters) != set(self.variable_args):
            other.attrs['variable_args'] = parameters
        # if only the posterior is desired, we'll flatten the results
        if not posterior_only and not self.posterior_only:
            kwargs['flatten'] = False
        samples = self.read_samples(parameters, **kwargs)
        logging.info("Copying {} samples".format(samples.size))
        # if different parameter names are desired, get them from the samples
        if parameter_names:
            arrs = {pname: samples[p] for p,pname in parameter_names.items()}
            arrs.update({p: samples[p] for p in parameters
                                        if p not in parameter_names})
            samples = FieldArray.from_kwargs(**arrs)
            other.attrs['variable_args'] = samples.fieldnames
        logging.info("Writing samples")
        other.samples_parser.write_samples_group(other, self.samples_group,
                                                 samples.fieldnames, samples)
        # do the same for the likelihood stats
        logging.info("Reading stats to copy")
        stats = self.read_likelihood_stats(**kwargs)
        logging.info("Writing stats")
        other.samples_parser.write_samples_group(other, self.stats_group,
                                                 stats.fieldnames, stats)
        # if any down selection was done, re-set the burn in iterations and
        # the acl, and the niterations.
        # The last dimension of the samples returned by the sampler should
        # be the number of iterations.
        if samples.shape[-1] != self.niterations:
            other.attrs['acl'] = 1
            other.attrs['burn_in_iterations'] = 0
            other.attrs['niterations'] = samples.shape[-1]
        return other


def check_integrity(filename):
    """Checks the integrity of an InferenceFile.

    Checks done are:

        * can the file open?
        * do all of the datasets in the samples group have the same shape?
        * can the first and last sample in all of the datasets in the samples
          group be read?

    If any of these checks fail, an IOError is raised.

    Parameters
    ----------
    filename: str
        Name of an InferenceFile to check.

    Raises
    ------
    ValueError
        If the given file does not exist.
    KeyError
        If the samples group does not exist.
    IOError
        If any of the checks fail.
    """
    # check that the file exists
    if not os.path.exists(filename):
        raise ValueError("file {} does not exist".format(filename))
    # if the file is corrupted such that it cannot be opened, the next line
    # will raise an IOError
    with InferenceFile(filename, 'r') as fp:
        # check that all datasets in samples have the same shape
        parameters = fp[fp.samples_group].keys()
        group = fp.samples_group + '/{}'
        # use the first parameter as a reference shape
        ref_shape = fp[group.format(parameters[0])].shape
        if not all(fp[group.format(param)].shape == ref_shape
                   for param in parameters):
            raise IOError("not all datasets in the samples group have the same "
                          "shape")
        # check that we can read the first/last sample
        firstidx = tuple([0]*len(ref_shape))
        lastidx = tuple([-1]*len(ref_shape))
        for param in parameters:
            _ = fp[group.format(param)][firstidx]
            _ = fp[group.format(param)][lastidx]
