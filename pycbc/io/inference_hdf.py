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

import h5py
import numpy
from pycbc import DYN_RANGE_FAC
from pycbc.types import FrequencySeries
from pycbc.waveform import parameters as wfparams
import pycbc.inference.sampler
import pycbc.inference.likelihood
import sys

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
    samples_group = 'samples'
    stats_group = 'likelihood_stats'
    sampler_group = 'sampler_states'

    def __init__(self, path, mode=None, **kwargs):
        super(InferenceFile, self).__init__(path, mode, **kwargs)

    @property
    def sampler_name(self):
        """Returns the name of the sampler that was used."""
        return self.attrs["sampler"]

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
    def burn_in_iterations(self):
        """Returns number of iterations in the burn in.
        """
        return self.attrs["burn_in_iterations"]

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
        """ Returns the saved command line.

        Returns
        -------
        cmd : {str}
            The command line that created this InferenceFile.
        """
        return self.attrs["cmd"]

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
        samples_group = samples_group if samples_group \
                             else InferenceFile.samples_group
        sclass = pycbc.inference.sampler.samplers[self.sampler_name]
        return sclass.read_samples(self, parameters,
                                   samples_group=samples_group, **kwargs)

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
        return self.read_samples(
                          parameters, samples_group=self.stats_group, **kwargs)

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
        # get the appropriate sampler class
        sclass = pycbc.inference.sampler.samplers[self.sampler_name]
        return sclass.read_acceptance_fraction(self, **kwargs)

    def read_acls(self):
        """Returns all of the individual chains' acls. See the `read_acls`
        function of this file's sampler for more details.
        """
        # get the appropriate sampler class
        sclass = pycbc.inference.sampler.samplers[self.sampler_name]
        return sclass.read_acls(self)

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
            A dictionary of psds. If None, psds will be written.
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


    def write_command_line(self):
        """Writes command line to attributes.

        Parameters
        ----------
        opts : argparse.ArgumentParser
            The parsed command line instance.
        """
        self.attrs["cmd"] = " ".join(sys.argv)

    def write_random_state(self, group=None):
        """ Writes the state of the random number generator from the file.

        Parameters
        ----------
        group : str
            Name of group to read random state to.
        """
        group = self.sampler_group if group is None else group
        dataset_name = "/".join([group, "random_state"])
        s, arr, pos, has_gauss, cached_gauss = numpy.random.get_state()
        if group in self:
            self[dataset_name][:] = arr
        else:
            self[dataset_name] = arr
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
            except KeyError:
                pass

        # default is to use stored ACL and accept every i-th sample
        if thin_interval is None:
            try:
                thin_interval = int(numpy.ceil(self.acl))
            except KeyError:
                pass
        return slice(thin_start, thin_end, thin_interval)
