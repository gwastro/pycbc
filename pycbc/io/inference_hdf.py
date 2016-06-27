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
from pycbc import pnutils
from pycbc.results import str_utils
from pycbc.io.record import WaveformArray
from pycbc.waveform import parameters as wfparams

def read_label_from_config(cp, variable_arg, section="labels", html=False):
    """ Returns the label for the variable_arg.

    Parameters
    ----------
    cp : WorkflowConfigParser
        A WorkflowConfigParser instance with [labels] section
    variable_arg : str
        The parameter to get label.
    section : str
        Name of section in configuration file to get label.
    html : bool
        If True then replace LaTeX substrings with HTML substrings.

    Returns
    -------
    label : str
        The label for the parameter.
    """

    # get label from configuration file if it exists
    if cp.has_option(section, variable_arg):
        label = cp.get(section, variable_arg)
    else:
        # try looking in pycbc.waveform.parameters
        try:
            label = getattr(wfparams, varable_arg).label
        except AttributeError:
            # just use the parameter name
            label = variable_arg

    # replace LaTeX with HTML
    if html:
        label = str_utils.latex_to_html(label)

    return label

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

    def __init__(self, path, mode=None, **kwargs):
        super(InferenceFile, self).__init__(path, mode, **kwargs)
        # create a result class to return values as. This is a sub-class
        # WaveformArray, with the _staticfields set to the variable args
        # that are in the file. We can only do this if there are actual
        # results in the file
        try:
            self._arraycls = self._create_arraycls()
        except KeyError:
            self._arraycls = None

    def _create_arraycls(self):
        """Returns a sub-class of WaveformArray, with the _staticfields
        set to the variable args that are in the file.
        """
        # we'll need the name of a walker to get the dtypes
        refwalker = self[self.variable_args[0]].keys()[0]
        # get the names, dtypes of the variable args
        fields = dict([[name, self[name][refwalker].dtype]
            for name in self.variable_args])
        class ResultArray(WaveformArray):
            _staticfields = fields
        return ResultArray

    @property
    def variable_args(self):
        """ Returns list of variable_args.

        Returns
        -------
        variable_args : {list, str}
            List of str that contain variable_args keys.
        """
        return self.attrs["variable_args"]

    @property
    def nwalkers(self):
        """ Returns number of walkers used.

        Returns
        -------
        nwalkesr : int
            Number of walkers used.
        """
        return self.attrs["nwalkers"]

    @property
    def niterations(self):
        """ Returns number of iterations performed.

        Returns
        -------
        niterations : int
            Number of iterations performed.
        """
        return self.attrs["niterations"]

    @property
    def acl(self):
        """ Returns the saved autocorelation length (ACL).

        Returns
        -------
        acl : {int, float}
            The ACL.
        """
        return self.attrs["acl"]

    def read_samples_from_walkers(self, parameters, walkers=None,
                                 thin_start=None, thin_interval=None):
        """Reads samples from the specified walker(s) for the given
        parameter(s).

        Parameters
        -----------
        parameters : (list of) strings
            The parameter(s) to retrieve. A parameter can be the name of a
            variable argument, a virtual field or method of WaveformArray
            (as long as the file contains the necessary variables to derive
            the virtual field or method), and/or a function of these.
        walkers : {None, (list of) int}
            The walker index (or a list of indices) to retrieve. If None,
            samples from all walkers will be obtained.
        thin_start : int
            Index of the sample to begin returning samples. Default is to read
            samples after burn in. To start from the beginning set thin_start
            to 0.
        thin_interval : int
            Interval to accept every i-th sample. Default is to use the
            self.acl attribute. If self.acl is not set, then use all samples
            (set thin_interval to 1).

        Returns
        -------
        WaveformArray
            Samples for the given parameters, as an instance of a
            WaveformArray.
        """
        if walkers is None:
            walkers = range(self.nwalkers)
        if isinstance(walkers, int):
            walkers = [walkers]
        # default is to skip burn in samples
        thin_start = self.attrs["burn_in_iterations"] if thin_start is None \
            else thin_start

        # default is to use stored ACL and accept every i-th sample
        if "acl" in self.attrs.keys():
            thin_interval = self.acl if thin_interval is None \
                else thin_interval
        else:
            thin_interval = 1 if thin_interval is None else thin_interval

        # figure out the size of the output array to create
        n_per_walker = \
            self[self.variable_args[0]]['walker0'][thin_start::thin_interval].size
        arrsize = len(walkers) * n_per_walker

        # create an array to store the results
        arr = self._arraycls(arrsize, names=parameters)
        # populate
        for name in arr.fieldnames:
            for ii,walker in enumerate(walkers):
                arr[name][ii*n_per_walker:(ii+1)*n_per_walker] = \
                    self[name]['walker%i' % walker][thin_start::thin_interval]
        return arr

    def read_samples(self, parameters, thin_start=None, thin_interval=None):
        """ Reads samples from all of the walkers for the given
        parameter(s). See read_samples_from_walkers for more details.
        """
        return self.read_samples_from_walkers(parameters, walkers=None,
            thin_start=thin_start, thin_interval=thin_interval)

    def read_acceptance_fraction(self, thin_start=None, thin_interval=None):
        """ Returns a numpy.array of the fraction of samples acceptanced at
        each iteration in the sampler..

        Parameters
        -----------
        thin_start : int
            Index of the sample to begin returning samples.
        thin_interval : int
            Interval to accept every i-th sample.

        Returns
        -------
        numpy.array
            The acceptance fraction.
        """
        return self["acceptance_fraction"][thin_start::thin_interval]

    def read_label(self, parameter, html=False, error_on_none=False):
        """Returns the label for the parameter.

        Parameters
        -----------
        parameter : str
            Name of parameter to get a label for. Will first try to retrieve
            a label from this file's "label" attributes. If the parameter
            is not found there, will look for a label from
            pycbc.waveform.parameters.
        html : bool
            If true then escape LaTeX formatting for HTML rendering.
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
        if label == None:
            if error_on_none:
                raise ValueError("Cannot find a label for paramter %s" %(
                    parameter))
            else:
                return parameter
        # replace LaTeX with HTML
        if html:
            label = str_utils.latex_to_html(label)

        return label

    def write_psds(self, psds, low_frequency_cutoff):
        """ Writes PSD for each IFO to file.

        Parameters
        -----------
        psds : {dict, FrequencySeries}
            A dict of FrequencySeries where the key is the IFO.
        low_frequency_cutoff : {dict, float}
            A dict of the low-frequency cutoff where the key is the IFO. The
            minimum value will be stored as an attr in the File.
        """
        self.attrs["low_frequency_cutoff"] = min(low_frequency_cutoff.values())
        for key in psds.keys():
            psd_dim = self.create_dataset(key+"/psds/0",
                                          data=psds[key])
            psd_dim.attrs["delta_f"] = psds[key].delta_f

    def write_sampler_attrs(self, sampler):
        """ Write information about the sampler to the file attrs. If sampler
        will run burn in, this should be called after sampler has already
        run burn in.

        Parameters
        -----------
        sampler : pycbc.inference._BaseSampler
            An instance of a sampler class from pycbc.inference.sampler.
        """

        # get attributes from waveform generator
        variable_args = sampler.likelihood_evaluator.waveform_generator.variable_args
        ifo_list = sampler.likelihood_evaluator.waveform_generator.detector_names

        # get shape of data from a (niterations,nwalkers,ndim) array
        niterations, nwalkers, _ = sampler.chain.shape

        # save parameters
        self.attrs["variable_args"] = variable_args
        self.attrs["ifo_list"] = ifo_list
        self.attrs["nwalkers"] = nwalkers
        self.attrs["burn_in_iterations"] = sampler.burn_in_iterations

        # save number of iterations so far
        if "niterations" not in self.attrs.keys():
            self.attrs["niterations"] = niterations

    def write_samples(self, variable_args, data=None, start=None, end=None,
                      nwalkers=0, niterations=0, labels=None):
        """ Writes samples to the file. To write to a subsample of the array
        use start and end. To initialize an empty array of length niterations
        for each walker use nwalkers and niterations.

        Parameters
        -----------
        variable_args : {list, str}
            List of names of parameters.
        data : numpy.array
            Data to be saved with shape (niterations,nwalkers,ndim). if data is
            None then create an array fo zeros for each walker with
            length niterations.
        start : int
            If given then begin inserting this data at this index.
        end : int
            If given then stop inserting data at this index.
        nwalkers : int
            Number of walkers should be given if data is None.
        niterations : int
            Number of iterations should be given if data is None.
        labels : {list, str}
            A list of str to use as dislay by downsteam executables.
        """

        # transpose past samples to get an (ndim,nwalkers,niteration) array
        if data is not None:
            samples = numpy.transpose(data)
            ndim, nwalkers, niterations = samples.shape

        # sanity check options
        elif nwalkers == 0 and niterations != 0:
            raise ValueError("If nwalkers is 0 then niterations must be 0")

        # if no data is given then initialize to array of numpy.NAN
        # with shape (ndim,nwalkers,niterations)
        else:
            ndim = len(variable_args)
            shape = (ndim, nwalkers, niterations)
            samples = numpy.zeros(shape)

        # save number of iterations so far
        if "niterations" in self.attrs.keys():
            self.attrs["niterations"] += niterations - self.attrs["niterations"]
        else:
            self.attrs["niterations"] = niterations

        # loop over number of dimensions
        for i,dim_name in zip(range(ndim), variable_args):

            # create a group in the output file for this dimension
            if dim_name not in self.keys():
                group_dim = self.create_group(dim_name)
                if "label" in group_dim.attrs.keys():
                    pass
                elif labels != None:
                    group_dim.attrs["label"] = labels[i]
                else:
                    group_dim.attrs["label"] = dim_name

            # loop over number of walkers
            for j in range(nwalkers):

                # create dataset with shape (ndim,nwalkers,niterations)
                dataset_name = "walker%d"%j
                if dataset_name not in self[dim_name].keys():
                    samples_subset = numpy.zeros(niterations)
                    if data is not None:
                        samples_subset[start:end] = samples[i,j,start:end]
                    group_dim.create_dataset(dataset_name,
                                             data=samples_subset)

                # write all samples in range from walker for this dimension
                else:
                    if end > len(samples[i,j,:]):
                        end = None
                    samples_subset = samples[i,j,start:end]
                    self[dim_name+"/"+dataset_name][start:end] = samples_subset

    def write_acceptance_fraction(self, data=None, start=None, end=None,
                                  niterations=None):
        """ Write acceptance_fraction data to file. To write to a subsample
        of the array use start and end. To initialize an array of zeros, set
        data to None and specify niterations.

        Parameters
        -----------
        data : numpy.array
            Data to be saved with shape (niterations,nwalkers,ndim). if data is
            None then create an array fo zeros for each walker with
            length niterations.
        start : int
            If given then begin inserting this data at this index.
        end : int
            If given then stop inserting data at this index.
        niterations : int
            Number of iterations should be given if data is None.
        """

        # sanity check options and if data is not given then make empty array
        if data is None and niterations is None:
            raise ValueError("Must specify either data or niterations")
        elif data is None:
            data = numpy.zeros(niterations)

        # write data
        if "acceptance_fraction" not in self.keys():
            self.create_dataset("acceptance_fraction",
                                data=data)
        else:
            self["acceptance_fraction"][start:end] = data[start:end]

    def write_samples_from_sampler(self, sampler, start=None, end=None,
                      nwalkers=0, niterations=0, labels=None):
        """ Write data from sampler to file.

        Parameters
        -----------
        sampler : pycbc.inference._BaseSampler
            An instance of a sampler class from pycbc.inference.sampler.
        start : int
            If given then begin inserting this data at this index.
        end : int
            If given then stop inserting data at this index.
        nwalkers : int
            Number of walkers should be given if data is None.
        niterations : int
            Number of iterations should be given if data is None.
        labels : {list, str}
            A list of str to use as dislay by downsteam executables.
        """
        variable_args = sampler.likelihood_evaluator.waveform_generator.variable_args
        self.write_samples(variable_args, data=sampler.chain,
                           start=start, end=end,
                           nwalkers=nwalkers, niterations=niterations,
                           labels=labels)

    def write_acl(self, acl):
        """ Writes the autocorrelation length (ACL) to file.

        Parameters
        ----------
        acl : float
            The ACL.
        """
        self.attrs["acl"] = acl
