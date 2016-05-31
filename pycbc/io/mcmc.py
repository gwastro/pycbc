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
"""This modules defines functions for reading and writing samples that the MCMC
samplers generate.
"""

import h5py
import numpy

class MCMCFile(h5py.File):
    """ A subclass of the h5py.File object that has extra functions for
    handling reading and writing the samples from the MCMC samplers.

    Parameters
    -----------
    path : str
        The path to the HDF file.
    mode : {None, str}
        The mode to open the file, eg. "w" for write and "r" for read.
    """

    def __init__(self, path, mode=None, **kwargs):
        super(MCMCFile, self).__init__(path, mode, **kwargs)

    def read_samples(self, variable_arg, thin_start=0, thin_interval=1):
        """ Reads independent samples from all walkers for a parameter.

        Parameters
        -----------
        variable_arg : str
            Name of parameter to get independent samples.
        thin_start : int
            Index of the sample to begin returning samples.
        thin_interval : int
            Interval to accept every i-th sample.

        Returns
        -------
        numpy.array
            All independent samples from all walkers for a parameter.
        """

        nwalkers = self.attrs["nwalkers"]
        return numpy.array([self.read_samples_from_walker(variable_arg, j, thin_start, thin_interval) for j in range(nwalkers)])

    def read_samples_from_walker(self, variable_arg, nwalker,
                                 thin_start=0, thin_interval=1):
        """ Reads all samples from a specific walker for a parameter.

        Parameters
        -----------
        variable_arg : str
            Name of parameter to get independent samples.
        nwalker : int
            Index of the walker to get samples.
        thin_start : int
            Index of the sample to begin returning samples.
        thin_interval : int
            Interval to accept every i-th sample.

        Returns
        -------
        numpy.array
            Samples from a specific walker for a parameter.
        """
        return self[variable_arg]["walker%d"%nwalker][thin_start::thin_interval]

    def read_label(self, variable_arg):
        """ Returns the label for the parameter.

        Parameters
        -----------
        variable_arg : str
            Name of parameter to get label.

        Returns
        -------
        label : str
            A formatted string for the name of the paramter.
        """
        return self[variable_arg].attrs["label"]

    def write(self, variable_args, ifo_list, samples, acceptance_fraction=None,
              labels=None, low_frequency_cutoff=None, psds=None):
        """ Writes the output from pycbc.io.sampler to a file.

        Parameters
        -----------
        variable_args : list
        ifo_list : list
        samples : numpy.Array
            An array with shape (ndim,nwalker,niterations) where ndim is the
            number of dimensions, nwalker is the number of walkers, and
            niterations is the number of iterations.
        acceptance_fraction : numpy.Array
            A 1-dimensional array that contains the number of samples accepted
            for each iteration.
        labels : list
            A list of str that have formatted names for parameter.
        low_frequency_cutoff : dict
            The low-frequency cutoff values each IFO PSD.
        psds : dict
            A dict with the IFO name as the key and a FreqeuncySeries as the
            value.
        """

        # get number of dimensions, walkers, and iterations
        ndim, nwalkers, niterations = samples.shape

        # save MCMC parameters
        self.attrs["variable_args"] = variable_args
        self.attrs["ifo_list"] = ifo_list
        self.attrs["nwalkers"] = nwalkers
        self.attrs["niterations"] = niterations
        self.attrs["low_frequency_cutoff"] = min(low_frequency_cutoff.values())

        # loop over number of dimensions
        for i,dim_name in zip(range(ndim), variable_args):

            # create a group in the output file for this dimension
            group_dim = self.create_group(dim_name)

            # add label
            if labels:
                group_dim.attrs["label"] = labels[i]
            else:
                group_dim.attrs["label"] = dim_name

            # loop over number of walkers
            for j in range(nwalkers):

                # get all samples from this walker for this dimension
                samples_subset = samples[i,j,:]

                # write to output file
                dataset_name = "walker%d"%j
                group_dim.create_dataset(dataset_name, data=samples_subset)

        # create a dataset for the acceptance fraction
        if acceptance_fraction is not None:
            self.create_dataset("acceptance_fraction",
                                data=acceptance_fraction)

        # create datasets for each PSD
        if psds and low_frequency_cutoff:
            for key in psds.keys():
                psd_dim = self.create_dataset(key+"/psds/0",
                                              data=psds[key])
                psd_dim.attrs["delta_f"] = psds[key].delta_f


