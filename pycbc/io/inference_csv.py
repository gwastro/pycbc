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

"""This modules defines functions for reading and samples that the
inference samplers generate and are stored in CSV file.
"""

import numpy
from pycbc.io import record
from pycbc.waveform import parameters as wfparams

COMMENT = ""
DELIMITER = " "

class InferenceCSVFile(object):
    """ A class that has extra functions for handling reading the samples
     from posterior-only CSV files.

    Parameters
    -----------
    path : str
        The path to the HDF file.
    mode : {None, str}
        The mode to open the file. Only accepts "r" or "rb" for reading.
    delimiter : str
        Delimiter to use for CSV file. Default is space-delimited.
    """

    def __init__(self, path, mode=None, delimiter=DELIMITER):
        self.path = path
        self.delimiter = delimiter
        if mode in ["r", "rb"]:
            self.mode = mode
        else:
            raise ValueError("Mode for InferenceCSVFile must be 'r' or 'rb'.")

        with open(self.path, self.mode) as fp:
            header = fp.readline().rstrip("\n")
        self._variable_args = header.split(self.delimiter)

    @property
    def variable_args(self):
        """ List of parameters from CSV file.
        """
        return self._variable_args

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
        try:
            label = getattr(wfparams, parameter).label
        except AttributeError:
            if error_on_none:
                raise ValueError(
                            "Cannot find a label for paramter %s" %(parameter))
            else:
                return parameter
        return label

    def read_samples(self, parameters, **kwargs):
        """Reads posterior samples from a posterior-only CSV file.
        """
        idxs = [self.variable_args.index(p) for p in parameters]
        samples = numpy.loadtxt(self.path, delimiter=self.delimiter,
                                comments=COMMENT, skiprows=1, usecols=idxs)
        dtype = zip(parameters, len(idxs) * [float])
        size = len(samples) if len(idxs) > 1 else samples.shape[0]
        rec_samples = record.FieldArray(size, dtype=dtype)
        if len(idxs) > 1:
            for i, p in zip(idxs, parameters):
                rec_samples[p] = samples[:, i]
        else:
            rec_samples[parameters[0]] = samples
        return rec_samples

    def read_likelihood_stats(self, parameters=None, **kwargs):
        """Reads likelihood stats from self.

        Parameters
        -----------
        parameters : list
            A list of parameters to return. Default is `None` which returns all
            `likelihood_stats` parameters.
        \**kwargs :
            The keyword args are passed to the `read_samples` method.

        Returns
        -------
        stats : {FieldArray, None}
            Likelihood stats in the file, as a FieldArray. The fields of the
            array are the names of the stats that are in the `likelihood_stats`
            group.
        """
        parameters = self.variable_args if parameters is None else parameters
        return self.read_samples(parameters, **kwargs)

    def close(self):
        """ A dummy function so that `InferenceCSVFile` has analogous method
        as `InferenceFile`.
        """
        pass
