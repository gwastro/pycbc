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
        return self._variable_args

    def read_label(self, parameter, error_on_none=False):
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
        parameters = self.variable_args if parameters is None else parameters
        return self.read_samples(parameters)

    def close(self):
        pass
