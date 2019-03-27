# Copyright (C) 2017 Christopher M. Biwer, Collin Capano
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
""" This modules defines functions for reading and samples that the
inference samplers generate and are stored in an ASCII TXT file.
"""

import numpy
import os
from pycbc.io import FieldArray


class InferenceTXTFile(object):
    """A class that has extra functions for handling reading the samples
    from posterior-only TXT files.

    Parameters
    -----------
    path : str
        The path to the TXT file.
    mode : {None, str}
        The mode to open the file. Only accepts "r" or "rb" for reading.
    delimiter : str
        Delimiter to use for TXT file. Default is space-delimited.
    """
    name = "txt"
    comments = ""
    delimiter = " "

    def __init__(self, path, mode=None, delimiter=None):
        self.path = path
        self.delimiter = delimiter if delimiter is not None else self.delimiter
        self.mode = mode

    def close(self):
        """Dummy function to make this class more like an HDF file."""
        pass

    @property
    def variable_params(self):
        """The variable parameters in the file.

        If the first line of the file starts with ``#``, then the parameter
        names are assumed to be given there as ``delimiter``-separated strings.

        If no comment string exists (i.e., the first line does not start with
        '#'), then will just return a list of 'pX', where X enumerates the
        columns.
        """
        with open(self.path, 'r') as fp:
            firstline = fp.readline().strip().rstrip('\n')
        if firstline.startswith('#'):
            variable_params = firstline.lstrip('#').strip().split(
                self.delimiter)
        else:
            nparams = len(firstline.split(self.delimiter))
            variable_params = ['p{}'.format(ii) for ii in range(nparams)]
        return variable_params

    def read_raw_samples(self, parameters=None):
        """Loads samples as a dictionary of arrays."""
        all_params = self.variable_params
        if parameters is None:
            parameters = all_params
        # load the file
        data = numpy.loadtxt(self.path)
        samples = {}
        for (pi, param) in enumerate(all_params):
            if param in parameters:
                samples[param] = data[:, pi]
        return samples

    def read_samples(self, parameters=None):
        """Loads samples as a ``FieldArray``."""
        return FieldArray.from_kwargs(**self.read_raw_samples(parameters))

    def write_samples(self, samples):
        """Writes the given samples to ``path``.

        Parameters
        ----------
        samples : dict
            Dictionary of numpy arrays to write.
        """
        params = samples.keys()
        numpy.savetxt(self.path, [samples[p] for p in params],
                      header=self.delimiter.join(params),
                      delimiter=self.delimiter)

    @staticmethod
    def extra_args_parser(parser=None, **kwargs):
        """Not used for this class."""
        return parser, []

    def parse_parameters(self, parameters, array_class=None):
        """Parses a parameters arg to figure out what fields need to be loaded.

        Parameters
        ----------
        parameters : (list of) strings
            The parameter(s) to retrieve. A parameter can be the name of any
            field in the ``variable_params``, and/or a function of these.
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
        possible_fields = self.variable_params
        return array_class.parse_parameters(parameters, possible_fields)

    def samples_from_cli(self, opts, parameters=None):
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
            parameters = self.variable_args
        elif parameters is None:
            parameters = opts.parameters
        # parse optional arguments
        return self.read_samples(parameters)
