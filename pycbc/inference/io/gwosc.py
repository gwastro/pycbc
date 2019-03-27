# Copyright (C) 2019 Collin Capano
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

"""Provides features to read GWOSC posterior files.
"""

from __future__ import absolute_import

import numpy
import h5py
import argparse

from pycbc.io import FieldArray
from .base_hdf import get_optional_args


class GWOSCPosteriorFile(h5py.File):
    """A class to be able to read and write GWOSC-style posterior HDF files.
    """
    name = 'gwosc_posterior_file'
    samples_dataset = 'Overall_posterior'

    def __init__(self, path, mode=None, **kwargs):
        super(GWOSCPosteriorFile, self).__init__(path, mode, **kwargs)
        if mode == 'w':
            # first time creating the file, add this class's name
            self.attrs['filetype'] = self.name

    def params_in_dataset(self, dataset=None):
        """Gets the parameters in the given dataset.
        
        Parameters
        ----------
        dataset : str, optional
            The dataset to read from. Default (None) is to read from
            ``self.samples_dataset``.
        """
        if dataset is None:
            dataset = self.samples_dataset
        return self[dataset].dtype.names

    def collection(self, collection):
        """A list of parameter names that are in the given collection.
        
        Parameters
        ----------
        collection : {'all', 'variable_params'}
            Name of collection of parameters to get.

        Returns
        -------
        list of str :
            List of parameter names that are in the given collection.
        """
        if collection == 'all':
            return self.all_params
        elif collection == 'variable_params':
            return self.variable_params
        else:
            raise ValueError("unknown parameter collection {}"
                             .format(collection))

    @property
    def variable_params(self):
        """Gets the parameters in the default samples dataset."""
        return self.params_in_dataset()

    @property
    def all_params(self):
        """Gets a set of all of the parameters across all datasets."""
        return list(set([p for d in self.keys()
                         for p in self.params_in_dataset(d)]))

    def write_samples(self, samples, samples_dataset=None):
        """Writes the given samples.
        
        Samples are combined into a single structured array before being
        written to file.

        Parameters
        ----------
        samples : dict
            Dictionary mapping parameter names to numpy arrays.
        samples_dataset : str, optional
            Which dataset to write the samples to. If ``None`` (the default),
            will write to "Overall_posterior".
        """
        # make sure all arrays have the same shape
        shapes = set([d.shape for d in samples.values()])
        if len(shapes) != 1:
            raise ValueError('all arrays must have the same shape')
        # combine samples into a structured array
        dtype = [(p, d.dtype) for (p, d) in samples.items()]
        data = numpy.zeros(shapes[0], dtype=dtype)
        for p in samples:
            data[p] = samples[p]
        if samples_dataset is None:
            samples_dataset = self.samples_dataset
        self[samples_dataset] = data

    def read_samples(self, parameters=None, samples_dataset=None):
        """Loads samples as a ``FieldArray``.

        Parameters
        ----------
        parameters : list of str, optional
            Only return the given parameters. Default is to load all of the
            parameters stored.
        samples_dataset : str, optional
            Which dataset in the file from which to read the samples. Default
            is ``"Overall_posterior"``.

        Returns
        -------
        FieldArray :
            The samples.
        """
        if samples_dataset is None:
            samples_dataset = self.samples_dataset
        samples = self[samples_dataset][:].view(type=FieldArray)
        # drop unwanted parameters
        if parameters is not None:
            if isinstance(parameters, (str, unicode)):
                parameters = [parameters]
            # check that we have the requested parameters
            unknown = set(parameters) - set(samples.fieldnames)
            if any(unknown):
                raise ValueError("unknown parameters {}"
                                 .format(', '.join(unknown)))
            samples = samples[parameters].view(type=FieldArray)
        return samples

    @classmethod
    def extra_args_parser(cls, parser=None, skip_args=None, **kwargs):
        """Adds ``--samples-dataset`` to an argument parser.
        
        Parameters
        ----------
        parser : argparse.ArgumentParser, optional
            Instead of creating a parser, add arguments to the given one. If
            none provided, will create one.
        skip_args : list, optional
            Don't include the given options. Options should be given as the
            option string, minus the '--'.
        \**kwargs :
            All other keyword arguments are passed to the parser that is
            created.

        Returns
        -------
        parser : argparse.ArgumentParser
            An argument parser with the extra arguments added.
        actions : list of argparse.Action
            A list of the actions that were added.
        """
        if parser is None:
            parser = argparse.ArgumentParser(**kwargs)
        if skip_args is None:
            skip_args = []
        actions = []
        if 'samples_dataset' not in skip_args:
            act = parser.add_argument(
                "--samples-dataset", type=str, default=cls.samples_dataset,
                help="Which dataset to read samples from. Default is {}."
                     .format(cls.samples_dataset))
            actions.append(act)
        return parser, actions

    def parse_parameters(self, parameters, array_class=None, dataset=None):
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
        possible_fields = self.params_in_dataset(dataset)
        return array_class.parse_parameters(parameters, possible_fields)

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
        # parse optional arguments
        _, extra_actions = self.extra_args_parser()
        extra_args = [act.dest for act in extra_actions]
        kwargs = get_optional_args(extra_args, opts, **kwargs)
        if parameters is None and opts.parameters is None:
            try:
                samples_dataset = kwargs['samples_dataset']
            except KeyError:
                samples_dataset = self.samples_dataset
            parameters = self.params_in_dataset(samples_dataset)
        elif parameters is None:
            parameters = opts.parameters
        parameters = self.parse_parameters(parameters)
        return self.read_samples(parameters, **kwargs)
