# Copyright (C) 2019 Collin Capano, Sumit Kumar
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
"""Provides IO for the dynesty sampler.
"""

import argparse
import numpy
from pycbc.io.hdf import (dump_state, load_state)

from .base_nested_sampler import BaseNestedSamplerFile
from .posterior import write_samples_to_file, read_raw_samples_from_file

class CommonNestedMetadataIO(object):
    """Provides functions for reading/writing dynesty metadata to file.
    """

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

        if 'raw_samples' not in skip_args:
            act = parser.add_argument(
                "--raw-samples", action='store_true', default=False,
                help="Extract raw samples rather than a posterior. "
                     "Raw samples are the unweighted samples obtained from "
                     "the nested sampler. Default value is False, which means "
                     "raw samples are weighted by the log-weight array "
                     "obtained from the sampler, giving an estimate of the "
                     "posterior.")
            actions.append(act)
        if 'seed' not in skip_args:
            act = parser.add_argument(
                "--seed", type=int, default=0,
                help="Set the random-number seed used for extracting the "
                     "posterior samples. This is needed because the "
                     "unweighted samples are randomly shuffled to produce "
                     "a posterior. Default is 0. Ignored if raw-samples are "
                     "extracted instead.")
        return parser, actions

class DynestyFile(CommonNestedMetadataIO, BaseNestedSamplerFile):
    """Class to handle file IO for the ``dynesty`` sampler."""

    name = 'dynesty_file'

    def read_raw_samples(self, fields, raw_samples=False, seed=0):
        """Reads samples from a dynesty file and constructs a posterior.

        Parameters
        ----------
        fields : list of str
            The names of the parameters to load. Names must correspond to
            dataset names in the file's ``samples`` group.
        raw_samples : bool, optional
            Return the raw (unweighted) samples instead of the estimated
            posterior samples. Default is False.
        seed : int, optional
            When extracting the posterior, samples are randomly shuffled. To
            make this reproduceable, numpy's random generator seed is set with
            the given value prior to the extraction. Default is 0.

        Returns
        -------
        dict :
            Dictionary of parameter names -> samples.
        """
        samples = read_raw_samples_from_file(self, fields)
        logwt = read_raw_samples_from_file(self, ['logwt'])['logwt']
        loglikelihood = read_raw_samples_from_file(
            self, ['loglikelihood'])['loglikelihood']
        logz = self.attrs.get('log_evidence')
        if not raw_samples:
            weights = numpy.exp(logwt - logz)
            N = len(weights)
            positions = (numpy.random.random() + numpy.arange(N)) / N
            idx = numpy.zeros(N, dtype=int)
            cumulative_sum = numpy.cumsum(weights)
            cumulative_sum /= cumulative_sum[-1]
            i, j = 0, 0
            while i < N:
                if positions[i] < cumulative_sum[j]:
                    idx[i] = j
                    i += 1
                else:
                    j += 1
            try:
                rng = numpy.random.default_rng(seed)
            except AttributeError:
                # numpy pre-1.17 uses RandomState
                # Py27: delete this after we drop python 2.7 support
                rng = numpy.random.RandomState(seed)
            rng.shuffle(idx)
            post = {'loglikelihood': loglikelihood[idx]}
            for i, param in enumerate(fields):
                post[param] = samples[param][idx]
            return post
        else:
            return samples

    def write_pickled_data_into_checkpoint_file(self, state):
        """Dump the sampler state into checkpoint file
        """
        if 'sampler_info/saved_state' not in self:
            self.create_group('sampler_info/saved_state')
        dump_state(state, self, path='sampler_info/saved_state')

    def read_pickled_data_from_checkpoint_file(self):
        """Load the sampler state (pickled) from checkpoint file
        """
        return load_state(self, path='sampler_info/saved_state')

    def write_raw_samples(self, data, parameters=None):
        """Write the nested samples to the file
        """
        if 'samples' not in self:
            self.create_group('samples')
        write_samples_to_file(self, data, parameters=parameters,
                              group='samples')

    def validate(self):
        """Runs a validation test.
        This checks that a samples group exist, and that pickeled data can
        be loaded.

        Returns
        -------
        bool :
            Whether or not the file is valid as a checkpoint file.
        """
        try:
            if 'sampler_info/saved_state' in self:
                load_state(self, path='sampler_info/saved_state')
            checkpoint_valid = True
        except KeyError:
            checkpoint_valid = False
        return checkpoint_valid
