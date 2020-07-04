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
from pycbc.io.hdf import (dump_state, load_state)
from .base_nested_sampler import BaseNestedSamplerFile

class DynestyFile(BaseNestedSamplerFile):
    """Class to handle file IO for the ``dynesty`` sampler."""

    name = 'dynesty_file'

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
