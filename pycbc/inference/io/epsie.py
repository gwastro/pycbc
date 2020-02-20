# Copyright (C) 2019  Collin Capano
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

"""This module provides IO classes for epsie samplers.
"""

from __future__ import absolute_import

from epsie.samplers import load_state

from .base_sampler import BaseSamplerFile
from .base_multitemper import (MultiTemperedMCMCIO, MultiTemperedMetadataIO)

class EpsieFile(MultiTemperedMCMCIO, MultiTemperedMetadataIO,
                BaseSamplerFile):
    """Class to handle IO for Epsie's parallel-tempered sampler."""
    
    name = 'epsie_file'

    @property
    def betas(self):
        """The betas that were used."""
        return self[self.sampler_group]['betas'][()]

    def write_sampler_metadata(self, sampler):
        """Adds writing betas to MultiTemperedMCMCIO.
        """
        super(EpsieFile, self).write_sampler_metadata(sampler)
        try:
            self[self.sampler_group]["betas"][:] = sampler.betas
        except KeyError:
            self[self.sampler_group]["betas"] = sampler.betas

    def write_acceptance_ratio(self, acceptance_ratios, last_iteration=None):
        """Writes the acceptance ratios to the sampler info group.
        
        Parameters
        ----------
        acceptance_ratios : array
            The acceptance ratios to write. Should be a have shape
            ``ntemps x nchains x niterations``.
        """
        # we'll use the write_samples machinery to write the acceptance ratios
        self.write_samples({'acceptance_ratios': acceptance_ratios},
                           last_iteration=last_iteration,
                           group=self.sampler_group)

    def write_temperature_data(self, temperature_data, last_iteration=None):
        """Writes temperature swaps and acceptance ratios.

        Parameters
        ----------
        temperature_data : dict
            Dictionary mapping ``'acceptance_ratio'`` and ``'swap_index'`` to
            arrays.
        """
        group = '/'.join([self.sampler_group, 'temperature_swaps'])
        self.write_samples(temperature_data, last_iteration=last_iteration,
                           group=group)

    def validate(self):
        """Adds looking for checkpoint group to validation test."""
        valid = super(EpsieFile, self).validate()
        return valid
        if valid:
            valid = self.state_path in self
        # try to load the checkpoint
        if valid:
            try:
                load_state(self, self.state_path)
            except:
                valid = False
        return valid
