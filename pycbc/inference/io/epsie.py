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

from __future__ import (absolute_import, division)

from pickle import UnpicklingError
from epsie import load_state

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

    @property
    def swap_interval(self):
        """The interval that temperature swaps occurred at."""
        return self[self.sampler_group].attrs['swap_interval']

    @swap_interval.setter
    def swap_interval(self, swap_interval):
        """Stores the swap interval to the sampler group's attrs."""
        self[self.sampler_group].attrs['swap_interval'] = swap_interval

    def write_sampler_metadata(self, sampler):
        """Adds writing seed and betas to MultiTemperedMCMCIO.
        """
        super(EpsieFile, self).write_sampler_metadata(sampler)
        self[self.sampler_group].attrs['seed'] = sampler.seed
        try:
            self[self.sampler_group]["betas"][:] = sampler.betas
        except KeyError:
            self[self.sampler_group]["betas"] = sampler.betas

    def thin(self, thin_interval):
        """Thins the samples on disk to the given thinning interval.

        Also thins the acceptance ratio and the temperature data, both of
        which are stored in the ``sampler_info`` group.
        """
        # thin the samples
        super(EpsieFile, self).thin(thin_interval)
        # thin the acceptance ratio
        new_interval = thin_interval // self.thinned_by
        self._thin_data(self.sampler_group, ['acceptance_ratio'],
                        new_interval)
        # thin the temperature swaps; since these may not happen every
        # iteration, the thin interval we use for these is different
        ts_group = '/'.join([self.sampler_group, 'temperature_swaps'])
        ts_thin_interval = new_interval // self.swap_interval
        if ts_thin_interval > 1:
            self._thin_data(ts_group, ['swap_index', 'acceptance_ratio'],
                            ts_thin_interval)

    def write_acceptance_ratio(self, acceptance_ratio, last_iteration=None):
        """Writes the acceptance ratios to the sampler info group.

        Parameters
        ----------
        acceptance_ratio : array
            The acceptance ratios to write. Should have shape
            ``ntemps x nchains x niterations``.
        """
        # we'll use the write_samples machinery to write the acceptance ratios
        self.write_samples({'acceptance_ratio': acceptance_ratio},
                           last_iteration=last_iteration,
                           samples_group=self.sampler_group)

    def write_temperature_data(self, swap_index, acceptance_ratio,
                               swap_interval, last_iteration):
        """Writes temperature swaps and acceptance ratios.

        Parameters
        ----------
        swap_index : array
            The indices indicating which temperatures were swapped. Should have
            shape ``ntemps x nchains x (niterations/swap_interval)``.
        acceptance_ratio : array
            The array of acceptance ratios between temperatures. Should
            have shape ``(ntemps-1) x nchains x (niterations/swap_interval)``.
            arrays.
        swap_interval : int
            The number of iterations between temperature swaps.
        last_iteration : int
            The iteration of the last sample.
        """
        self.swap_interval = swap_interval
        group = '/'.join([self.sampler_group, 'temperature_swaps'])
        # we'll use the write_samples machinery to write the acceptance ratios;
        # if temperature swaps didn't happen every iteration, then a smaller
        # thinning interval than what is used for the samples should be used
        thin_by = self.thinned_by // swap_interval
        # we'll also tell the write samples that the last "iteration" is the
        # last iteration / the swap interval, to get the spacing correct
        last_iteration = last_iteration // swap_interval
        # we need to write the two arrays separately, since they have different
        # dimensions in temperature
        self.write_samples({'swap_index': swap_index},
                           last_iteration=last_iteration,
                           samples_group=group, thin_by=thin_by)
        self.write_samples({'acceptance_ratio': acceptance_ratio},
                           last_iteration=last_iteration,
                           samples_group=group, thin_by=thin_by)

    def validate(self):
        """Adds attemp to load checkpoint to validation test."""
        valid = super(EpsieFile, self).validate()
        # try to load the checkpoint
        if valid:
            try:
                load_state(self, self.sampler_group)
            except (KeyError, UnpicklingError):
                # will get this if the state wasn't written, or it was
                # corrupted for some reason
                valid = False
        return valid
