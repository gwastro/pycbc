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

from .base_sampler import BaseSamplerFile

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
        super(EmceePTFile, self).write_sampler_metadata(sampler)
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

    def write_state(self, state):
        """Writes the given state to file.

        Parameters
        ----------
        state : dict
            Dictionary (of possibly other dicts) giving all of the information
            needed to restart the epsie sampler.
        """
        group = '/'.join([self.sampler_group, 'state'])
        write_recursive_dict(self, group, state)

    def read_state(self):
        """Reads the sampler state from the given file.

        Returns
        -------
        dict :
            The state as a dictionary. This can be passed to the sampler to
            set its state.
        """
        group = '/'.join([self.sampler_group, 'state'])
        return read_recursive_dict(self, group)


def write_recursive_dict(fp, group, write_dict):
    """Writes a dictionary to the given group.

    If an element in the given dict is itself a dict, the element will be
    written as a sub-group. This continues until all sub dictionaries have
    been written.

    Parameters
    ----------
    fp : h5py.File instance
        An hdf file to write the dictionary to.
    group : str
        Name of the top-level group to write to.
    write_dict : dict
        The (possibly nested) dictionary to write.
    """
    for key, val in write_dict.items():
        writeto = '/'.join([group, key])
        if isinstance(val, dict):
            # call myself with the new sub group
            write_recursive_dict(fp, writeto, val)
        else:
            fp[writeto] = val


def read_recursive_dict(fp, group):
    """Creates a dictionary of all datasets in a given group.

    If the group contains sub-groups, the information in that group is returned
    as a sub-dictionary.

    Parameters
    ----------
    fp : h5py.File instance
        An hdf file to read the dictionary from.
    group : str
        Name of the top-level group to read from.

    Returns
    -------
    dict :
        Dictionary representing all elements in the group.
    """
    out = {}
    for key in fp[group].keys():
        readfrom = '/'.join([group, key])
        if isinstance(fp[readfrom], h5py.Group):
            # call myself to get the sub info
            out[key] = read_recursive_dict(fp, readfrom)
        else:
            out[key] = fp[readfrom][()]
    return out
