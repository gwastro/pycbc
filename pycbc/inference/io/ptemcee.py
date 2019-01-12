# Copyright (C) 2018 Collin Capano
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


"""Provides I/O support for ptemcee.
"""

from __future__ import absolute_import
from .base_hdf import BaseInferenceFile
from .base_multitemper import (MultiTemperedMetadataIO, MultiTemperedMCMCIO)
from .posterior import PosteriorFile


class PTEmceeFile(MultiTemperedMCMCIO, MultiTemperedMetadataIO,
                  BaseInferenceFile):
    """Class to handle file IO for the ``ptemcee`` sampler."""

    name = 'ptemcee_file'

    # attributes for setting up an ensemble from file
    _ensemble_attrs = ['jumps_proposed', 'jumps_accepted', 'swaps_proposed',
                      'swaps_accepted', 'logP', 'logl']
    @property
    def starting_betas(self):
        """The starting betas that were used."""
        return self[self.sampler_group].attrs["starting_betas"]

    def write_sampler_metadata(self, sampler):
        """Adds writing ptemcee-specific metadata to MultiTemperedMCMCIO.
        """
        super(PTEmceeFile, self).write_sampler_metadata(sampler)
        group = self[self.sampler_group]
        group.attrs["starting_betas"] = sampler.starting_betas
        group.attrs["adaptive"] = sampler.adaptive
        group.attrs["adaptation_lag"] = sampler.adaptation_lag
        group.attrs["adaptation_time"] = sampler.adaptation_time
        group.attrs["scale_factor"] = sampler.scale_factor

    def write_betas(self, betas):
        """Writes the betas to sampler group.

        As the betas may change with iterations, this writes the betas as
        a ntemps x niterations array to the file.
        """
        group = self[self.sampler_group]
        ntemps, niterations = betas.shape
        try:
            fp_ntemps, fp_niterations = group['betas'].shape
            # check that the number of temps matches
            assert ntemps == fp_ntemps, ("length of betas' first axis must "
                                         "match the betas array in the file")
            istart = fp_niterations
            istop = istart + niterations
            # resize the dataset to accomodate the new temps
            group['betas'].resize(istop, axis=1)
        except KeyError:
            # dataset doesn't exist yet
            istart = 0
            istop = niterations
            group.create_dataset('betas', (ntemps, istop),
                                 maxshape=(ntemps, None),
                                 dtype=float, fletcher32=True)
        group['betas'][:, istart:istop] = betas

    def read_betas(self, thin_start=None, thin_interval=None, thin_end=None,
                   iteration=None):
        """Reads betas from the file.

        Parameters
        ----------
        thin_start : int, optional
            Start reading from the given iteration. Default is start from the
            first iteration.
        thin_interval : int, optional
            Ony read every ``thin_interval`` -th sample. Default is 1.
        thin_end : int, optional
            Stop reading at the given iteration. Default is to end at the last
            iteration.
        iteration : int, optional
            Only read betas from the given iteration. If provided, overrides
            the ``thin_(start|interval|end)`` options.

        Returns
        -------
        array
            A ntemps x niterations array of the betas.
        """
        group = self[self.sampler_group]
        if iteration is not None:
            get_index = int(iteration)
        else:
            get_index = self.get_slice(thin_start=thin_start,
                                       thin_end=thin_end,
                                       thin_interval=thin_interval)
        return group['betas'][:, get_index]

    def write_ensemble_attrs(self, ensemble):
        """Writes ensemble attributes necessary to restart from checkpoint.

        Paramters
        ---------
        ensemble : ptemcee.Ensemble
            The ensemble to write attributes for.
        """
        group = self[self.sampler_group]
        for attr in self._ensemble_attrs:
            vals = getattr(ensemble, attr)
            try:
                group[attr][:] = vals
            except KeyError:
                group[attr] = vals

    def read_ensemble_attrs(self):
        """Reads ensemble attributes from the file.

        Returns
        -------
        dict :
            Dictionary of the ensemble attributes.
        """
        group = self[self.sampler_group]
        return {attr: group[attr][:] for attr in self._ensemble_attrs}

    def write_posterior(self, filename, **kwargs):
        """Write posterior only file

        Parameters
        ----------
        filename : str
            Name of output file to store posterior
        """
        f = h5py.File(filename, 'w')

        # Preserve top-level metadata
        for key in self.attrs:
            f.attrs[key] = self.attrs[key]

        f.attrs['filetype'] = PosteriorFile.name
        s = f.create_group('samples')
        fields = self[self.samples_group].keys()

        # Copy and squash fields into one dimensional arrays
        for field_name in fields:
            fvalue = self[self.samples_group][field_name][:]
            thin = fvalue[0,:,self.thin_start:self.thin_end:self.thin_interval]
            s[field_name] = thin.flatten()
