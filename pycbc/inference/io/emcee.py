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


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""Provides IO for the emcee sampler.
"""
import h5py, numpy
from .base_hdf import BaseInferenceFile
from .base_mcmc import (MCMCMetadataIO, SingleTempMCMCIO)
from .posterior import PosteriorFile

class EmceeFile(SingleTempMCMCIO, MCMCMetadataIO, BaseInferenceFile):
    """Class to handle file IO for the ``emcee`` sampler."""

    name = 'emcee_file'

    def read_acceptance_fraction(self, walkers=None):
        """Reads the acceptance fraction.

        Parameters
        -----------
        walkers : (list of) int, optional
            The walker index (or a list of indices) to retrieve. If None,
            samples from all walkers will be obtained.

        Returns
        -------
        array
            Array of acceptance fractions with shape (requested walkers,).
        """
        group = self.sampler_group + '/acceptance_fraction'
        if walkers is None:
            wmask = numpy.ones(self.nwalkers, dtype=bool)
        else:
            wmask = numpy.zeros(self.nwalkers, dtype=bool)
            wmask[walkers] = True
        return self[group][wmask]

    def write_acceptance_fraction(self, acceptance_fraction):
        """Write acceptance_fraction data to file. Results are written to
        the ``[sampler_group]/acceptance_fraction``.

        Parameters
        -----------
        acceptance_fraction : numpy.ndarray
            Array of acceptance fractions to write.
        """
        group = self.sampler_group + '/acceptance_fraction'
        try:
            self[group][:] = acceptance_fraction
        except KeyError:
            # dataset doesn't exist yet, create it
            self[group] = acceptance_fraction

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
            thin = fvalue[:,self.thin_start:self.thin_end:self.thin_interval]
            s[field_name] = thin.flatten()
