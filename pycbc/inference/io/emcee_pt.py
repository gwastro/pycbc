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


"""Provides I/O support for emcee_pt.
"""

from __future__ import absolute_import
import h5py, numpy
from .base_sampler import BaseSamplerFile
from .base_multitemper import (MultiTemperedMetadataIO, MultiTemperedMCMCIO)
from .posterior import PosteriorFile


class EmceePTFile(MultiTemperedMCMCIO, MultiTemperedMetadataIO,
                  BaseSamplerFile):
    """Class to handle file IO for the ``emcee`` sampler."""

    name = 'emcee_pt_file'

    @property
    def betas(self):
        """The betas that were used."""
        return self[self.sampler_group].attrs["betas"]

    def write_sampler_metadata(self, sampler):
        """Adds writing betas to MultiTemperedMCMCIO.
        """
        super(EmceePTFile, self).write_sampler_metadata(sampler)
        self[self.sampler_group].attrs["betas"] = sampler.betas

    def read_acceptance_fraction(self, temps=None, walkers=None):
        """Reads the acceptance fraction.

        Parameters
        -----------
        temps : (list of) int, optional
            The temperature index (or a list of indices) to retrieve. If None,
            acfs from all temperatures and all walkers will be retrieved.
        walkers : (list of) int, optional
            The walker index (or a list of indices) to retrieve. If None,
            samples from all walkers will be obtained.

        Returns
        -------
        array
            Array of acceptance fractions with shape (requested temps,
            requested walkers).
        """
        group = self.sampler_group + '/acceptance_fraction'
        if walkers is None:
            wmask = numpy.ones(self.nwalkers, dtype=bool)
        else:
            wmask = numpy.zeros(self.nwalkers, dtype=bool)
            wmask[walkers] = True
        if temps is None:
            tmask = numpy.ones(self.ntemps, dtype=bool)
        else:
            tmask = numpy.zeros(self.ntemps, dtype=bool)
            tmask[temps] = True
        return self[group][:][numpy.ix_(tmask, wmask)]

    def write_acceptance_fraction(self, acceptance_fraction):
        """Write acceptance_fraction data to file.

        Results are written to ``[sampler_group]/acceptance_fraction``; the
        resulting dataset has shape (ntemps, nwalkers).

        Parameters
        -----------
        acceptance_fraction : numpy.ndarray
            Array of acceptance fractions to write. Must have shape
            ntemps x nwalkers.
        """
        # check
        assert acceptance_fraction.shape == (self.ntemps, self.nwalkers), (
            "acceptance fraction must have shape ntemps x nwalker")
        group = self.sampler_group + '/acceptance_fraction'
        try:
            self[group][:] = acceptance_fraction
        except KeyError:
            # dataset doesn't exist yet, create it
            self[group] = acceptance_fraction

    def read_posterior_samples(self, parameters):
        """Read posterior samples.

        These are taken from the coldest temperature, starting from the
        burn in iteration, and thinned by the ACL.

        Parameters
        ----------
        parameters : list of str
            The names of the parameters to read.

        Returns
        -------
        FieldArray :
            The posterior samples, as a 1D ``FieldArray``.
        """
        # the burn in and ACL should be stored by the thin_start and
        # thin_interval attributes, respectively
        return self.read_samples(parameters, temps=0, 
                                 thin_start=self.thin_start,
                                 thin_interval=self.thin_interval,
                                 thin_end=self.thin_end,
                                 flatten=True)
