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


import numpy

from .base_sampler import BaseSamplerFile
from .base_mcmc import EnsembleMCMCMetadataIO
from .base_multitemper import (CommonMultiTemperedMetadataIO,
                               write_samples,
                               ensemble_read_raw_samples)


class EmceePTFile(EnsembleMCMCMetadataIO, CommonMultiTemperedMetadataIO,
                  BaseSamplerFile):
    """Class to handle file IO for the ``emcee`` sampler."""

    name = 'emcee_pt_file'

    @property
    def betas(self):
        """The betas that were used."""
        return self[self.sampler_group].attrs["betas"]

    def write_samples(self, samples, **kwargs):
        r"""Writes samples to the given file.

        Calls :py:func:`base_multitemper.write_samples`. See that function for
        details.

        Parameters
        ----------
        samples : dict
            The samples to write. Each array in the dictionary should have
            shape ntemps x nwalkers x niterations.
        \**kwargs :
            All other keyword arguments are passed to
            :py:func:`base_multitemper.write_samples`.
        """
        write_samples(self, samples, **kwargs)

    def read_raw_samples(self, fields, **kwargs):
        r"""Base function for reading samples.

        Calls :py:func:`base_multitemper.ensemble_read_raw_samples`. See that
        function for details.

        Parameters
        -----------
        fields : list
            The list of field names to retrieve.
        \**kwargs :
            All other keyword arguments are passed to
            :py:func:`base_multitemper.ensemble_read_raw_samples`.

        Returns
        -------
        dict
            A dictionary of field name -> numpy array pairs.
        """
        return ensemble_read_raw_samples(self, fields, **kwargs)

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
