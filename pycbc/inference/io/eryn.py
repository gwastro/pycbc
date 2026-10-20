# Copyright (C) 2026 Alexander Nitz
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


"""Provides I/O support for the eryn sampler.

The on-disk layout mirrors that of the ``ptemcee`` sampler: samples and model
stats are stored as ``ntemps x nwalkers x niterations`` arrays, and the
(possibly adapting) inverse temperatures are stored as an
``ntemps x niterations`` ``betas`` dataset in the sampler-info group. This lets
all of the standard multi-tempered reading/plotting tools work unchanged.
"""


from . import base_mcmc
from .base_mcmc import EnsembleMCMCMetadataIO
from .base_multitemper import (
    CommonMultiTemperedMetadataIO,
    ensemble_read_raw_samples,
    write_samples,
)
from .base_sampler import BaseSamplerFile


class ErynFile(EnsembleMCMCMetadataIO, CommonMultiTemperedMetadataIO,
               BaseSamplerFile):
    """Class to handle file IO for the ``eryn`` sampler."""

    name = 'eryn_file'

    def write_sampler_metadata(self, sampler):
        """Adds writing eryn-specific metadata to MultiTemperedMCMCIO."""
        super(ErynFile, self).write_sampler_metadata(sampler)
        group = self[self.sampler_group]
        group.attrs["starting_betas"] = sampler.starting_betas
        group.attrs["adaptive"] = sampler.adaptive
        # record the move schedule that was used
        group.attrs["moves"] = sampler.moves_labels

    @property
    def starting_betas(self):
        """The starting betas that were used."""
        return self[self.sampler_group].attrs["starting_betas"]

    def write_betas(self, betas, last_iteration=None):
        """Writes the betas to the sampler group.

        As the betas may change with iterations (adaptive tempering), this
        writes the betas as a ``ntemps x niterations`` array to the file.
        """
        # reuse the single-temperature write_samples for the thinning settings
        base_mcmc.write_samples(self, {'betas': betas},
                                last_iteration=last_iteration,
                                samples_group=self.sampler_group)

    def read_betas(self, thin_start=None, thin_interval=None, thin_end=None,
                   iteration=None):
        """Reads betas from the file.

        Parameters
        -----------
        thin_start : int, optional
            Start reading from the given iteration. Default is to start from
            the first iteration.
        thin_interval : int, optional
            Only read every ``thin_interval`` -th sample. Default is 1.
        thin_end : int, optional
            Stop reading at the given iteration. Default is to end at the last
            iteration.
        iteration : int, optional
            Only read the given iteration. If this provided, it overrides
            the ``thin_(start|interval|end)`` options.

        Returns
        -------
        array
            A ntemps x niterations array of the betas.
        """
        slc = base_mcmc._ensemble_get_index(self, thin_start=thin_start,
                                            thin_interval=thin_interval,
                                            thin_end=thin_end,
                                            iteration=iteration)
        betas = self[self.sampler_group]['betas'][:]
        return betas[:, slc]

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
        ----------
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
