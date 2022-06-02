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


import numpy
from pickle import UnpicklingError
from epsie import load_state

from .base_sampler import BaseSamplerFile
from .base_mcmc import MCMCMetadataIO
from .base_multitemper import (CommonMultiTemperedMetadataIO,
                               write_samples,
                               read_raw_samples)


class EpsieFile(MCMCMetadataIO, CommonMultiTemperedMetadataIO,
                BaseSamplerFile):
    """Class to handle IO for Epsie's parallel-tempered sampler."""

    name = 'epsie_file'

    @property
    def nchains(self):
        """Alias for nwalkers."""
        return self.nwalkers

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

    @property
    def seed(self):
        """The sampler's seed."""
        # convert seed from str back to int (see setter below for reason)
        return int(self[self.sampler_group].attrs['seed'])

    @seed.setter
    def seed(self, seed):
        """Store the sampler's seed."""
        # epsie uses the numpy's new random generators, which use long integers
        # for seeds. hdf5 doesn't know how to handle long integers, so we'll
        # store it as a string
        self[self.sampler_group].attrs['seed'] = str(seed)

    def write_sampler_metadata(self, sampler):
        """Adds writing seed and betas to MultiTemperedMCMCIO.
        """
        super(EpsieFile, self).write_sampler_metadata(sampler)
        self.seed = sampler.seed
        self.write_data("betas", sampler.betas, path=self.sampler_group)

    def thin(self, thin_interval):
        """Thins the samples on disk to the given thinning interval.

        Also thins the acceptance ratio and the temperature data, both of
        which are stored in the ``sampler_info`` group.
        """
        # We'll need to know what the new interval to thin by will be
        # so we can properly thin the acceptance ratio and temperatures swaps.
        # We need to do this before calling the base thin, as we need to know
        # what the current thinned by is.
        new_interval = thin_interval // self.thinned_by
        # now thin the samples
        super(EpsieFile, self).thin(thin_interval)
        # thin the acceptance ratio
        self._thin_data(self.sampler_group, ['acceptance_ratio'],
                        new_interval)
        # thin the temperature swaps; since these may not happen every
        # iteration, the thin interval we use for these is different
        ts_group = '/'.join([self.sampler_group, 'temperature_swaps'])
        ts_thin_interval = new_interval // self.swap_interval
        if ts_thin_interval > 1:
            self._thin_data(ts_group, ['swap_index'],
                            ts_thin_interval)
            self._thin_data(ts_group, ['acceptance_ratio'],
                            ts_thin_interval)

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

        Calls :py:func:`base_multitemper.read_raw_samples`. See that
        function for details.

        Parameters
        -----------
        fields : list
            The list of field names to retrieve.
        \**kwargs :
            All other keyword arguments are passed to
            :py:func:`base_multitemper.read_raw_samples`.

        Returns
        -------
        dict
            A dictionary of field name -> numpy array pairs.
        """
        return read_raw_samples(self, fields, **kwargs)

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

    def read_acceptance_ratio(self, temps=None, chains=None):
        """Reads the acceptance ratios.

        Ratios larger than 1 are set back to 1 before returning.

        Parameters
        -----------
        temps : (list of) int, optional
            The temperature index (or a list of indices) to retrieve. If None,
            acceptance ratios from all temperatures and all chains will be
            retrieved.
        chains : (list of) int, optional
            The chain index (or a list of indices) to retrieve. If None,
            ratios from all chains will be obtained.

        Returns
        -------
        array
            Array of acceptance ratios with shape (requested temps,
            requested chains, niterations).
        """
        group = self.sampler_group + '/acceptance_ratio'
        if chains is None:
            wmask = numpy.ones(self.nchains, dtype=bool)
        else:
            wmask = numpy.zeros(self.nchains, dtype=bool)
            wmask[chains] = True
        if temps is None:
            tmask = numpy.ones(self.ntemps, dtype=bool)
        else:
            tmask = numpy.zeros(self.ntemps, dtype=bool)
            tmask[temps] = True
        all_ratios = self[group][:]
        # make sure values > 1 are set back to 1
        all_ratios[all_ratios > 1] = 1.
        return all_ratios[numpy.ix_(tmask, wmask)]

    def read_acceptance_rate(self, temps=None, chains=None):
        """Reads the acceptance rate.

        This calls :py:func:`read_acceptance_ratio`, then averages the ratios
        over all iterations to get the average rate.

        Parameters
        -----------
        temps : (list of) int, optional
            The temperature index (or a list of indices) to retrieve. If None,
            acceptance rates from all temperatures and all chains will be
            retrieved.
        chains : (list of) int, optional
            The chain index (or a list of indices) to retrieve. If None,
            rates from all chains will be obtained.

        Returns
        -------
        array
            Array of acceptance ratios with shape (requested temps,
            requested chains).
        """
        all_ratios = self.read_acceptance_ratio(temps, chains)
        # average over the number of iterations
        all_ratios = all_ratios.mean(axis=-1)
        return all_ratios

    def read_acceptance_fraction(self, temps=None, walkers=None):
        """Alias for :py:func:`read_acceptance_rate`.
        """
        return self.read_acceptance_rate(temps=temps, chains=walkers)

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

    @staticmethod
    def _get_optional_args(args, opts, err_on_missing=False, **kwargs):
        # need this to make sure options called "walkers" are renamed to
        # "chains"
        parsed = BaseSamplerFile._get_optional_args(
            args, opts, err_on_missing=err_on_missing, **kwargs)
        try:
            chains = parsed.pop('walkers')
            parsed['chains'] = chains
        except KeyError:
            pass
        return parsed
