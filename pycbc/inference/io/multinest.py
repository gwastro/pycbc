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


"""Provides I/O support for multinest.
"""

from __future__ import absolute_import
import h5py, numpy
from .base_hdf import BaseInferenceFile
from .base_multitemper import (MultiTemperedMetadataIO, MultiTemperedMCMCIO)
from .posterior import PosteriorFile


class MultinestFile(BaseInferenceFile):
    """Class to handle file IO for the ``multinest`` sampler."""

    name = 'multinest_file'

    def write_samples(self, samples, parameters=None):
        niterations = len(samples.values()[0])
        assert all(len(p) == niterations
                   for p in samples.values()), (
               "all samples must have the same shape")
        group = self.samples_group + '/{name}'
        if parameters is None:
            parameters = samples.keys()
        # loop over number of dimensions
        for param in parameters:
            dataset_name = group.format(name=param)
            try:
                fp_niterations = len(self[dataset_name])
                if niterations != fp_niterations:
                    # resize the dataset
                    self[dataset_name].resize(niterations, axis=0)
            except KeyError:
                # dataset doesn't exist yet
                self.create_dataset(dataset_name, (niterations,),
                                    maxshape=(None,),
                                    dtype=samples[param].dtype,
                                    fletcher32=True)
            self[dataset_name][:] = samples[param]

    def write_stats(self, stats):
        group = 'likelihood_stats' + '/{name}'
        for k, v in stats.items():
            dataset_name = group.format(name=k)
            try:
                self[dataset_name] = v
            except KeyError:
                self.create_dataset(dataset_name, (1,))

    def write_logevidence(self, lnz, dlnz, importance_lnz, importance_dlnz):
        """Writes the given log evidence and its error.

        Results are saved to file's 'log_evidence' and 'dlog_evidence'
        attributes, as well as the importance-weighted versions of these
        stats if they exist.

        Parameters
        ----------
        lnz : float
            The log of the evidence.
        dlnz : float
            The error in the estimate of the log evidence.
        importance_lnz : float
            The importance-weighted log of the evidence.
        importance_dlnz : float
            The error in the importance-weighted estimate of the log evidence.
        """
        self.attrs['log_evidence'] = lnz
        self.attrs['dlog_evidence'] = dlnz
        if all([e is not None for e in [importance_lnz, importance_dlnz]]):
            self.attrs['importance_log_evidence'] = importance_lnz
            self.attrs['importance_dlog_evidence'] = importance_dlnz

    def read_raw_samples(self, fields, iteration=None):
        if isinstance(fields, (str, unicode)):
            fields = [fields]
        # load
        group = self.samples_group + '/{name}'
        arrays = {}
        for name in fields:
            if iteration is not None:
                arr = self[group.format(name=name)][int(iteration)]
            else:
                arr = self[group.format(name=name)][:]
            arrays[name] = arr
        return arrays

    def write_resume_point(self):
        """Keeps a list of the number of iterations that were in a file when a
        run was resumed from a checkpoint."""
        try:
            resume_pts = self.attrs["resume_points"].tolist()
        except KeyError:
            resume_pts = []
        try:
            niterations = self.niterations
        except KeyError:
            niterations = 0
        resume_pts.append(niterations)
        self.attrs["resume_points"] = resume_pts

    @property
    def niterations(self):
        """Returns the number of iterations the sampler was run for."""
        return self[self.sampler_group].attrs['niterations']

    def write_niterations(self, niterations):
        """Writes the given number of iterations to the sampler group."""
        self[self.sampler_group].attrs['niterations'] = niterations

    def write_sampler_metadata(self, sampler):
        """Writes the sampler's metadata."""
        self.attrs['sampler'] = sampler.name
        if self.sampler_group not in self.keys():
            # create the sampler group
            self.create_group(self.sampler_group)
        self[self.sampler_group].attrs['nwalkers'] = sampler.nlivepoints
        # write the model's metadata
        sampler.model.write_metadata(self)

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
        #assert acceptance_fraction.shape == (self.ntemps, self.nwalkers), (
        #    "acceptance fraction must have shape ntemps x nwalker")
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
            thin = fvalue[0,:,self.thin_start:self.thin_end:self.thin_interval]
            s[field_name] = thin.flatten()
