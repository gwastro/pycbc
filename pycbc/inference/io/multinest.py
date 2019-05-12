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

from .base_sampler import BaseSamplerFile


class MultinestFile(BaseSamplerFile):
    """Class to handle file IO for the ``multinest`` sampler."""

    name = 'multinest_file'

    def write_samples(self, samples, parameters=None):
        """Writes samples to the given file.

        Results are written to ``samples_group/{vararg}``, where ``{vararg}``
        is the name of a model params. The samples are written as an
        array of length ``niterations``.

        Parameters
        ----------
        samples : dict
            The samples to write. Each array in the dictionary should have
            length niterations.
        parameters : list, optional
            Only write the specified parameters to the file. If None, will
            write all of the keys in the ``samples`` dict.
        """
        niterations = len(samples.values()[0])
        assert all(len(p) == niterations for p in samples.values()), (
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
        importance_lnz : float, optional
            The importance-weighted log of the evidence.
        importance_dlnz : float, optional
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
        self[self.sampler_group].attrs['nlivepoints'] = sampler.nlivepoints
        # write the model's metadata
        sampler.model.write_metadata(self)
