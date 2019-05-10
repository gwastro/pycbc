# Copyright (C) 2018 Alex Nitz
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
"""Provides simplified standard format just for posterior data
"""

from .base_hdf import BaseInferenceFile


class PosteriorFile(BaseInferenceFile):
    """Class to handle file IO for the simplified Posterior file"""

    name = 'posterior_file'

    def read_raw_samples(self, fields, **kwargs):
        samples = self[self.samples_group]
        return {field: samples[field][:] for field in fields}

    def write_posterior(self, filename, **kwargs):
        """Write me."""
        raise NotImplementedError

    def write_resume_point(self):
        raise NotImplementedError

    def write_sampler_metadata(self, sampler):
        raise NotImplementedError

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
