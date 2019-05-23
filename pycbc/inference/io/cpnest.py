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
from .base_sampler import BaseSamplerFile
from .posterior import PosteriorFile


class CPNestFile(BaseSamplerFile):
    """Class to handle file IO for the ``cpnest`` sampler."""

    name = 'cpnest_file'

    def write_resume_point(self):
        pass

    def write_niterations(self, niterations):
        """
        Writes the given number of iterations to the sampler group.
        """
        self[self.sampler_group].attrs['niterations'] = niterations

    def write_sampler_metadata(self, sampler):
        """
        Adds writing betas to MultiTemperedMCMCIO.
        """
        self.attrs['sampler'] = sampler.name
        if self.sampler_group not in self.keys():
            # create the sampler group
            self.create_group(self.sampler_group)
        self[self.sampler_group].attrs['nlivepoints'] = sampler.nlive
        # write the model's metadata
        sampler.model.write_metadata(self)

    def write_samples(self, samples, parameters=None):
        """Writes samples to the given file.

        Results are written to ``samples_group/{vararg}``, where ``{vararg}``
        is the name of a model params. The samples are written as an
        array of length ``niterations``.

        Parameters
        -----------
        samples : dict
            The samples to write. Each array in the dictionary should have
            length niterations.
        parameters : list, optional
            Only write the specified parameters to the file. If None, will
            write all of the keys in the ``samples`` dict.
        """
        # since we're just writing a posterior use
        # PosteriorFile's write_samples
        PosteriorFile.write_samples(self, samples, parameters=parameters)
