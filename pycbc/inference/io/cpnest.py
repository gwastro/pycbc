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
from .posterior import PosteriorFile


class CPNestFile(PosteriorFile):
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
        self[self.sampler_group].attrs['nwalkers'] = sampler.nlive
        # write the model's metadata
        sampler.model.write_metadata(self)
