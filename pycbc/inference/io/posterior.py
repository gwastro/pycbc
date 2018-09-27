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

from .base_hdf import BaseInferenceFile
from .base_mcmc import MCMCIO


class PosteriorFile(BaseInferenceFile):
    """Class to handle file IO for the ``emcee`` sampler."""

    name = 'posterior_file'

    def read_raw_samples(self, fields, **kwargs):
        return {field: field[:] for field in self[self.samples_group]}        

    def write_posterior(self, filename, **kwargs):
        """Write me."""
        raise RuntimeError("This function should not be needed. "
                           "Why was this called????")
