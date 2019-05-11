# Copyright (C) 2019 Collin Capano
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

"""Provides abstract base class for all samplers."""

from __future__ import absolute_import

from abc import (ABCMeta, abstractmethod)

from .base_hdf import BaseInferenceFile


class BaseSamplerFile(BaseInferenceFile):
    """Base HDF class for all samplers.

    This adds abstract methods ``write_resume_point`` and
    ``write_sampler_metadata`` to :py:class:`BaseInferenceFile`.
    """
    __metaclass__ = ABCMeta

    @abstractmethod
    def write_resume_point(self):
        """Should write the point that a sampler starts up.

        How the resume point is indexed is up to the sampler. For example,
        MCMC samplers use the number of iterations that are stored in the
        checkpoint file.
        """
        pass

    @abstractmethod
    def write_sampler_metadata(self, sampler):
        """This should write the given sampler's metadata to the file.

        This should also include the model's metadata.
        """
        pass
