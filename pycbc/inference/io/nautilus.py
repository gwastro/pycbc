# Copyright (C) 2026  Alex Nitz
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
"""Provides IO for the nautilus sampler"""
from .dynesty import DynestyFile


class NautilusFile(DynestyFile):
    """Class to handle file IO for the ``nautilus`` sampler.

    Nautilus returns weighted samples in the same form as dynesty (samples
    plus log-weights), so turning them into a posterior is identical.
    """

    name = 'nautilus_file'
    # nautilus' own state is kept verbatim under this path
    state_path = 'sampler_info/nautilus'

    def validate(self):
        """A nautilus checkpoint is only resumable if nautilus' own state is
        in it, since that is what nautilus reads to pick a run back up.
        """
        return self.state_path in self
