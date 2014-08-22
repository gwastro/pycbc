# Copyright (C) 2013 Ian Harry
#
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

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""
This package provides the utilities to construct an inspiral workflow for
performing a coincident CBC matched-filter analysis on gravitational-wave
interferometer data
"""

from configuration import *
from core import *
from legacy_ihope import *
from jobsetup import *
from matched_filter import *
from datafind import *
from segment import *
from tmpltbank import *
from splittable import *
from coincidence import *
from injection import *
from timeslides import *
from postprocessing_prep import *
from postprocessing import *
from analysislogging import *
