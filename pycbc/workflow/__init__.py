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
import os.path

from pycbc.workflow.configuration import *
try:
    from pycbc.workflow.core import *
    from pycbc.workflow.legacy_ihope import *
    from pycbc.workflow.grb_utils import *
    from pycbc.workflow.jobsetup import *
    from pycbc.workflow.psd import *
    from pycbc.workflow.matched_filter import *
    from pycbc.workflow.datafind import *
    from pycbc.workflow.segment import *
    from pycbc.workflow.tmpltbank import *
    from pycbc.workflow.psdfiles import *
    from pycbc.workflow.splittable import *
    from pycbc.workflow.coincidence import *
    from pycbc.workflow.injection import *
    from pycbc.workflow.postprocessing_cohptf import *
    from pycbc.workflow.plotting import *
    from pycbc.workflow.minifollowups import *
except ImportError:
    pass

# Set the configuration file base directory
INI_FILE_DIRECTORY = os.path.join(os.path.dirname(__file__), 'ini_files')

# Set the pycbc workflow specific pegasus configuration and planning files
PEGASUS_FILE_DIRECTORY = os.path.join(os.path.dirname(__file__), 'pegasus_files')
