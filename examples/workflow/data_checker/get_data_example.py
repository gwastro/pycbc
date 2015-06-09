#!/usr/bin/env python

# Copyright (C) 2013 Ian W. Harry
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

"""
An example of how to query science + CAT_X segments and do datafind with pycbc.workflow.
"""

import pycbc
import pycbc.version
__author__  = "Ian Harry <ian.harry@astro.cf.ac.uk>"
__version__ = pycbc.version.git_verbose_msg
__date__    = pycbc.version.date
__program__ = "get_data_example"


import os, sys
import argparse
import logging
from glue import segments
import pycbc.workflow as _workflow

logging.basicConfig(format='%(asctime)s:%(levelname)s : %(message)s', 
                    level=logging.INFO, datefmt='%I:%M:%S')

# command line options
_desc = __doc__[1:]
parser = argparse.ArgumentParser(description=_desc)
parser.add_argument('--version', action='version', version=__version__)
_workflow.add_workflow_command_line_group(parser)
args = parser.parse_args()

workflow = _workflow.Workflow(args)
currDir = os.getcwd()
segDir = os.path.join(currDir,"segments")
if not os.path.exists(segDir+'/logs'):
    os.makedirs(segDir+'/logs')
dfDir = os.path.join(currDir,"datafind")
if not os.path.exists(dfDir+'/logs'):
    os.makedirs(dfDir+'/logs')

scienceSegs, segsList = _workflow.setup_segment_generation(workflow, segDir)
datafinds, scienceSegs = _workflow.setup_datafind_workflow(workflow, scienceSegs, 
                           dfDir, segsList)

# scienceSegs is a glue.segmentlist of the times you should analyse.
# It contains science times, that are present on disk, with CAT_1 times
# removed.
