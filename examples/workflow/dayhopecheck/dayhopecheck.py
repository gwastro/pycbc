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
Program for determining the following segment lists and dumping these to a
single summary XML file:

1) Science time, according to the segment database (X1:CBC_DAYHOPE_SCIENCE)
2) Science time - CAT_1, using a supplied veto-definer file to define CAT_1.
   If you define a minimum segment length, this is also removed here.
   (X1:CBC_DAYHOPE_SCIENCE_OK)
3) SCIENCE_OK - times not available for datafind. Depending on how you set the
   .ini file options this can include SCIENCE_OK times that are not found on
   the local LDR server, and files that the local LDR server returns but are
   not actually visible (using os.path.isfile)
   (X1:CBC_DAYHOPE_SCIENCE_AVAILABLE)
4) Time analysable by daily [i,a]hope. This will compute the times that the
   daily CBC analysis can analyse, which must be a subset of SCIENCE_AVAILABLE.
   The options given in the config file tell the code what data can be 
   analysed by the daily analysis. (X1:CBC_DAYHOPE_ANALYSABLE)
"""

import pycbc
import pycbc.version
__author__  = "Ian Harry <ian.harry@astro.cf.ac.uk>"
__version__ = pycbc.version.git_verbose_msg
__date__    = pycbc.version.date
__program__ = "dayhopecheck"


import os
import copy
import logging
import argparse
from glue import segments
import pycbc.workflow as _workflow

from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import utils as ligolw_utils
from glue.ligolw.utils import segments as ligolw_segments
from glue.ligolw.utils import process as ligolw_process

from glue.segmentdb import segmentdb_utils

logging.basicConfig(format='%(asctime)s:%(levelname)s : %(message)s', \
                    level=logging.INFO,datefmt='%I:%M:%S')

# command line options
_desc = __doc__[1:]
parser = argparse.ArgumentParser(description=_desc)
parser.add_argument('--version', action='version', version=__version__)
parser.add_argument("-d", "--output-dir", required=True,\
                    help="Path to output directory.")
_workflow.add_workflow_command_line_group(parser)
args = parser.parse_args()

workflow = _workflow.Workflow(args)
currDir = os.path.abspath(args.output_dir)
segDir = os.path.join(currDir,"segments")
dfDir = os.path.join(currDir,"datafind")

print "BEGIN BY GENERATING SCIENCE AND CAT_X VETOES"

def segment_report(sSegs):
    fullLen = 0
    fullNum = 0
    shortLen = 0
    shortNum = 0
    longLen = 0
    longNum = 0
    for ifo in sSegs.keys():
        for seg in sSegs[ifo]:
            fullLen += abs(seg)
            fullNum += 1
            if abs(seg) > 500:
                shortLen+=abs(seg)
                shortNum+=1
            if abs(seg) > 2000:
                longLen+=abs(seg)
                longNum+=1
        print "For ifo %s there is %d seconds of data in %d segments, %d seconds (%d unique segments) in segments longer than 500s and %d seconds (%d unique segments) longer than 2000s." %(ifo, fullLen, fullNum, shortLen, shortNum, longLen, longNum)


scienceSegs, segsList = _workflow.setup_segment_generation(workflow, segDir)

segment_report(scienceSegs)

print
print

print "RUNNING DATAFIND"
datafinds, scienceSegs = _workflow.setup_datafind_workflow(workflow, scienceSegs,
                     dfDir, segsList)

# This is needed to know what times will be analysed by daily ahope
# Template bank stuff
banks = _workflow.setup_tmpltbank_workflow(workflow, scienceSegs, datafinds,
                                       dfDir)
# Do matched-filtering
insps = _workflow.setup_matchedfltr_workflow(workflow, scienceSegs, datafinds,
                                         banks, dfDir)

# Now construct the summary XML file

outdoc = ligolw.Document()
outdoc.appendChild(ligolw.LIGO_LW())
# FIXME: PROGRAM NAME and dictionary of opts should be variables defined up above
proc_id = ligolw_process.register_to_xmldoc(outdoc, 'dayhopetest',
                                            vars(args) ).process_id
for ifo in workflow.ifos:
    # Lets get the segment lists we need
    segIfoFiles = segsList.find_output_with_ifo(ifo)
    # SCIENCE
    sciSegFile = segIfoFiles.find_output_with_tag('SCIENCE')
    assert(len(sciSegFile) == 1)
    sciSegFile = sciSegFile[0]
    sciSegs = sciSegFile.segmentList
    # SCIENCE_OK
    sciokSegFile = segIfoFiles.find_output_with_tag('SCIENCE_OK')
    assert(len(sciokSegFile) == 1)
    sciokSegFile = sciokSegFile[0]
    sciokSegs = sciokSegFile.segmentList
    # SCIENCE_AVAILABLE
    sciavailableSegFile = segIfoFiles.find_output_with_tag('SCIENCE_AVAILABLE')
    assert(len(sciavailableSegFile) == 1)
    sciavailableSegFile = sciavailableSegFile[0]
    sciavailableSegs = sciavailableSegFile.segmentList
    # ANALYSABLE - This one needs to come from inspiral outs
    analysableSegs = insps.get_times_covered_by_files()
   
    # And add these to the output file
    # Start with the segment summary
    summSegs = segments.segmentlist([workflow.analysis_time])
    sci_def_id = segmentdb_utils.add_to_segment_definer(outdoc, proc_id, ifo,
                                                      "CBC_DAYHOPE_SCIENCE", 0)
    sciok_def_id = segmentdb_utils.add_to_segment_definer(outdoc, proc_id, ifo,
                                                   "CBC_DAYHOPE_SCIENCE_OK", 0)
    sciavailable_def_id = segmentdb_utils.add_to_segment_definer(outdoc,
                              proc_id, ifo, "CBC_DAYHOPE_SCIENCE_AVAILABLE", 0)
    analysable_def_id = segmentdb_utils.add_to_segment_definer(outdoc, proc_id,
                                              ifo, "CBC_DAYHOPE_ANALYSABLE", 0)
    
    segmentdb_utils.add_to_segment(outdoc, proc_id, sci_def_id, sciSegs)
    segmentdb_utils.add_to_segment(outdoc, proc_id, sciok_def_id, sciokSegs)
    segmentdb_utils.add_to_segment(outdoc, proc_id, sciavailable_def_id,
                                   sciavailableSegs)
    segmentdb_utils.add_to_segment(outdoc, proc_id, analysable_def_id,
                                   analysableSegs)

    segmentdb_utils.add_to_segment_summary(outdoc, proc_id, sci_def_id,
                                           summSegs, comment='')
    segmentdb_utils.add_to_segment_summary(outdoc, proc_id, sciok_def_id,
                                           summSegs, comment='')
    segmentdb_utils.add_to_segment_summary(outdoc, proc_id, sciavailable_def_id,
                                           summSegs, comment='')
    segmentdb_utils.add_to_segment_summary(outdoc, proc_id, analysable_def_id,
                                           summSegs, comment='')

ligolw_utils.write_filename(outdoc, "SUMMARY.xml")

