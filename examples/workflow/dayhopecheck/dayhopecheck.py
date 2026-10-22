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

__author__ = "Ian Harry <ian.harry@astro.cf.ac.uk>"
__version__ = pycbc.version.git_verbose_msg
__date__ = pycbc.version.date
__program__ = "dayhopecheck"


import argparse
import logging
import os

import igwn_segments as segments
from glue.segmentdb import segmentdb_utils
from igwn_ligolw import ligolw
from igwn_ligolw import utils as ligolw_utils

import pycbc.workflow as _workflow

logging.basicConfig(
    format="%(asctime)s:%(levelname)s : %(message)s",
    level=logging.INFO,
    datefmt="%I:%M:%S",
)

# command line options
_desc = __doc__[1:]
parser = argparse.ArgumentParser(description=_desc)
parser.add_argument("--version", action="version", version=__version__)
parser.add_argument(
    "-d", "--output-dir", required=True, help="Path to output directory."
)
_workflow.add_workflow_command_line_group(parser)
args = parser.parse_args()

workflow = _workflow.Workflow(args)
curr_dir = os.path.abspath(args.output_dir)
seg_dir = os.path.join(curr_dir, "segments")
df_dir = os.path.join(curr_dir, "datafind")

print("BEGIN BY GENERATING SCIENCE AND CAT_X VETOES")


def segment_report(s_segs):
    full_len = 0
    full_num = 0
    short_len = 0
    short_num = 0
    long_len = 0
    long_num = 0
    for ifo in s_segs.keys():
        for seg in s_segs[ifo]:
            full_len += abs(seg)
            full_num += 1
            if abs(seg) > 500:
                short_len += abs(seg)
                short_num += 1
            if abs(seg) > 2000:
                long_len += abs(seg)
                long_num += 1
        print(
            "For ifo %s there is %d seconds of data in %d segments, %d seconds (%d unique segments) in segments longer than 500s and %d seconds (%d unique segments) longer than 2000s."
            % (ifo, full_len, full_num, short_len, short_num, long_len, long_num)
        )


science_segs, segs_list = _workflow.setup_segment_generation(workflow, seg_dir)

segment_report(science_segs)

print()
print()

print("RUNNING DATAFIND")
datafinds, science_segs = _workflow.setup_datafind_workflow(
    workflow, science_segs, df_dir, segs_list
)

# This is needed to know what times will be analysed by daily ahope
# Template bank stuff
banks = _workflow.setup_tmpltbank_workflow(workflow, science_segs, datafinds, df_dir)
# Do matched-filtering
insps = _workflow.setup_matchedfltr_workflow(
    workflow, science_segs, datafinds, banks, df_dir
)

# Now construct the summary XML file

outdoc = ligolw.Document()
outdoc.appendChild(ligolw.LIGO_LW())
# FIXME: PROGRAM NAME and dictionary of opts should be variables defined up above
proc_id = outdoc.register_process("dayhopetest", vars(args)).process_id
for ifo in workflow.ifos:
    # Lets get the segment lists we need
    seg_ifo_files = segs_list.find_output_with_ifo(ifo)
    # SCIENCE
    sci_seg_file = seg_ifo_files.find_output_with_tag("SCIENCE")
    assert len(sci_seg_file) == 1
    sci_seg_file = sci_seg_file[0]
    sci_segs = sci_seg_file.segmentList
    # SCIENCE_OK
    sci_ok_seg_file = seg_ifo_files.find_output_with_tag("SCIENCE_OK")
    assert len(sci_ok_seg_file) == 1
    sci_ok_seg_file = sci_ok_seg_file[0]
    sci_ok_segs = sci_ok_seg_file.segmentList
    # SCIENCE_AVAILABLE
    sci_available_seg_file = seg_ifo_files.find_output_with_tag("SCIENCE_AVAILABLE")
    assert len(sci_available_seg_file) == 1
    sci_available_seg_file = sci_available_seg_file[0]
    sci_available_segs = sci_available_seg_file.segmentList
    # ANALYSABLE - This one needs to come from inspiral outs
    analysable_segs = insps.get_times_covered_by_files()

    # And add these to the output file
    # Start with the segment summary
    summ_segs = segments.segmentlist([workflow.analysis_time])
    sci_def_id = segmentdb_utils.add_to_segment_definer(
        outdoc, proc_id, ifo, "CBC_DAYHOPE_SCIENCE", 0
    )
    sciok_def_id = segmentdb_utils.add_to_segment_definer(
        outdoc, proc_id, ifo, "CBC_DAYHOPE_SCIENCE_OK", 0
    )
    sciavailable_def_id = segmentdb_utils.add_to_segment_definer(
        outdoc, proc_id, ifo, "CBC_DAYHOPE_SCIENCE_AVAILABLE", 0
    )
    analysable_def_id = segmentdb_utils.add_to_segment_definer(
        outdoc, proc_id, ifo, "CBC_DAYHOPE_ANALYSABLE", 0
    )

    segmentdb_utils.add_to_segment(outdoc, proc_id, sci_def_id, sci_segs)
    segmentdb_utils.add_to_segment(outdoc, proc_id, sciok_def_id, sci_ok_segs)
    segmentdb_utils.add_to_segment(
        outdoc, proc_id, sciavailable_def_id, sci_available_segs
    )
    segmentdb_utils.add_to_segment(outdoc, proc_id, analysable_def_id, analysable_segs)

    segmentdb_utils.add_to_segment_summary(
        outdoc, proc_id, sci_def_id, summ_segs, comment=""
    )
    segmentdb_utils.add_to_segment_summary(
        outdoc, proc_id, sciok_def_id, summ_segs, comment=""
    )
    segmentdb_utils.add_to_segment_summary(
        outdoc, proc_id, sciavailable_def_id, summ_segs, comment=""
    )
    segmentdb_utils.add_to_segment_summary(
        outdoc, proc_id, analysable_def_id, summ_segs, comment=""
    )

ligolw_utils.write_filename(outdoc, "SUMMARY.xml")
