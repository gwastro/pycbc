# Copyright (C) 2014  Ian Harry
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
This module is responsible for creating any files that log what the workflow
has done. This includes command line options sent to the workflow, analysis
times used within the workflow, etc.
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/NOTYETCREATED.html
"""

import logging
from pycbc.workflow.core import File, FileList, make_analysis_dir
from glue.ligolw import ligolw, utils
from glue.ligolw.utils import process
from glue.segmentdb import segmentdb_utils
from glue.segments import segmentlist

def setup_analysislogging(workflow, segs_list, insps, args, output_dir,
                          program_name="workflow", tags=[]):
    """
    This module sets up the analysis logging xml file that contains the
    following information:

    * Command line arguments that the code was run with
    * Segment list of times marked as SCIENCE
    * Segment list of times marked as SCIENCE and "OK" ie. not CAT_1 vetoed
    * Segment list of times marked as SCIENCE_OK and present on the cluster
    * The times that will be analysed by the matched-filter jobs

    Parameters
    -----------
    workflow : pycbc.workflow.core.Workflow
        The Workflow instance.
    segs_list : pycbc.workflow.core.FileList
        A list of Files containing the information needed to generate the
        segments above. For segments generated at run time the associated
        segmentlist is a property of this object.
    insps : pycbc.workflow.core.FileList
        The output files from the matched-filtering module. Used to identify
        what times have been analysed in this workflow.
    output_dir : path
        Directory to output any files to.
    program_name : string (optional, default = "workflow")
        The program name to stick in the process/process_params tables.
    tags : list (optional, default = [])
        If given restrict to considering inspiral and segment files that
        are tagged with all tags in this list.
    """
    logging.info("Entering analysis logging module.")
    make_analysis_dir(output_dir)

    # Construct the summary XML file
    outdoc = ligolw.Document()
    outdoc.appendChild(ligolw.LIGO_LW())

    # Add process and process_params tables
    proc_id = process.register_to_xmldoc(outdoc, program_name,
                                            vars(args) ).process_id

    # Now add the various segment lists to this file
    summ_segs = segmentlist([workflow.analysis_time])
    
    # If tags is given filter by tags
    if tags:
        for tag in tags:
            segs_list = segs_list.find_output_with_tag(tag)
            insps = insps.find_output_with_tag(tag)

    for ifo in workflow.ifos:
        # Lets get the segment lists we need
        seg_ifo_files = segs_list.find_output_with_ifo(ifo)
        # SCIENCE
        sci_seg_file = seg_ifo_files.find_output_with_tag('SCIENCE')
        if len(sci_seg_file) == 1:
            sci_seg_file = sci_seg_file[0]
            sci_segs = sci_seg_file.segmentList
            sci_def_id = segmentdb_utils.add_to_segment_definer(outdoc, proc_id,
                                                   ifo, "CBC_WORKFLOW_SCIENCE", 0)
            segmentdb_utils.add_to_segment(outdoc, proc_id, sci_def_id,
                                                                      sci_segs)
            segmentdb_utils.add_to_segment_summary(outdoc, proc_id, sci_def_id,
                                                         summ_segs, comment='')
        elif sci_seg_file:
            # FIXME: While the segment module is still fractured (#127) this
            #        may not work. Please update when #127 is resolved
            pass
            #err_msg = "Got %d files matching %s and %s. Expected 1 or 0." \
            #          %(len(sci_seg_file), ifo, 'SCIENCE')
            #raise ValueError(err_msg)

        # SCIENCE_OK
        sci_ok_seg_file = seg_ifo_files.find_output_with_tag('SCIENCE_OK')
        if len(sci_ok_seg_file) == 1:
            sci_ok_seg_file = sci_ok_seg_file[0]
            sci_ok_segs = sci_ok_seg_file.segmentList
            sci_ok_def_id = segmentdb_utils.add_to_segment_definer(outdoc,
                                       proc_id, ifo, "CBC_WORKFLOW_SCIENCE_OK", 0)
            segmentdb_utils.add_to_segment(outdoc, proc_id, sci_ok_def_id,
                                                                   sci_ok_segs)
            segmentdb_utils.add_to_segment_summary(outdoc, proc_id,
                                          sci_ok_def_id, summ_segs, comment='')
        elif sci_ok_seg_file:
            # FIXME: While the segment module is still fractured (#127) this
            #        may not work. Please update when #127 is resolved
            pass
            #err_msg = "Got %d files matching %s and %s. Expected 1 or 0." \
            #          %(len(sci_ok_seg_file), ifo, 'SCIENCE_OK')
            #raise ValueError(err_msg)


        # SCIENCE_AVAILABLE
        sci_available_seg_file = seg_ifo_files.find_output_with_tag(\
                                                           'SCIENCE_AVAILABLE')
        if len(sci_available_seg_file) == 1:
            sci_available_seg_file = sci_available_seg_file[0]
            sci_available_segs = sci_available_seg_file.segmentList
            sci_available_def_id = segmentdb_utils.add_to_segment_definer(\
                        outdoc, proc_id, ifo, "CBC_WORKFLOW_SCIENCE_AVAILABLE", 0)
            segmentdb_utils.add_to_segment(outdoc, proc_id,
                                      sci_available_def_id, sci_available_segs)
            segmentdb_utils.add_to_segment_summary(outdoc, proc_id,
                                   sci_available_def_id, summ_segs, comment='')
        elif sci_available_seg_file:
            # FIXME: While the segment module is still fractured (#127) this
            #        may not work. Please update when #127 is resolved
            pass
            #err_msg = "Got %d files matching %s and %s. Expected 1 or 0." \
            #          %(len(sci_available_seg_file), ifo, 'SCIENCE_AVAILABLE')
            #raise ValueError(err_msg)

        # ANALYSABLE - This one needs to come from inspiral outs
        ifo_insps = insps.find_output_with_ifo(ifo)
        analysable_segs = ifo_insps.get_times_covered_by_files()

        analysable_def_id = segmentdb_utils.add_to_segment_definer(outdoc,
                                     proc_id, ifo, "CBC_WORKFLOW_ANALYSABLE", 0)
        segmentdb_utils.add_to_segment(outdoc, proc_id, analysable_def_id,
                                                               analysable_segs)
        segmentdb_utils.add_to_segment_summary(outdoc, proc_id,
                                      analysable_def_id, summ_segs, comment='')

    summ_file = File(workflow.ifos, "WORKFLOW_SUMMARY",
                     workflow.analysis_time, extension=".xml",
                     directory=output_dir)
    summ_file.PFN(summ_file.storage_path, site='local')
    utils.write_filename(outdoc, summ_file.storage_path)

    return FileList([summ_file])

