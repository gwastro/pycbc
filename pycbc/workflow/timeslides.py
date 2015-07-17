# Copyright (C) 2013  Ian Harry
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
This module is responsible for setting up the time slide files for
workflows. For details about this module and its capabilities see here:
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/NOTYETCREATED.html
"""

import logging
import urllib
from pycbc.workflow.core import File, FileList, make_analysis_dir
from pycbc.workflow.jobsetup import select_generic_executable

def setup_timeslides_workflow(workflow, output_dir=None, tags=[],
                              timeSlideSectionName='ligolw_tisi'):
    '''
    Setup generation of time_slide input files in the workflow.
    Currently used
    only with ligolw_tisi to generate files containing the list of slides to be
    performed in each time slide job.

    Parameters
    -----------
    workflow : pycbc.workflow.core.Workflow
        The Workflow instance that the coincidence jobs will be added to.
    output_dir : path
        The directory in which output files will be stored.
    tags : list of strings (optional, default = [])
        A list of the tagging strings that will be used for all jobs created
        by this call to the workflow. This will be used in output names.
    timeSlideSectionName : string (optional, default='injections')
        The string that corresponds to the option describing the exe location
        in the [executables] section of the .ini file and that corresponds to
        the section (and sub-sections) giving the options that will be given to
        the code at run time.
    Returns
    --------
    timeSlideOuts : pycbc.workflow.core.FileList
        The list of time slide files created by this call.
    '''
    logging.info("Entering time slides setup module.")
    make_analysis_dir(output_dir)
    # Get ifo list and full analysis segment for output file naming
    ifoList = workflow.ifos
    ifo_string = workflow.ifo_string
    fullSegment = workflow.analysis_time

    # Identify which time-slides to do by presence of sub-sections in the
    # configuration file
    all_sec = workflow.cp.sections()
    timeSlideSections = [sec for sec in all_sec if sec.startswith('tisi-')]
    timeSlideTags = [(sec.split('-')[-1]).upper() for sec in timeSlideSections]

    timeSlideOuts = FileList([])

    # FIXME: Add ability to specify different exes

    # Make the timeSlideFiles
    for timeSlideTag in timeSlideTags:
        currTags = tags + [timeSlideTag]

        timeSlideMethod = workflow.cp.get_opt_tags("workflow-timeslides",
                                                 "timeslides-method", currTags)

        if timeSlideMethod in ["IN_WORKFLOW", "AT_RUNTIME"]:
            timeSlideExeTag = workflow.cp.get_opt_tags("workflow-timeslides",
                                                    "timeslides-exe", currTags)
            timeSlideExe = select_generic_executable(workflow, timeSlideExeTag)
            timeSlideJob = timeSlideExe(workflow.cp, timeSlideExeTag, ifos=ifo_string,
                                             tags=currTags, out_dir=output_dir)
            timeSlideNode = timeSlideJob.create_node(fullSegment)
            if timeSlideMethod == "AT_RUNTIME":
                workflow.execute_node(timeSlideNode)
            else:
                workflow.add_node(timeSlideNode)
            tisiOutFile = timeSlideNode.output_files[0]
        elif timeSlideMethod == "PREGENERATED":
            timeSlideFilePath = workflow.cp.get_opt_tags("workflow-timeslides",
                                      "timeslides-pregenerated-file", currTags)
            file_url = urlparse.urljoin('file:', urllib.pathname2url(\
                                                  timeSlideFilePath))
            tisiOutFile = File(ifoString, 'PREGEN_TIMESLIDES',
                               fullSegment, file_url, tags=currTags)

        timeSlideOuts.append(tisiOutFile)

    return timeSlideOuts
