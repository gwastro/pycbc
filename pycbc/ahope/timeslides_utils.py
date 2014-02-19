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
This module is responsible for setting up the time slide files for ahope
workflows. For details about this module and its capabilities see here:
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/NOTYETCREATED.html
"""

import os
import logging
import urllib
from glue import segments
from pycbc.ahope.ahope_utils import *

def setup_timeslides_workflow(workflow, output_dir=None, tags=[],
                              timeSlideSectionName='ligolw_tisi'):
    '''
    Setup generation of time_slide input files in the ahope workflow.
    Currently used
    only with ligolw_tisi to generate files containing the list of slides to be
    performed in each time slide job.

    Parameters
    -----------
    Workflow : ahope.Workflow
        The ahope workflow instance that the coincidence jobs will be added to.
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
    timeSlideOuts : ahope.AhopeFileList
        The list of time slide files created by this call.
    '''
    logging.info("Entering time slides setup module.")
    make_analysis_dir(output_dir)
    # Get ifo list and full analysis segment for output file naming
    ifoList = workflow.ifos
    ifoString = workflow.ifoString
    fullSegment = workflow.analysis_time

    # Identify which time-slides to do by presence of sub-sections in the
    # configuration file
    all_sec = workflow.cp.sections()
    timeSlideSections = [sec for sec in all_sec if sec.startswith('tisi-')]
    timeSlideTags = [sec.split('-')[-1] for sec in timeSlideSections]

    timeSlideOuts = AhopeFileList([])

    # FIXME: Add ability to specify different exes

    # Make the timeSlideFiles
    for timeSlideTag in timeSlideTags:
        currTags = tags + [timeSlideTag]

        timeSlideMethod = workflow.cp.get_opt_tags("ahope-timeslides",
                                                 "timeslides-method", currTags)

        if timeSlideMethod in ["IN_WORKFLOW", "AT_RUNTIME"]:
            # FIXME: Currently hacks because ligolw_tisi's options cannot be
            # supplied cleanly from a config file. This needs fixing and this
            # part cleaned up. This will then allow AT_WORKFLOW as a method
            if timeSlideMethod == "IN_WORKFLOW":
                errMsg = "Cannot yet put ligolw_tisi into the workflow. "
                errMsg += "Please see the item in the ahope to-do list."
                errMsg += " Try method=AT_RUNTIME instead for now."
                raise ValueError(errMsg)
            # FIXME: Add ability to specify different exes

            # First we need to run ligolw_tisi to make the necessary time slide
            # input xml files
            tisiOutFile = AhopeFile(ifoString, 'TIMESLIDES', fullSegment,
                                    directory=output_dir, extension=".xml.gz",
                                    tags=currTags)
            ligolw_tisi_call = [workflow.cp.get('executables', 'tisi'), "-v"]
            subString = 'tisi-%s' %(timeSlideTag.lower())
            if workflow.cp.has_option('tisi', 'inspiral-num-slides'):
                ligolw_tisi_call.append("--inspiral-num-slides")
                ligolw_tisi_call.append(\
                        workflow.cp.get('tisi', 'inspiral-num-slides'))
            elif workflow.cp.has_option(subString, 'inspiral-num-slides'):
                ligolw_tisi_call.append("--inspiral-num-slides")
                ligolw_tisi_call.append(\
                        workflow.cp.get(subString, 'inspiral-num-slides'))
            else:
                for ifo in ifoList:
                    ifoSlideStart = workflow.cp.get_opt_tag('tisi',
                                        '%s-slide-start' %(ifo), timeSlideTag)
                    ifoSlideEnd = workflow.cp.get_opt_tag('tisi',
                                        '%s-slide-end' %(ifo), timeSlideTag)
                    ifoSlideStep = workflow.cp.get_opt_tag('tisi',
                                        '%s-slide-step' %(ifo), timeSlideTag)
                    ligolw_tisi_call.append("-i")
                    optionString = ':'.join([ifoSlideStart, ifoSlideEnd, \
                                             ifoSlideStep])
                    optionString = '%s=%s' %(ifo.upper(), optionString)
                    ligolw_tisi_call.append(optionString)
            if workflow.cp.has_option('tisi', 'remove-zero-lag') or\
                       workflow.cp.has_option(subString, 'remove-zero-lag'):
                ligolw_tisi_call.append("--remove-zero-lag")
            ligolw_tisi_call.append(tisiOutFile.path)
            make_external_call(ligolw_tisi_call,
                            outDir=os.path.join(output_dir,'logs'),
                            outBaseName='%s-ligolw_tisi-call' %(timeSlideTag) )
        elif timeSlideMethod == "PREGENERATED":
            timeSlideFilePath = workflow.cp.get_opt_tags("ahope-timeslides",
                                      "timeslides-pregenerated-file", currTags)
            file_url = urlparse.urljoin('file:', urllib.pathname2url(\
                                                  timeSlideFilePath))
            tisiOutFile = AhopeFile(ifoString, 'PREGEN_TIMESLIDES',
                                    fullSegment, file_url, tags=currTags)

        timeSlideOuts.append(tisiOutFile)

    return timeSlideOuts
