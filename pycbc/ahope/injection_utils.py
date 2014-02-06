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
This module is responsible for setting up the part of an ahope workflow that
will generate the injection files to be used for assessing the workflow's
ability to detect predicted signals. (In ihope parlance, this sets up the
inspinj jobs). Full documentation for this module can be found here:
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/NOTYETCREATED.html
"""

import os
import logging
import pycbc.ahope
from pycbc.ahope.jobsetup_utils import *
from pycbc.ahope.matchedfltr_utils import *
from glue import segments

def setup_injection_workflow(workflow, science_segs, output_dir=None,
                             injSectionName='injections', tags =[]):
    '''
    This function is the gateway for setting up injection-generation jobs in an
    ahope workflow. It should be possible for this function to support a number
    of different ways/codes that could be used for doing this, however as this
    will presumably stay as a single call to a single code (which need not be
    inspinj) there are currently no subfunctions in this moudle. 

    Parameters
    -----------
    Workflow : ahope.Workflow
        The ahope workflow instance that the coincidence jobs will be added to.
    science_segs : ifo-keyed dictionary of glue.segments.segmentlist instances
        The list of times that are being analysed in this workflow. 
    output_dir : path
        The directory in which injection files will be stored.
    injSectionName : string (optional, default='injections')
        The string that corresponds to the option describing the exe location
        in the [executables] section of the .ini file and that corresponds to
        the section (and sub-sections) giving the options that will be given to
        the code at run time.
    tags : list of strings (optional, default = [])
        A list of the tagging strings that will be used for all jobs created
        by this call to the workflow. This will be used in output names.

    Returns
    --------
    inj_files : ahope.AhopeFileList
        The list of injection files created by this call.
    inj_tags : list of strings
        The tag corresponding to each injection file and used to uniquely
        identify them. The AhopeFileList class contains functions to search
        based on tags.
    '''
    logging.info("Entering injection module.")
    make_analysis_dir(output_dir)

    # Get full analysis segment for output file naming
    fullSegment = pycbc.ahope.get_full_analysis_chunk(science_segs)

    # FIXME: Add ability to specify different exes
    inj_exe = LalappsInspinjExec(injSectionName)
    all_sec = workflow.cp.sections()
    sections = [sec for sec in all_sec if sec.startswith(injSectionName +'-')]
    inj_tags = []
    inj_files = []
    for sec in sections:  
        inj_tag = sec.split('-')[1]
        curr_tags = tags + [inj_tag]
        inj_job = inj_exe.create_job(workflow.cp, tags=curr_tags,
                                     out_dir=output_dir)
        node = inj_job.create_node(fullSegment)
        workflow.add_node(node)
        injection_file = node.output_files[0]
        inj_files.append(injection_file)
        inj_tags.append(inj_tag)
    logging.info("Leaving injection module.")  
    return inj_files, inj_tags

