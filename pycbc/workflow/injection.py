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
This module is responsible for setting up the part of a pycbc workflow that
will generate the injection files to be used for assessing the workflow's
ability to detect predicted signals. (In ihope parlance, this sets up the
inspinj jobs). Full documentation for this module can be found here:
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/NOTYETCREATED.html
"""

import os
import logging
import urllib
import pycbc.workflow.core
import pycbc.workflow.jobsetup
from glue import segments

def setup_injection_workflow(workflow, output_dir=None,
                             injSectionName='injections', tags =[]):
    '''
    This function is the gateway for setting up injection-generation jobs in a
    workflow. It should be possible for this function to support a number
    of different ways/codes that could be used for doing this, however as this
    will presumably stay as a single call to a single code (which need not be
    inspinj) there are currently no subfunctions in this moudle. 

    Parameters
    -----------
    workflow : pycbc.workflow.core.Workflow
        The Workflow instance that the coincidence jobs will be added to.
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
    inj_files : pycbc.workflow.core.FileList
        The list of injection files created by this call.
    inj_tags : list of strings
        The tag corresponding to each injection file and used to uniquely
        identify them. The FileList class contains functions to search
        based on tags.
    '''
    logging.info("Entering injection module.")
    make_analysis_dir(output_dir)
    # Get full analysis segment for output file naming
    fullSegment = workflow.analysis_time

    # Identify which injections to do by presence of sub-sections in
    # the configuration file
    all_sec = workflow.cp.sections()
    sections = [sec for sec in all_sec if sec.startswith(injSectionName +'-')]

    inj_tags = []
    inj_files = pycbc.workflow.core.FileList([])   

    for section in sections:
        split_sec_name = section.split('-')
        # Sanity check the user has realised that the '-' is a special char.
        if len(split_sec_name) > 2:
            # This is unusual, but a format [injection-name-tag] is okay. Just
            # check that [injection-name] section exists. If not it is possible
            # the user is trying to use an injection name with '-' in it
            sect_check = "%s-%s" %(split_sec_name[0], split_sec_name[1])
            if not workflow.cp.has_section(sect_check):
                err_msg = "Injection section found with name %s. " %(section,)
                err_msg += "Workflow uses the '-' as a delimiter so this is "
                err_msg += "interpreted as exe_tag-inj_tag-other_tag. "
                err_msg += "No section with name %s was found. " %(sect_check,)
                err_msg += "If you did not intend to use tags in an "
                err_msg += "'advanced user' manner, or do not understand what "
                err_msg += "this means, don't use dashes in injection "
                err_msg += "names. So [injection-nsbhinj] is good. "
                err_msg += "[injection-nsbh-inj] is not."
                raise ValueError(err_msg)
            continue
 
        inj_tag = (split_sec_name[1]).upper()
        currTags = tags + [inj_tag]

        # FIXME: Remove once fixed in pipedown
        # TEMPORARILY we require inj tags to end in "INJ"
        if not inj_tag.endswith("INJ"):
            err_msg = "Currently workflow requires injection names to end with "
            err_msg += "a inj suffix. Ie. bnslininj or bbhinj. "
            err_msg += "%s is not good." %(inj_tag.lower())
            raise ValueError(err_msg)

        # Parse for options in ini file
        injectionMethod = workflow.cp.get_opt_tags("workflow-injections", 
                                                 "injections-method", currTags)

        if injectionMethod in ["IN_WORKFLOW", "AT_RUNTIME"]:
            # FIXME: Add ability to specify different exes
            inj_job = pycbc.workflow.jobsetup.LalappsInspinjExecutable(workflow.cp, injSectionName, tags=currTags,
                                         out_dir=output_dir, ifos='HL')
            node = inj_job.create_node(fullSegment)
            if injectionMethod == "AT_RUNTIME":
                workflow.execute_node(node)
            else:
                workflow.add_node(node)
            injFile = node.output_files[0]
        elif injectionMethod == "PREGENERATED":
            injectionFilePath = workflow.cp.get_opt_tags("workflow-injections",
                                      "injections-pregenerated-file", currTags)
            file_url = urlparse.urljoin('file:', urllib.pathname2url(\
                                                  injectionFilePath))
            injFile = pycbc.workflow.core.File('HL', 'PREGEN_INJFILE', fullSegment, file_url,
                                tags=currTags)
            injFile.PFN(injectionFilePath, site='local')
        else:
            errMsg = "Injection method must be one of IN_WORKFLOW, "
            errMsg += "AT_RUNTIME or PREGENERATED. Got %s." %(injectionMethod)
            raise ValueError(errMsg)

        inj_files.append(injFile)
        inj_tags.append(inj_tag)
    logging.info("Leaving injection module.")
    return inj_files, inj_tags

