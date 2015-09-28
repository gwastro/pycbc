# Copyright (C) 2015  Ian Harry, Alex Nitz
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

import logging, urllib, urlparse
from pycbc.workflow.core import File, FileList, make_analysis_dir, Executable
from pycbc.workflow.jobsetup import (LalappsInspinjExecutable,
        LigolwCBCJitterSkylocExecutable, LigolwCBCAlignTotalSpinExecutable,
        PycbcDarkVsBrightInjectionsExecutable)

def veto_injections(workflow, inj_file, veto_file, veto_name, out_dir, tags=None):
    tags = [] if tags is None else tags
    make_analysis_dir(out_dir)
    
    node = Executable(workflow.cp, 'strip_injections', ifos=workflow.ifos,
                          out_dir=out_dir, tags=tags).create_node()
    node.add_opt('--segment-name', veto_name)
    node.add_input_opt('--veto-file', veto_file)
    node.add_input_opt('--injection-file', inj_file)
    node.add_opt('--ifos', ' '.join(workflow.ifos))
    node.new_output_file_opt(workflow.analysis_time, '.xml', '--output-file')
    workflow += node
    return node.output_files[0]  

def setup_injection_workflow(workflow, output_dir=None,
                             inj_section_name='injections', exttrig_file=None,
                             tags =[]):
    """
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
    inj_section_name : string (optional, default='injections')
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
    """
    logging.info("Entering injection module.")
    make_analysis_dir(output_dir)
    
    # Get full analysis segment for output file naming
    full_segment = workflow.analysis_time
    ifos = workflow.ifos

    # Identify which injections to do by presence of sub-sections in
    # the configuration file
    inj_tags = []
    inj_files = FileList([])  

    for section in  workflow.cp.get_subsections(inj_section_name):
        inj_tag = section.upper()
        curr_tags = tags + [inj_tag]

        # FIXME: Remove once fixed in pipedown
        # TEMPORARILY we require inj tags to end in "INJ"
        if not inj_tag.endswith("INJ"):
            err_msg = "Currently workflow requires injection names to end with "
            err_msg += "a inj suffix. Ie. bnslininj or bbhinj. "
            err_msg += "%s is not good." %(inj_tag.lower())
            raise ValueError(err_msg)

        # Parse for options in ini file
        injection_method = workflow.cp.get_opt_tags("workflow-injections", 
                                                    "injections-method",
                                                    curr_tags)

        if injection_method in ["IN_WORKFLOW", "AT_RUNTIME"]:
            # FIXME: Add ability to specify different exes
            inj_job = LalappsInspinjExecutable(workflow.cp, inj_section_name,
                                               out_dir=output_dir, ifos='HL',
                                               tags=curr_tags)
            node = inj_job.create_node(full_segment)
            if injection_method == "AT_RUNTIME":
                workflow.execute_node(node)
            else:
                workflow.add_node(node)
            inj_file = node.output_files[0]
            inj_files.append(inj_file)
        elif injection_method == "PREGENERATED":
            injectionFilePath = workflow.cp.get_opt_tags("workflow-injections",
                                      "injections-pregenerated-file", curr_tags)
            file_url = urlparse.urljoin('file:',
                                        urllib.pathname2url(injectionFilePath))
            inj_file = File('HL', 'PREGEN_inj_file', full_segment, file_url,
                            tags=curr_tags)
            inj_file.PFN(injectionFilePath, site='local')
            inj_files.append(inj_file)
        elif injection_method in ["IN_COH_PTF_WORKFLOW", "AT_COH_PTF_RUNTIME"]:
            inj_job = LalappsInspinjExecutable(workflow.cp, inj_section_name,
                                               out_dir=output_dir, ifos=ifos,
                                               tags=curr_tags)
            node = inj_job.create_node(full_segment, exttrig_file)
            if injection_method == "AT_COH_PTF_RUNTIME":
                workflow.execute_node(node)
            else:
                workflow.add_node(node)
            inj_file = node.output_files[0]

            if workflow.cp.has_option("workflow-injections",
                                      "em-bright-only"):
                em_filter_job = PycbcDarkVsBrightInjectionsExecutable(
                                                 workflow.cp,
                                                 'em_bright_filter',
                                                 tags=curr_tags,
                                                 out_dir=output_dir,
                                                 ifos=ifos)
                node = em_filter_job.create_node(inj_file, full_segment,
                                                 curr_tags)
                if injection_method == "AT_COH_PTF_RUNTIME":
                    workflow.execute_node(node)
                else:
                    workflow.add_node(node)
                inj_file = node.output_files[0] 

            if workflow.cp.has_option("workflow-injections",
                                      "do-jitter-skyloc"):
                jitter_job = LigolwCBCJitterSkylocExecutable(workflow.cp, 
                                                             'jitter_skyloc',
                                                             tags=curr_tags,
                                                             out_dir=output_dir,
                                                             ifos=ifos)
                node = jitter_job.create_node(inj_file, full_segment, curr_tags)
                if injection_method == "AT_COH_PTF_RUNTIME":
                    workflow.execute_node(node)
                else:
                    workflow.add_node(node)
                inj_file = node.output_files[0]
            
            if workflow.cp.has_option("workflow-injections",
                                      "do-align-total-spin"):
                align_job = LigolwCBCAlignTotalSpinExecutable(workflow.cp,
                        'align_total_spin', tags=curr_tags, out_dir=output_dir,
                        ifos=ifos)
                node = align_job.create_node(inj_file, full_segment, curr_tags)
                
                if injection_method == "AT_COH_PTF_RUNTIME":
                    workflow.execute_node(node)
                else:
                    workflow.add_node(node)
                inj_file = node.output_files[0]
            
            inj_files.append(inj_file)
        else:
            err = "Injection method must be one of IN_WORKFLOW, "
            err += "AT_RUNTIME or PREGENERATED. Got %s." % (injection_method)
            raise ValueError(err)

        inj_tags.append(inj_tag)
        
    logging.info("Leaving injection module.")
    return inj_files, inj_tags

