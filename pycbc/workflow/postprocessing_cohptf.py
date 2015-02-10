# Copyright (C) 2015 Andrew Williamson
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
This module is used to prepare output from the coherent matched filter workflow
for post-processing (post-processing = calculating significance of candidates
and making any rates statements).
For details of this module and its capabilities see here:
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/NOTYETCREATED.html
"""

from __future__ import division

import os
import os.path
import logging
from glue import segments
from pycbc.workflow.core import FileList, make_analysis_dir
from pycbc.workflow.jobsetup import select_generic_executable
from pycbc.workflow.legacy_ihope import LegacyCohPTFTrigCombiner 

def setup_coh_PTF_postprocessing(workflow, trigger_files, injection_cache,
                                 output_dir, segment_dir, ifos=None, tags=[],
                                 **kwargs):
    """
    This function aims to be the gateway for running postprocessing in CBC
    offline workflows. Post-processing generally consists of calculating the
    significance of triggers and making any statements about trigger rates.
    Dedicated plotting jobs do not belong here.

    Properties
    -----------
    workflow : pycbc.workflow.core.Workflow
        The Workflow instance that the coincidence jobs will be added to.
    trigger_files : pycbc.workflow.core.FileList
        An FileList of the trigger files that are used as
        input at this stage.
    summary_xml_files : pycbc.workflow.core.FileList
        An FileList of the output of the analysislogging_utils module.
    output_dir : path
        The directory in which output files will be stored.
    tags : list of strings (optional, default = [])
        A list of the tagging strings that will be used for all jobs created
        by this call to the workflow. An example might be ['POSTPROC1'] or
        ['DENTYSNEWPOSTPROC']. This will be used in output names.

    Returns
    --------
    post_proc_files : pycbc.workflow.core.FileList
        A list of the output from this stage.

    """
    logging.info("Entering post-processing preparation stage.")
    make_analysis_dir(output_dir)

    # Parse for options in .ini file
    post_proc_method = workflow.cp.get_opt_tags("workflow-postproc",
                                                "postproc-method", tags)

    # Scope here for adding different options/methods here. For now we only
    # have the single_stage ihope method which consists of converting the
    # ligolw_thinca output xml into one file, clustering, performing injection
    # finding and putting everything into one SQL database.
    if post_proc_method == "PIPEDOWN_WORKFLOW":
        # If you want the intermediate output files, call this directly
        post_proc_files = setup_postproc_pipedown_workflow(workflow,
                           trigger_files, summary_xml_files, output_dir,
                           tags=tags, **kwargs)
    if post_proc_method == "COH_PTF_WORKFLOW":
        post_proc_files = setup_postproc_coh_PTF_workflow(workflow,
                           trigger_files, injection_cache, output_dir,
                           segment_dir, ifos=ifos, tags=tags, **kwargs)
    else:
        errMsg = "Post-processing method not recognized. Must be "
        errMsg += "one of PIPEDOWN_WORKFLOW or COH_PTF_WORKFLOW."
        raise ValueError(errMsg)

    logging.info("Leaving post-processing module.")

    return post_proc_files


def setup_postproc_coh_PTF_workflow(workflow, trigger_files, injection_files,
                                    output_dir, segment_dir, ifos,
                                    tags=[]):
    """
    This module sets up the post-processing stage in the workflow, using a
    coh_PTF style set up. This consists of running trig_combiner to find
    coherent triggers, and injfinder to look for injections. It then runs
    a horizon_dist job, trig_cluster to cluster triggers, and injcombiner to
    calculate injection statistics. Finally, efficiency and sbv_plotter jobs
    calculate efficiency and signal based veto statistics and make plots.
    
    workflow : pycbc.workflow.core.Workflow
        The Workflow instance that the jobs will be added to.
    trigger_files : pycbc.workflow.core.FileList
        A FileList containing the combined databases.
   
    Returns
    --------
    
    """
    cp = workflow.cp
    post_proc_outs = FileList([])

    # Setup needed exe classes
    trig_combiner_exe = os.path.basename(cp.get("executables",
                                                "trig_combiner"))
    trig_combiner_class = select_generic_executable(workflow,
                                                    "trig_combiner")

    """
    trig_cluster_exe = workflow.cp.get("executables", "trig_cluster")
    sbv_plotter_exe = workflow.cp.get("executables", "sbv_plotter")
    efficiency_exe = workflow.cp.get("executables", "efficiency")
    horizon_dist = workflow.cp.get("executables", "horizon_dist")
    """

    # Setup trig_combiner job
    trig_combiner_jobs = trig_combiner_class(workflow.cp, 'trig_combiner',
                                             ifo=ifos, 
                                             injection_file=injection_files,
                                             out_dir=output_dir, tags=tags)
    trig_combiner_node = trig_combiner_jobs.create_node(trigger_files,
                                                        segment_dir)
    workflow.add_node(trig_combiner_node)
    post_proc_outs.append(trig_combiner_node.output_files)

    return post_proc_outs
