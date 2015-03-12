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
import Pegasus.DAX3 as dax
from glue import segments
from pycbc.workflow.core import File, FileList, make_analysis_dir
from pycbc.workflow.jobsetup import select_generic_executable
from pycbc.workflow.legacy_ihope import LegacyCohPTFTrigCombiner 

def setup_coh_PTF_post_processing(workflow, trigger_files, injection_files,
                                  output_dir, segment_dir, run_dir=None,
                                  ifos=None, tags=[], **kwargs):
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
    logging.info("Entering post-processing stage.")
    make_analysis_dir(output_dir)

    # Parse for options in .ini file
    post_proc_method = workflow.cp.get_opt_tags("workflow-postproc",
                                                "postproc-method", tags)

    # Scope here for adding different options/methods here. For now we only
    # have the single_stage ihope method which consists of converting the
    # ligolw_thinca output xml into one file, clustering, performing injection
    # finding and putting everything into one SQL database.
    if post_proc_method == "COH_PTF_WORKFLOW":
        post_proc_files = setup_postproc_coh_PTF_workflow(workflow,
                           trigger_files, injection_files, output_dir,
                           segment_dir, ifos=ifos, tags=tags, **kwargs)
    else:
        errMsg = "Post-processing method not recognized. Must be "
        errMsg += "COH_PTF_WORKFLOW."
        raise ValueError(errMsg)

    logging.info("Leaving post-processing module.")

    return post_proc_files


def setup_old_coh_PTF_pp_workflow(workflow, analysis_files, injection_files,
                                  output_dir, segment_dir, run_dir, ifos,
                                  tags=[]):
    """
    Set up a job to run coh_PTF_post_processing and submit the resulting DAG.
    """
    cp = workflow.cp
    post_proc_outs = FileList([])

    coh_ptf_pp_exe = os.path.basename(cp.get("executables", "post_processing"))
    coh_ptf_pp_class = select_generic_executable(workflow, "post_processing")

    coh_ptf_pp_job = coh_ptf_pp_class(workflow.cp, "post_processing",
                                      ifo=ifos,
                                      injection_file=injection_files,
                                      out_dir=output_dir, tags=tags)
    coh_ptf_pp_node = coh_ptf_pp_job.create_node(analysis_files, run_dir,
                                                 segment_dir, output_dir)
    workflow.add_node(coh_ptf_pp_node)

    return post_proc_outs


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

    # Set up needed exe classes
    trig_combiner_exe = os.path.basename(cp.get("executables",
                                                "trig_combiner"))
    trig_combiner_class = select_generic_executable(workflow,
                                                    "trig_combiner")

    trig_cluster_exe = os.path.basename(cp.get("executables",
                                               "trig_cluster"))
    trig_cluster_class = select_generic_executable(workflow,
                                                   "trig_cluster")

    sbv_plotter_exe = os.path.basename(cp.get("executables",
                                              "sbv_plotter"))
    sbv_plotter_class = select_generic_executable(workflow,
                                                  "sbv_plotter")
    """
    efficiency_exe = os.path.basename(cp.get("executables",
                                             "efficiency"))
    efficiency_class = select_generic_executable(workflow,
                                                 "efficiency")

    horizon_dist_exe = os.path.basename(cp.get("executables",
                                               "horizon_dist"))
    horizon_dist_class = select_generic_executable(workflow,
                                                   "horizon_dist")
    """
    # Set up trig_combiner job
    trig_combiner_jobs = trig_combiner_class(workflow.cp, "trig_combiner",
                                             ifo=ifos, 
                                             injection_file=injection_files,
                                             out_dir=output_dir, tags=tags)
    trig_combiner_node = trig_combiner_jobs.create_node(trigger_files,
                                                        segment_dir, tags=tags)
    workflow.add_node(trig_combiner_node)
    #FIXME: trig_combiner_node.output_files is empty! This is a hack.
    #trig_combiner_outs = trig_combiner_node.output_files
    #post_proc_outs.append(trig_combiner_outs)

    # Initialise trig_cluster class
    trig_cluster_outs = FileList([])
    trig_cluster_jobs = trig_cluster_class(workflow.cp, "trig_cluster",
                                           ifo=ifos,
                                           out_dir=output_dir, tags=tags)

    # Initialise sbv_plotter class
    sbv_plotter_outs = FileList([])
    sbv_plotter_jobs = sbv_plotter_class(workflow.cp, "sbv_plotter", ifo=ifos,
                                         out_dir=output_dir, tags=tags)

    # Add trig_cluster jobs and their corresponding plotting jobs
    trig_name = cp.get("workflow", "trigger-name")
    cluster_out_tags = ["ALL_TIMES", "OFFSOURCE", "ONSOURCE"]

    for out_tag in cluster_out_tags:
        trig_cluster_node, \
        unclustered_trigs = trig_cluster_jobs.create_node(trigger_files[-1],
                tags=["GRB%s_%s" % (trig_name, out_tag)])
        trig_cluster_outs.extend(trig_cluster_node.output_files)

        if out_tag != "ONSOURCE":
            # Add memory requirememnt for jobs with potentially large files
            trig_cluster_node.set_memory(1300)
            workflow.add_node(trig_cluster_node)
            dep = dax.Dependency(parent=trig_combiner_node._dax_node,
                                 child=trig_cluster_node._dax_node)
            workflow._adag.addDependency(dep)

            # Add sbv_plotter job
            parent = trig_cluster_node.output_files[0]
            sbv_out_tags = [out_tag, "_clustered"]
            sbv_plotter_node = sbv_plotter_jobs.create_node(parent,
                                                            segment_dir,
                                                            tags=sbv_out_tags)
            workflow.add_node(sbv_plotter_node)

            # Also add sbv_plotter job for unclustered triggers
            parent = unclustered_trigs
            sbv_out_tags = [out_tag, "_unclustered"]
            sbv_plotter_node = sbv_plotter_jobs.create_node(parent,
                                                            segment_dir,
                                                            tags=sbv_out_tags)
            sbv_plotter_node.set_memory(1300)
            workflow.add_node(sbv_plotter_node)
            dep = dax.Dependency(parent=trig_combiner_node._dax_node,
                                 child=sbv_plotter_node._dax_node)
            workflow._adag.addDependency(dep)
        else:
            workflow.add_node(trig_cluster_node)
            dep = dax.Dependency(parent=trig_combiner_node._dax_node,
                                 child=trig_cluster_node._dax_node)
            workflow._adag.addDependency(dep)

    # Add further trig_cluster jobs for trials
    trials = int(cp.get("trig_combiner", "num-trials"))
    trial = 1

    while trial <= trials:
        trig_cluster_node, \
        unclustered_trigs = trig_cluster_jobs.create_node(trigger_files[-1],
                tags=["GRB%s_OFFTRIAL_%d" % (trig_name, trial)])
        trig_cluster_outs.extend(trig_cluster_node.output_files)
        workflow.add_node(trig_cluster_node)
        dep = dax.Dependency(parent=trig_combiner_node._dax_node,
                             child=trig_cluster_node._dax_node)
        workflow._adag.addDependency(dep)
        trial += 1

    post_proc_outs.extend(trig_cluster_outs)
    
    return post_proc_outs
