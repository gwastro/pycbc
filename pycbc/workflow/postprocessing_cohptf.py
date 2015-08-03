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

def setup_coh_PTF_post_processing(workflow, trigger_files, output_dir,
                                  segment_dir, injection_trigger_files=None,
                                  injection_files=None,
                                  injection_trigger_caches=None,
                                  injection_caches=None, config_file=None,
                                  run_dir=None, ifos=None, inj_tags=[],
                                  tags=[], **kwargs):
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
                trigger_files, injection_trigger_files, injection_files,
                injection_trigger_caches, injection_caches, config_file,
                output_dir, segment_dir, ifos=ifos, inj_tags=inj_tags,
                tags=tags, **kwargs)
    else:
        errMsg = "Post-processing method not recognized. Must be "
        errMsg += "COH_PTF_WORKFLOW."
        raise ValueError(errMsg)

    logging.info("Leaving post-processing module.")

    return post_proc_files


def setup_postproc_coh_PTF_workflow(workflow, trig_files, inj_trig_files,
                                    inj_files, inj_trig_caches, inj_caches,
                                    config_file, output_dir, segment_dir,
                                    ifos, inj_tags=[], tags=[]):
    """
    This module sets up the post-processing stage in the workflow, using a
    coh_PTF style set up. This consists of running trig_combiner to find
    coherent triggers, and injfinder to look for injections. It then runs
    a horizon_dist job, trig_cluster to cluster triggers, and injcombiner to
    calculate injection statistics. Finally, efficiency and sbv_plotter jobs
    calculate efficiency and signal based veto statistics and make plots.
    
    workflow : pycbc.workflow.core.Workflow
        The Workflow instance that the jobs will be added to.
    trig_files : pycbc.workflow.core.FileList
        A FileList containing the combined databases.
   
    Returns
    --------
    
    """
    cp = workflow.cp
    full_segment = trig_files[0].segment
    post_proc_outs = FileList([])
    pp_nodes = []

    # Set up needed exe classes
    trig_combiner_exe = os.path.basename(cp.get("executables",
                                                "trig_combiner"))
    trig_combiner_class = select_generic_executable(workflow, "trig_combiner")

    trig_cluster_exe = os.path.basename(cp.get("executables", "trig_cluster"))
    trig_cluster_class = select_generic_executable(workflow, "trig_cluster")


    sbv_plotter_exe = os.path.basename(cp.get("executables", "sbv_plotter"))
    sbv_plotter_class = select_generic_executable(workflow, "sbv_plotter")
    
    efficiency_exe = os.path.basename(cp.get("executables", "efficiency"))
    efficiency_class = select_generic_executable(workflow, "efficiency")
    """
    horizon_dist_exe = os.path.basename(cp.get("executables",
                                               "horizon_dist"))
    horizon_dist_class = select_generic_executable(workflow,
                                                   "horizon_dist")
    """
    html_summary_exe = os.path.basename(cp.get("executables", "html_summary"))
    html_summary_class = select_generic_executable(workflow, "html_summary")

    # Set up trig_combiner job
    trig_combiner_jobs = trig_combiner_class(cp, "trig_combiner", ifo=ifos, 
                                             out_dir=output_dir, tags=tags)
    trig_combiner_node = trig_combiner_jobs.create_node(trig_files,
                                                        segment_dir, tags=tags)
    pp_nodes.append(trig_combiner_node)
    workflow.add_node(trig_combiner_node)
    #FIXME: trig_combiner_node.output_files is empty! This is a hack.
    #trig_combiner_outs = trig_combiner_node.output_files
    #post_proc_outs.append(trig_combiner_outs)

    # Initialise trig_cluster class
    trig_cluster_outs = FileList([])
    trig_cluster_jobs = trig_cluster_class(cp, "trig_cluster", ifo=ifos,
                                           out_dir=output_dir, tags=tags)

    # Set up injfinder jobs
    if cp.has_section("workflow-injections"):
        injfinder_nodes = []

        injfinder_exe = os.path.basename(cp.get("executables", "injfinder"))
        injfinder_class = select_generic_executable(workflow, "injfinder")
        injfinder_jobs = injfinder_class(cp, "injfinder", ifo=ifos,
                                         out_dir=output_dir, tags=tags)

        injcombiner_exe = os.path.basename(cp.get("executables",
                                                  "injcombiner"))
        injcombiner_class = select_generic_executable(workflow, "injcombiner")
        injcombiner_jobs = injcombiner_class(cp, "injcombiner", ifo=ifos,
                                             out_dir=output_dir, tags=tags)

        injfinder_outs = FileList([])
        for inj_tag in inj_tags:
            triggers = FileList([file for file in inj_trig_files \
                                 if inj_tag in file.tag_str])
            injections = FileList([file for file in inj_files \
                                   if inj_tag in file.tag_str])
            trig_cache = [file for file in inj_trig_caches \
                          if inj_tag in file.tag_str][0]
            inj_cache = [file for file in inj_caches \
                         if inj_tag in file.tag_str][0]
            injfinder_node = injfinder_jobs.create_node(triggers, injections,
                    trig_cache, inj_cache, segment_dir, tags=tags)
            injfinder_nodes.append(injfinder_node)
            pp_nodes.append(injfinder_node)
            workflow.add_node(injfinder_node)

            # Add post-processing files as File objects and make cache
            name_string = injections[0].tagged_description.rsplit('_', 1)[0]
            found_file_name = name_string + "_FOUND"
            found_file = File(ifos, found_file_name, full_segment,
                              extension="xml", directory=output_dir, tags=[])
            missed_file_name = name_string + "_MISSED"
            missed_file = File(ifos, missed_file_name, full_segment,
                               extension="xml", directory=output_dir, tags=[])

            injfinder_outs.extend(FileList([found_file, missed_file]))

        fm_cache = File(ifos, "foundmissed", full_segment,
                        extension="lcf", directory=output_dir)
        fm_cache.PFN(fm_cache.cache_entry.path, site="local")
        fP = open(fm_cache.storage_path, "w")
        for entry in injfinder_outs:
            start = str(int(full_segment[0]))
            duration = str(int(abs(full_segment)))
            print >> fP, "%s %s %s %s file://localhost%s" \
                %(ifos, entry.description.upper(), start, duration,
                  entry.storage_path)
        fP.close()

        injcombiner_tags = [inj_tag for inj_tag in inj_tags \
                            if "DETECTION" not in inj_tag]
        for injcombiner_tag in injcombiner_tags:
            max_inc = cp.get_opt_tags("injections", "max-inc",
                                      [injcombiner_tag])
            inj_str = injcombiner_tag[:4]
            injcombiner_node = injcombiner_jobs.create_node(fm_cache, inj_str,
                    max_inc, workflow.analysis_time)
            pp_nodes.append(injcombiner_node)
            workflow.add_node(injcombiner_node)
            for injfinder_node in injfinder_nodes:
                dep = dax.Dependency(parent=injfinder_node._dax_node,
                                     child=injcombiner_node._dax_node)
                workflow._adag.addDependency(dep)

    # Initialise sbv_plotter class
    sbv_plotter_outs = FileList([])
    sbv_plotter_jobs = sbv_plotter_class(cp, "sbv_plotter", ifo=ifos,
                                         out_dir=output_dir, tags=tags)

    # Initialise efficiency class
    efficiency_outs = FileList([])
    efficiency_jobs = efficiency_class(cp, "efficiency", ifo=ifos,
                                       out_dir=output_dir, tags=tags)

    # Initialise html_summary class
    html_summary_jobs = html_summary_class(cp, "html_summary", ifo=ifos,
                                           out_dir=output_dir, tags=tags)

    # Add trig_cluster jobs and their corresponding plotting jobs
    trig_name = cp.get("workflow", "trigger-name")
    cluster_out_tags = ["ALL_TIMES", "OFFSOURCE", "ONSOURCE"]

    for out_tag in cluster_out_tags:
        trig_cluster_node, \
        unclustered_trigs = trig_cluster_jobs.create_node(trig_files[-1],
                tags=["GRB%s_%s" % (trig_name, out_tag)])
        trig_cluster_outs.extend(trig_cluster_node.output_files)

        if out_tag != "ONSOURCE":
            # Add memory requirememnt for jobs with potentially large files
            trig_cluster_node.set_memory(1300)
            pp_nodes.append(trig_cluster_node)
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
            pp_nodes.append(sbv_plotter_node)
            workflow.add_node(sbv_plotter_node)

            if out_tag == "OFFSOURCE":
                offsource_clustered = parent
                off_node = sbv_plotter_node

            # Also add sbv_plotter job for unclustered triggers
            parent = unclustered_trigs
            sbv_out_tags = [out_tag, "_unclustered"]
            sbv_plotter_node = sbv_plotter_jobs.create_node(parent,
                                                            segment_dir,
                                                            tags=sbv_out_tags)
            sbv_plotter_node.set_memory(1300)
            pp_nodes.append(sbv_plotter_node)
            workflow.add_node(sbv_plotter_node)
            dep = dax.Dependency(parent=trig_combiner_node._dax_node,
                                 child=sbv_plotter_node._dax_node)
            workflow._adag.addDependency(dep)
        else:
            pp_nodes.append(trig_cluster_node)
            workflow.add_node(trig_cluster_node)
            dep = dax.Dependency(parent=trig_combiner_node._dax_node,
                                 child=trig_cluster_node._dax_node)
            workflow._adag.addDependency(dep)

            # Add efficiency job for on/off
            parent = trig_cluster_node.output_files[0]
            efficiency_tag = [out_tag]
            efficiency_node = efficiency_jobs.create_node(parent,
                                                          offsource_clustered,
                                                          segment_dir,
                                                          tags=efficiency_tag)
            pp_nodes.append(efficiency_node)
            workflow.add_node(efficiency_node)
            dep = dax.Dependency(parent=off_node._dax_node,
                                 child=efficiency_node._dax_node)
            workflow._adag.addDependency(dep)

    # Add further trig_cluster jobs for trials
    trials = int(cp.get("trig_combiner", "num-trials"))
    trial = 1

    while trial <= trials:
        trig_cluster_node, \
        unclustered_trigs = trig_cluster_jobs.create_node(trig_files[-1],
                tags=["GRB%s_OFFTRIAL_%d" % (trig_name, trial)])
        trig_cluster_outs.extend(trig_cluster_node.output_files)
        pp_nodes.append(trig_cluster_node)
        workflow.add_node(trig_cluster_node)
        dep = dax.Dependency(parent=trig_combiner_node._dax_node,
                             child=trig_cluster_node._dax_node)
        workflow._adag.addDependency(dep)

        # Add efficiency job
        parent = trig_cluster_node.output_files[0]
        efficiency_tag = ["OFFTRIAL_%d" % trial]
        efficiency_node = efficiency_jobs.create_node(parent,
                                                      offsource_clustered,
                                                      segment_dir,
                                                      tags=efficiency_tag)
        pp_nodes.append(efficiency_node)
        workflow.add_node(efficiency_node)
        dep = dax.Dependency(parent=off_node._dax_node,
                             child=efficiency_node._dax_node)
        workflow._adag.addDependency(dep)

        trial += 1

    # Initialise html_summary class and set up job
    #FIXME: We may want this job to run even if some jobs fail
    html_summary_jobs = html_summary_class(cp, "html_summary", ifo=ifos,
                                           out_dir=output_dir, tags=tags)
    html_summary_node = html_summary_jobs.create_node(config_file=config_file)
    workflow.add_node(html_summary_node)
    for pp_node in pp_nodes:
        dep = dax.Dependency(parent=pp_node._dax_node,
                             child=html_summary_node._dax_node)
        workflow._adag.addDependency(dep)

    post_proc_outs.extend(trig_cluster_outs)

    return post_proc_outs
