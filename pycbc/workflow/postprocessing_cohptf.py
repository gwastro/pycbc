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

def setup_coh_PTF_post_processing(workflow, trigger_files, trigger_cache, 
        output_dir, segment_dir, injection_trigger_files=None,
        injection_files=None, injection_trigger_caches=None,
        injection_caches=None, config_file=None, run_dir=None, ifos=None,
        web_dir=None, inj_tags=[], tags=[], **kwargs):
    """
    This function aims to be the gateway for running postprocessing in CBC
    offline workflows. Post-processing generally consists of calculating the
    significance of triggers and making any statements about trigger rates.
    Dedicated plotting jobs do not belong here.

    Parameters
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
                trigger_files, trigger_cache, injection_trigger_files,
                injection_files, injection_trigger_caches, injection_caches,
                config_file, output_dir, web_dir, segment_dir, ifos=ifos,
                inj_tags=inj_tags, tags=tags, **kwargs)
    else:
        errMsg = "Post-processing method not recognized. Must be "
        errMsg += "COH_PTF_WORKFLOW."
        raise ValueError(errMsg)

    logging.info("Leaving post-processing module.")

    return post_proc_files


def setup_postproc_coh_PTF_workflow(workflow, trig_files, trig_cache,
                                    inj_trig_files, inj_files, inj_trig_caches,
                                    inj_caches, config_file, output_dir,
                                    html_dir, segment_dir, ifos, inj_tags=[],
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
    trig_files : pycbc.workflow.core.FileList
        A FileList containing the combined databases.
   
    Returns
    --------
    
    """
    cp = workflow.cp
    full_segment = trig_files[0].segment
    trig_name = cp.get("workflow", "trigger-name")
    grb_string = "GRB" + trig_name
    num_trials = int(cp.get("trig_combiner", "num-trials"))

    pp_outs = FileList([])
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
    trig_combiner_out_tags = ["OFFSOURCE", "ONSOURCE", "ALL_TIMES"]
    trig_combiner_jobs = trig_combiner_class(cp, "trig_combiner", ifo=ifos, 
                                             out_dir=output_dir, tags=tags)
    trig_combiner_node, trig_combiner_outs = trig_combiner_jobs.create_node(\
            trig_files, segment_dir, out_tags=trig_combiner_out_tags,
            tags=tags)
    pp_nodes.append(trig_combiner_node)
    workflow.add_node(trig_combiner_node)
    pp_outs.extend(trig_combiner_outs)

    # Initialise trig_cluster class
    trig_cluster_outs = FileList([])
    trig_cluster_jobs = trig_cluster_class(cp, "trig_cluster", ifo=ifos,
                                           out_dir=output_dir, tags=tags)

    # Set up injfinder jobs
    if cp.has_section("workflow-injections"):
        injfinder_nodes = []
        injcombiner_parent_nodes = []
        inj_sbv_plotter_parent_nodes = []

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
            injfinder_node, curr_outs = injfinder_jobs.create_node(\
                    triggers, injections, segment_dir, tags=[inj_tag])
            injfinder_nodes.append(injfinder_node)
            pp_nodes.append(injfinder_node)
            workflow.add_node(injfinder_node)
            injfinder_outs.extend(curr_outs)
            if "DETECTION" not in curr_outs[0].tagged_description:
                injcombiner_parent_nodes.append(injfinder_node)
            else:
                inj_sbv_plotter_parent_nodes.append(injfinder_node)

        pp_outs.extend(injfinder_outs)

        # Make injfinder output cache
        fm_cache = File(ifos, "foundmissed", full_segment,
                        extension="lcf", directory=output_dir)
        fm_cache.PFN(fm_cache.cache_entry.path, site="local")
        injfinder_outs.convert_to_lal_cache().tofile(\
                open(fm_cache.storage_path, "w"))
        pp_outs.extend(FileList([fm_cache]))

        # Set up injcombiner jobs
        injcombiner_outs = FileList([file for file in injfinder_outs \
                                     if "DETECTION" in file.tag_str])
        injcombiner_tags = [inj_tag for inj_tag in inj_tags \
                            if "DETECTION" not in inj_tag]
        injcombiner_out_tags = [injcombiner_outs[0].tag_str.rsplit('_', 1)[0]]
        injcombiner_nodes = []

        for injcombiner_tag in injcombiner_tags:
            max_inc = cp.get_opt_tags("injections", "max-inc",
                                      [injcombiner_tag])
            inj_str = injcombiner_tag[:4]
            inputs = FileList([file for file in injfinder_outs \
                               if injcombiner_tag in file.tagged_description])
            #                   if any(tag in file.tagged_description \
            #                          for tag in injcombiner_tags)])
            injcombiner_node, curr_outs = injcombiner_jobs.create_node(\
                    fm_cache, inputs, inj_str, max_inc, workflow.analysis_time)
            injcombiner_nodes.append(injcombiner_node)
            injcombiner_out_tags.append("%s_FILTERED_%s" % (inj_str, max_inc))
            injcombiner_outs.extend(curr_outs)
            pp_outs.extend(curr_outs)
            pp_nodes.append(injcombiner_node)
            workflow.add_node(injcombiner_node)
            for parent_node in injcombiner_parent_nodes:
                dep = dax.Dependency(parent=parent_node._dax_node,
                                     child=injcombiner_node._dax_node)
                workflow._adag.addDependency(dep)

        # Initialise injection_efficiency class
        inj_efficiency_jobs = efficiency_class(cp, "inj_efficiency", ifo=ifos,
                                               out_dir=output_dir, tags=tags)

    # Initialise sbv_plotter class
    sbv_plotter_outs = FileList([])
    sbv_plotter_jobs = sbv_plotter_class(cp, "sbv_plotter", ifo=ifos,
                                         out_dir=output_dir, tags=tags)

    # Initialise efficiency class
    efficiency_outs = FileList([])
    efficiency_jobs = efficiency_class(cp, "efficiency", ifo=ifos,
                                       out_dir=output_dir, tags=tags)

    # Add trig_cluster jobs and their corresponding plotting jobs
    for out_tag in trig_combiner_out_tags:
        unclust_file = [file for file in trig_combiner_outs \
                        if out_tag in file.tag_str][0]
        trig_cluster_node, curr_outs = trig_cluster_jobs.create_node(\
                unclust_file)
        trig_cluster_outs.extend(curr_outs)
        clust_file = curr_outs[0]
        if out_tag != "ONSOURCE":
            # Add memory requirememnt for jobs with potentially large files
            trig_cluster_node.set_memory(1300)
            pp_nodes.append(trig_cluster_node)
            workflow.add_node(trig_cluster_node)
            dep = dax.Dependency(parent=trig_combiner_node._dax_node,
                                 child=trig_cluster_node._dax_node)
            workflow._adag.addDependency(dep)

            # Add sbv_plotter job
            sbv_out_tags = [out_tag, "_clustered"]
            sbv_plotter_node = sbv_plotter_jobs.create_node(clust_file,
                                                            segment_dir,
                                                            tags=sbv_out_tags)
            pp_nodes.append(sbv_plotter_node)
            workflow.add_node(sbv_plotter_node)
            dep = dax.Dependency(parent=trig_cluster_node._dax_node,
                                 child=sbv_plotter_node._dax_node)
            workflow._adag.addDependency(dep)

            # Add injection sbv_plotter nodes if appropriate
            if out_tag == "OFFSOURCE" and \
                    cp.has_section("workflow-injections"):
                offsource_clustered = clust_file
                off_node = sbv_plotter_node

                found_inj_files = FileList([file for file in injcombiner_outs \
                                            if "FOUND" in file.tag_str])
                for curr_injs in found_inj_files:
                    curr_tags = [tag for tag in injcombiner_out_tags \
                                 if tag in curr_injs.name]
                    curr_tags.append("_clustered")
                    sbv_plotter_node = sbv_plotter_jobs.create_node(clust_file,
                            segment_dir, inj_file=curr_injs, tags=curr_tags)
                    pp_nodes.append(sbv_plotter_node)
                    workflow.add_node(sbv_plotter_node)
                    dep = dax.Dependency(parent=trig_cluster_node._dax_node,
                                         child=sbv_plotter_node._dax_node)
                    workflow._adag.addDependency(dep)
                    if "DETECTION" in curr_injs.tagged_description:
                        for parent_node in inj_sbv_plotter_parent_nodes:
                            dep = dax.Dependency(parent=parent_node._dax_node,
                                    child=sbv_plotter_node._dax_node)
                            workflow._adag.addDependency(dep)
                    else:
                        for parent_node in injcombiner_nodes:
                            dep = dax.Dependency(parent=parent_node._dax_node,
                                    child=sbv_plotter_node._dax_node)
                            workflow._adag.addDependency(dep)

            # Also add sbv_plotter job for unclustered triggers
            sbv_plotter_node = sbv_plotter_jobs.create_node(unclust_file,
                    segment_dir, tags=[out_tag, "_unclustered"])
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
            efficiency_node = efficiency_jobs.create_node(clust_file,
                    offsource_clustered, segment_dir, tags=[out_tag])
            pp_nodes.append(efficiency_node)
            workflow.add_node(efficiency_node)
            dep = dax.Dependency(parent=off_node._dax_node,
                                 child=efficiency_node._dax_node)
            workflow._adag.addDependency(dep)

            if cp.has_section("workflow-injections"):
                for tag in injcombiner_out_tags:
                    if "_FILTERED_" in tag:
                        inj_set_tag = [t for t in inj_tags if \
                                       str(tag).replace("_FILTERED_", "") \
                                       in t][0]
                    else:
                        inj_set_tag = str(tag)
                    
                    found_file = [file for file in injcombiner_outs \
                                  if tag + "_FOUND" in file.tag_str][0]
                    missed_file = [file for file in injcombiner_outs \
                                   if tag + "_MISSED" in file.tag_str][0]
                    inj_efficiency_node = inj_efficiency_jobs.create_node(\
                            clust_file, offsource_clustered, segment_dir,
                            found_file, missed_file, tags=[out_tag, tag,
                                                           inj_set_tag])
                    pp_nodes.append(inj_efficiency_node)
                    workflow.add_node(inj_efficiency_node)
                    dep = dax.Dependency(parent=off_node._dax_node,
                                         child=inj_efficiency_node._dax_node)
                    workflow._adag.addDependency(dep)
                    for injcombiner_node in injcombiner_nodes:
                        dep = dax.Dependency(parent=injcombiner_node._dax_node,
                                child=inj_efficiency_node._dax_node)
                        workflow._adag.addDependency(dep)
                    for injfinder_node in injfinder_nodes:
                        dep = dax.Dependency(parent=injfinder_node._dax_node,
                                child=inj_efficiency_node._dax_node)
                        workflow._adag.addDependency(dep)

    # Add further trig_cluster jobs for trials
    trial = 1

    while trial <= num_trials:
        trial_tag = "OFFTRIAL_%d" % trial
        unclust_file = [file for file in trig_combiner_outs \
                        if trial_tag in file.tag_str][0]
        trig_cluster_node, clust_outs = trig_cluster_jobs.create_node(\
                unclust_file)
        clust_file = clust_outs[0]
        trig_cluster_outs.extend(clust_outs)
        pp_nodes.append(trig_cluster_node)
        workflow.add_node(trig_cluster_node)
        dep = dax.Dependency(parent=trig_combiner_node._dax_node,
                             child=trig_cluster_node._dax_node)
        workflow._adag.addDependency(dep)

        # Add efficiency job
        efficiency_node = efficiency_jobs.create_node(clust_file,
                offsource_clustered, segment_dir, tags=[trial_tag])
        pp_nodes.append(efficiency_node)
        workflow.add_node(efficiency_node)
        dep = dax.Dependency(parent=off_node._dax_node,
                             child=efficiency_node._dax_node)
        workflow._adag.addDependency(dep)
        dep = dax.Dependency(parent=trig_cluster_node._dax_node,
                             child=efficiency_node._dax_node)
        workflow._adag.addDependency(dep)

        # Adding inj_efficiency job
        if cp.has_section("workflow-injections"):
            for tag in injcombiner_out_tags:
                if "_FILTERED_" in tag:
                    inj_set_tag = [t for t in inj_tags if \
                                   str(tag).replace("_FILTERED_", "") in t][0]
                else:
                    inj_set_tag = str(tag)

                found_file = [file for file in injcombiner_outs \
                              if tag + "_FOUND" in file.tag_str][0]
                missed_file = [file for file in injcombiner_outs \
                               if tag + "_MISSED" in file.tag_str][0]
                inj_efficiency_node = inj_efficiency_jobs.create_node(\
                        clust_file, offsource_clustered, segment_dir,
                        found_file, missed_file, tags=[trial_tag, tag,
                                                       inj_set_tag])
                pp_nodes.append(inj_efficiency_node)
                workflow.add_node(inj_efficiency_node)
                dep = dax.Dependency(parent=off_node._dax_node,
                                     child=inj_efficiency_node._dax_node)
                workflow._adag.addDependency(dep)
                for injcombiner_node in injcombiner_nodes:
                    dep = dax.Dependency(parent=injcombiner_node._dax_node,
                                         child=inj_efficiency_node._dax_node)
                    workflow._adag.addDependency(dep)
                for injfinder_node in injfinder_nodes:
                    dep = dax.Dependency(parent=injfinder_node._dax_node,
                                         child=inj_efficiency_node._dax_node)
                    workflow._adag.addDependency(dep)

        trial += 1

    # Initialise html_summary class and set up job
    #FIXME: We may want this job to run even if some jobs fail
    html_summary_jobs = html_summary_class(cp, "html_summary", ifo=ifos,
                                           out_dir=output_dir, tags=tags)
    if cp.has_section("workflow-injections"):
        tuning_tags = [inj_tag for inj_tag in injcombiner_out_tags \
                       if "DETECTION" in inj_tag]
        exclusion_tags = [inj_tag for inj_tag in injcombiner_out_tags \
                          if "DETECTION" not in inj_tag]
        html_summary_node = html_summary_jobs.create_node(c_file=config_file,
                tuning_tags=tuning_tags, exclusion_tags=exclusion_tags,
                html_dir=html_dir)
    else:
        html_summary_node = html_summary_jobs.create_node(c_file=config_file,
                                                          html_dir=html_dir)
    workflow.add_node(html_summary_node)
    for pp_node in pp_nodes:
        dep = dax.Dependency(parent=pp_node._dax_node,
                             child=html_summary_node._dax_node)
        workflow._adag.addDependency(dep)

    # Make the open box shell script
    open_box_cmd = html_summary_node.executable.get_pfn() + " "
    open_box_cmd += ' '.join(html_summary_node._args + \
                             html_summary_node._options)
    open_box_cmd += " --open-box"
    open_box_path = "%s/open_the_box.sh" % output_dir
    f = open(open_box_path, "w")
    f.write("#!/bin/sh\n%s" % open_box_cmd)
    f.close()
    os.chmod(open_box_path, 0500)

    pp_outs.extend(trig_cluster_outs)

    return pp_outs
