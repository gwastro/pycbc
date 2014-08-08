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
This module is used to prepare output from the coincidence and/or matched filter
stages of the workflow for post-processing (post-processing = calculating
significance of candidates and making any rates statements). For details of this
module and its capabilities see here:
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/NOTYETCREATED.html
"""

from __future__ import division

import os
import os.path
import logging
from glue import segments
from pycbc.workflow.core import FileList, make_analysis_dir
from pycbc.workflow.jobsetup import select_generic_executable

def setup_postprocessing_preparation(workflow, triggerFiles, output_dir,
                                     tags=[], **kwargs):
    """
    This function aims to be the gateway for preparing the output of the
    coincidence and/or matched-filtering stages of the workflow for calculation 
    of the significance of triggers and any rate statements that are to made. In
    practice this normally means combining output files, performing any
    clustering and performing mapping between triggers and simulations where
    needed.

    Properties
    -----------
    workflow : pycbc.workflow.core.Workflow
        The Workflow instance that the coincidence jobs will be added to.
    triggerFiles : pycbc.workflow.core.FileList
        An FileList of the trigger files that are used as
        input at this stage.
    output_dir : path
        The directory in which output files will be stored.
    tags : list of strings (optional, default = [])
        A list of the tagging strings that will be used for all jobs created
        by this call to the workflow. An example might be ['POSTPROC1'] or
        ['DENTYSNEWPOSTPROC']. This will be used in output names.

    Returns
    --------
    postProcPreppedFiles : pycbc.workflow.core.FileList
        A list of files that can be used as input for the post-processing stage.
    """
    logging.info("Entering post-processing preparation stage.")
    make_analysis_dir(output_dir)

    # Parse for options in .ini file
    postProcPrepMethod = workflow.cp.get_opt_tags("workflow-postprocprep",
                                        "postprocprep-method", tags)

    # Scope here for adding different options/methods here. For now we only
    # have the single_stage ihope method which consists of converting the
    # ligolw_thinca output xml into one file, clustering, performing injection
    # finding and putting everything into one SQL database.
    if postProcPrepMethod == "PIPEDOWN_WORKFLOW":
        # If you want the intermediate output files, call this directly
        postPostPreppedFiles,_,_,_ = setup_postprocprep_pipedown_workflow(\
                                       workflow, triggerFiles, output_dir,
                                       tags=tags, **kwargs) 
    elif postProcPrepMethod == "GSTLAL_POSTPROCPREP":
        postPostPreppedFiles = setup_postprocprep_gstlal_workflow(workflow,
                                 triggerFiles, output_dir, tags=tags, **kwargs)
    else:
        errMsg = "Post-processing preparation method not recognized. Must be "
        errMsg += "one of PIPEDOWN_AHOPE or GSTLAL_POSTPROCPREP."
        raise ValueError(errMsg)

    logging.info("Leaving post-processing preparation module.")

    return postPostPreppedFiles

def setup_postprocprep_pipedown_workflow(workflow, coincFiles, output_dir,
                                      tags=[], injectionFiles=None,
                                      vetoFiles=None, injLessTag=None,
                                      injectionTags=[], vetoCats=[]):
    """
    Properties
    -----------
    workflow : pycbc.workflow.core.Workflow
        The Workflow instance that the coincidence jobs will be added to.
    coincFiles : pycbc.workflow.core.FileList
        An FileList of the coincident trigger files that are used as
        input at this stage.
    output_dir : path
        The directory in which output files will be stored.
    tags : list of strings (optional, default = [])
        A list of the tagging strings that will be used for all jobs created
        by this call to the workflow. An example might be ['POSTPROC1'] or
        ['DENTYSNEWPOSTPROC']. This will be used in output names.
    injectionFiles : pycbc.workflow.core.FileList (optional, default=None)
        The injection files to be used in this stage. An empty list (or any
        other input that evaluates as false) is valid and will imply that no
        injections are being done.
    vetoFiles : pycbc.workflow.core.FileList (required)
        The data quality files to be used in this stage. This is required and
        will be used to determine the analysed times when doing post-processing.
    injLessTag : string (required)
        The tag that identifies files that do not have simulations in them.
        Ie. the primary search results.
    injectionTags : list of strings (optional, default = [])
        Each injection file has a unique tag. If used in the method, this
        tells the post-processing preparation code which injection tags it
        should include when creating the combined output.
    vetoCats : list of integers (optional, default = [])
        Decide which set of veto files should be used in the post-processing
        preparation. This is used, for example, to tell the workflow that you
        are only interested in quoting results at categories 2, 3 and 4. In
        which case just supply [2,3,4] for those veto files here.

    Returns
    --------
    finalFiles : pycbc.workflow.core.FileList
        A list of the single SQL database storing the clustered, injection
        found, triggers for all injections, time slid and zero lag analyses.
    initialSqlFiles : pycbc.workflow.core.FileList
        The SQL files before clustering is applied and injection finding
        performed.
    clusteredSqlFiles : pycbc.workflow.core.FileList
        The clustered SQL files before injection finding performed.
    combinedSqlFiles : pycbc.workflow.core.FileList
        A combined file containing all triggers after clustering, including
        the injection and veto tables, but before injection finding performed.
        Probably there is no need to ever keep this file and it will be a
        temporary file in most cases.
    """
    if not vetoCats:
        raise ValueError("Veto cats is required.")

    # Setup needed exe classes
    sqliteCombine1ExeTag = workflow.cp.get_opt_tags("workflow-postprocprep",
                                   "postprocprep-combiner1-exe", tags)
    sqliteCombine1Exe = select_generic_executable(workflow, 
                                                  sqliteCombine1ExeTag)
    sqliteCombine2ExeTag = workflow.cp.get_opt_tags("workflow-postprocprep",
                                   "postprocprep-combiner2-exe", tags)
    sqliteCombine2Exe = select_generic_executable(workflow, 
                                                  sqliteCombine2ExeTag)
    clusterCoincsExeTag = workflow.cp.get_opt_tags("workflow-postprocprep",
                                   "postprocprep-cluster-exe", tags)
    clusterCoincsExe = select_generic_executable(workflow, clusterCoincsExeTag)
    injFindExeTag = workflow.cp.get_opt_tags("workflow-postprocprep",
                                   "postprocprep-injfind-exe", tags)
    injFindExe = select_generic_executable(workflow, injFindExeTag)

    sqliteCombine1Outs = FileList([])
    clusterCoincsOuts = FileList([])
    injFindOuts = FileList([])
    sqliteCombine2Outs = FileList([])
    for vetoCat in vetoCats:
        # FIXME: Some hacking is still needed while we support pipedown
        # FIXME: There are currently 3 names to say cumulative cat_3
        vetoTag = 'CUMULATIVE_CAT_%d' %(vetoCat,)
        dqSegFile = vetoFiles.find_output_with_tag(vetoTag)
        if not len(dqSegFile) == 1:
            errMsg = "Did not find exactly 1 data quality file."
            raise ValueError(errMsg)
        # Don't think this is used here, this is the tag *in* the file
        dqVetoName = 'VETO_CAT%d_CUMULATIVE' %(vetoCat,)
        # FIXME: Here we set the dqVetoName to be compatible with pipedown
        pipedownDQVetoName = 'CAT_%d_VETO' %(vetoCat,)

        sqliteCombine2Inputs = FileList([])
        # Do injection-less jobs first.

        # Combine trig files first
        currTags = tags + [injLessTag, vetoTag]
        trigVetoInpFiles = coincFiles.find_output_with_tag(pipedownDQVetoName)
        trigInpFiles = trigVetoInpFiles.find_output_with_tag(injLessTag)
        trigInpFiles.append(dqSegFile[0])
        sqliteCombine1Job = sqliteCombine1Exe(workflow.cp, 
                                                    sqliteCombine1ExeTag,
                                                    ifo=workflow.ifo_string, 
                                                    out_dir=output_dir,
                                                    tags=currTags)
        sqliteCombine1Node = sqliteCombine1Job.create_node(\
                                          workflow.analysis_time, trigInpFiles)
        workflow.add_node(sqliteCombine1Node)
        # Node has only one output file
        sqliteCombine1Out = sqliteCombine1Node.output_files[0]
        sqliteCombine1Outs.append(sqliteCombine1Out)

        # Cluster coincidences
        clusterCoincsJob = clusterCoincsExe(workflow.cp, 
                                                   clusterCoincsExeTag,
                                                   ifo=workflow.ifo_string, 
                                                   out_dir=output_dir, 
                                                   tags=currTags)
        clusterCoincsNode = clusterCoincsJob.create_node(\
                                     workflow.analysis_time, sqliteCombine1Out)
        workflow.add_node(clusterCoincsNode)
        # Node has only one output file
        clusterCoincsOut = clusterCoincsNode.output_files[0]
        clusterCoincsOuts.append(clusterCoincsOut)
        sqliteCombine2Inputs.append(clusterCoincsOut)

        # Do injection jobs
        for injTag in injectionTags:
            # Combine trig files first
            currTags = tags + [injTag, vetoTag]
            trigInpFiles = trigVetoInpFiles.find_output_with_tag(injTag)
            trigInpFiles.append(dqSegFile[0])
            injFile = injectionFiles.find_output_with_tag(injTag)
            assert (len(injFile) == 1)
            sqliteCombine1Job = sqliteCombine1Exe(workflow.cp, 
                                                  sqliteCombine1ExeTag,
                                                  ifo=workflow.ifo_string, 
                                                  out_dir=output_dir,
                                                  tags=currTags)
            sqliteCombine1Node = sqliteCombine1Job.create_node(\
                                          workflow.analysis_time, trigInpFiles,
                                          injFile=injFile[0], injString=injTag)
            workflow.add_node(sqliteCombine1Node)
            # Node has only one output file
            sqliteCombine1Out = sqliteCombine1Node.output_files[0]
            sqliteCombine1Outs.append(sqliteCombine1Out)

            # Cluster coincidences
            clusterCoincsJob = clusterCoincsExe(workflow.cp, 
                                                clusterCoincsExeTag,
                                                ifo=workflow.ifo_string, 
                                                out_dir=output_dir, 
                                                tags=currTags)
            clusterCoincsNode = clusterCoincsJob.create_node(\
                                     workflow.analysis_time, sqliteCombine1Out)
            workflow.add_node(clusterCoincsNode)
            # Node has only one output file
            clusterCoincsOut = clusterCoincsNode.output_files[0]
            clusterCoincsOuts.append(clusterCoincsOut)
            sqliteCombine2Inputs.append(clusterCoincsOut)

        # Combine everything together and add veto file
        currTags = tags + [vetoTag]
        sqliteCombine2Job = sqliteCombine2Exe(workflow.cp, 
                                              sqliteCombine2ExeTag,
                                              ifo=workflow.ifo_string, 
                                              out_dir=output_dir,
                                              tags=currTags)
        sqliteCombine2Node = sqliteCombine2Job.create_node(\
                                  workflow.analysis_time, sqliteCombine2Inputs)
        workflow.add_node(sqliteCombine2Node)
        sqliteCombine2Out = sqliteCombine2Node.output_files[0]
        sqliteCombine2Outs.append(sqliteCombine2Out)

        # Inj finding
        injFindJob = injFindExe(workflow.cp, injFindExeTag, 
                                          ifo=workflow.ifo_string,
                                          out_dir=output_dir,tags=currTags)
        injFindNode = injFindJob.create_node(workflow.analysis_time,
                                                         sqliteCombine2Out)
        workflow.add_node(injFindNode)
        injFindOut = injFindNode.output_files[0]
        injFindOuts.append(injFindOut)


    return injFindOuts, sqliteCombine1Outs, clusterCoincsOuts,\
           sqliteCombine2Outs


def setup_postprocprep_gstlal_workflow(workflow, coinc_files, output_dir,
                                       tags=[], injection_files=None,
                                       veto_files=None, inj_less_tag=None,
                                       injection_tags=[], veto_cats=[],
                                       summary_xml_files=None):
    """
    Properties
    -----------
    workflow : ahope.Workflow
        The ahope workflow instance that the coincidence jobs will be added to.
    coinc_files : ahope.AhopeFileList
        An AhopeFileList of the coincident trigger files that are used as
        input at this stage.
    output_dir : path
        The directory in which output files will be stored.
    tags : list of strings (optional, default = [])
        A list of the tagging strings that will be used for all jobs created
        by this call to the workflow. An example might be ['POSTPROC1'] or
        ['DENTYSNEWPOSTPROC']. This will be used in output names.
    injection_files : ahope.AhopeFileList (optional, default=None)
        The injection files to be used in this stage. An empty list (or any
        other input that evaluates as false) is valid and will imply that no
        injections are being done.
    veto_files : ahope.AhopeFileList (required)
        The data quality files to be used in this stage. This is required and
        will be used to determine the analysed times when doing post-processing.
    inj_less_tag : string (required)
        The tag that identifies files that do not have simulations in them.
        Ie. the primary search results.
    injection_tags : list of strings (optional, default = [])
        Each injection file has a unique tag. If used in the method, this
        tells the post-processing preparation code which injection tags it
        should include when creating the combined output.
    veto_cats : list of integers (optional, default = [])
        Decide which set of veto files should be used in the post-processing
        preparation. This is used, for example, to tell the workflow that you
        are only interested in quoting results at categories 2, 3 and 4. In
        which case just supply [2,3,4] for those veto files here.
    summary_xml_files : ahope.AhopeFileList
        An AhopeFileList of the output of the analysislogging_utils module.
        Here, this will be one file that includes the segments analysed by the
        workflow.

    Returns
    --------
    finalFiles : ahope.AhopeFileList
        A list of the single SQL database storing the clustered, injection
        found, triggers for all injections, time slid and zero lag analyses.
    initialSqlFiles : ahope.AhopeFileList
        The SQL files before clustering is applied and injection finding
        performed.
    clusteredSqlFiles : ahope.AhopeFileList
        The clustered SQL files before injection finding performed.
    combinedSqlFiles : ahope.AhopeFileList
        A combined file containing all triggers after clustering, including
        the injection and veto tables, but before injection finding performed.
        Probably there is no need to ever keep this file and it will be a
        temporary file in most cases.
    """
    # Sanity checks
    if not len(summary_xml_files) == 1:
        errMsg = "I need exactly one summaryXML file, got %d." \
                                                     %(len(summary_xml_files),)
        raise ValueError(errMsg)

    # Setup needed exe classes
    run_sqlite_exe_name = workflow.cp.get_opt_tags("ahope-postprocprep",
                                   "postprocprep-runsqlite-exe", tags)
    ligolw_sqlite_exe_name = workflow.cp.get_opt_tags("ahope-postprocprep",
                                   "postprocprep-ligolwsqlite-exe", tags) 
    inspinjfind_exe_name = workflow.cp.get_opt_tags("ahope-postprocprep",
                                   "postprocprep-inspinjfind-exe", tags)
    sql_to_xml_exe_name = workflow.cp.get_opt_tags("ahope-postprocprep",
                                   "postprocprep-sqltoxml-exe", tags)
    run_sqlite_exe = select_generic_executable(workflow, run_sqlite_exe_name)
    ligolw_sqlite_exe = select_generic_executable(workflow,
                                                        ligolw_sqlite_exe_name)
    inspinjfind_exe = select_generic_executable(workflow, inspinjfind_exe_name)
    sql_to_xml_exe = select_generic_executable(workflow, sql_to_xml_exe_name)

    # SETUP
    # FIXME: Some hacking is still needed while we support pipedown
    # FIXME: How does gstlal deal with veto categories?
    #         Hardcode to CAT1 for now.
    veto_cat = 1
    veto_tag = 'CUMULATIVE_CAT_%d' %(veto_cat,)
    #dqSegFile = vetoFiles.find_output_with_tag(vetoTag)
    #if not len(dqSegFile) == 1:
    #    errMsg = "Did not find exactly 1 data quality file."
    #    raise ValueError(errMsg)
    # FIXME: Here we set the dqVetoName to be compatible with pipedown
    pipedown_dq_veto_name = 'CAT_%d_VETO' %(veto_cat,)
                                                  
    # STAGE 1: Run ligolw_sqlite on all input files
    # STAGE 2: Run run_sqlite on all outputs of stage 1
    # STAGE 3: Combine all files into one sqlite file
    # STAGE 4: Run run_sqlite on outputs of stage 3
    # STAGE 5: Add segments.xml and inj.xml
    # STAGE 6: Run run_sqlite (cluster an simplify) on outputs of stage 5
    # STAGE 7: Dump SQL database to xml
    # STAGE 8: Run injfind on the xml document
    # STAGE 9: Convert back to SQL

    stage1_outputs = {}
    stage2_outputs = {}
    stage3_outputs = {}
    stage4_outputs = {}
    stage5_outputs = {}
    stage6_outputs = {}
    stage7_outputs = {}
    stage8_outputs = {}
    stage9_outputs = {}
    final_outputs = AhopeFileList([])
    # Do for all injection runs and zero lag
    for inj_tag in [inj_less_tag] + injection_tags:
        curr_tags = tags + [inj_tag, veto_tag]
        trig_veto_inp_files = \
                        coinc_files.find_output_with_tag(pipedown_dq_veto_name)
        trig_inp_files = trig_veto_inp_files.find_output_with_tag(inj_tag)
        stage1_job = ligolw_sqlite_exe(workflow.cp, ligolw_sqlite_exe_name,
                                      ifo=workflow.ifo_string,
                                      out_dir=output_dir,
                                      tags=['STAGE1'] + curr_tags)
        stage2_job = run_sqlite_exe(workflow.cp, run_sqlite_exe_name,
                                      ifo=workflow.ifo_string,
                                      out_dir=output_dir,
                                      tags=['STAGE2'] + curr_tags)
        stage3_job = ligolw_sqlite_exe(workflow.cp, ligolw_sqlite_exe_name,
                                      ifo=workflow.ifo_string,
                                      out_dir=output_dir,
                                      tags=['STAGE3'] + curr_tags)
        stage4_job = run_sqlite_exe(workflow.cp, run_sqlite_exe_name,
                                      ifo=workflow.ifo_string,
                                      out_dir=output_dir,
                                      tags=['STAGE4'] + curr_tags)
        stage5_job = ligolw_sqlite_exe(workflow.cp, ligolw_sqlite_exe_name,
                                      ifo=workflow.ifo_string,
                                      out_dir=output_dir,
                                      tags=['STAGE5'] + curr_tags)
        if inj_tag == inj_less_tag:
            # For zero-lag we stop here, so use the FINAL tag to indicate this
            stage6_zl_job = run_sqlite_exe(workflow.cp, run_sqlite_exe_name,
                                          ifo=workflow.ifo_string,
                                          out_dir=output_dir,
                                          tags=['FINAL'] + curr_tags)
        else:
            stage6_job = run_sqlite_exe(workflow.cp, run_sqlite_exe_name,
                                          ifo=workflow.ifo_string,
                                          out_dir=output_dir,
                                          tags=['STAGE6'] + curr_tags)
            stage7_job = sql_to_xml_exe(workflow.cp, sql_to_xml_exe_name,
                                          ifo=workflow.ifo_string,
                                          out_dir=output_dir,
                                          tags=['STAGE7'] + curr_tags)
            stage8_job = inspinjfind_exe(workflow.cp, inspinjfind_exe_name,
                                          ifo=workflow.ifo_string,
                                          out_dir=output_dir,
                                          tags=['STAGE8'] + curr_tags)
            stage9_job = ligolw_sqlite_exe(workflow.cp, run_sqlite_exe_name,
                                          ifo=workflow.ifo_string,
                                          out_dir=output_dir,
                                          tags=['FINAL'] + curr_tags)

        stage1_outputs[inj_tag] = AhopeFileList([])
        stage2_outputs[inj_tag] = AhopeFileList([])
        for file in trig_inp_files:
            stage1_node = stage1_job.create_node(file.segment, [file])
            workflow.add_node(stage1_node)
            # Node has only one output file
            stage1_out = stage1_node.output_files[0]
            stage1_outputs[inj_tag].append(stage1_out)
            stage2_node = stage2_job.create_node(stage1_out.segment,
                                                                    stage1_out)
            workflow.add_node(stage2_node)
            # Node has only one output file
            stage2_out = stage2_node.output_files[0]
            stage2_outputs[inj_tag].append(stage2_out)

        stage3_node = stage3_job.create_node(workflow.analysis_time,
                                                       stage2_outputs[inj_tag])
        workflow.add_node(stage3_node)
        # Node has only one output file
        stage3_out = stage3_node.output_files[0]
        stage3_outputs[inj_tag] = stage3_out
        stage4_node = stage4_job.create_node(workflow.analysis_time,
                                                                    stage3_out)
        workflow.add_node(stage4_node)
        # Node has only one output file
        stage4_out = stage4_node.output_files[0]
        stage4_outputs[inj_tag] = stage4_out

        stage5_inputs = [stage4_out]
        stage5_inputs.append(summary_xml_files[0])
        if inj_tag != inj_less_tag:
            inj_file = injection_files.find_output_with_tag(inj_tag)
            assert (len(inj_file) == 1)
            stage5_inputs.append(inj_file[0])
        stage5_node = stage5_job.create_node(workflow.analysis_time,
                                                                 stage5_inputs)
        workflow.add_node(stage5_node)
        # Node has only one output file
        stage5_out = stage5_node.output_files[0]
        stage5_outputs[inj_tag] = stage5_out
  
        if inj_tag == inj_less_tag:
            stage6_node = stage6_zl_job.create_node(workflow.analysis_time,
                                                                    stage5_out)
            workflow.add_node(stage6_node)
            stage6_out = stage6_node.output_files[0]
            stage6_outputs[inj_tag] = stage6_out
            final_outputs.append(stage6_out)
        else:
            stage6_node = stage6_job.create_node(workflow.analysis_time,
                                                                    stage5_out)
            workflow.add_node(stage6_node)
            stage6_out = stage6_node.output_files[0]
            stage6_outputs[inj_tag] = stage6_out
            stage7_node = stage7_job.create_node(workflow.analysis_time,
                                                                    stage6_out)
            workflow.add_node(stage7_node)
            stage7_out = stage7_node.output_files[0]
            stage7_outputs[inj_tag] = stage7_out
            stage8_node = stage8_job.create_node(workflow.analysis_time,
                                                                    stage7_out)
            workflow.add_node(stage8_node)
            stage8_out = stage8_node.output_files[0]
            stage8_outputs[inj_tag] = stage8_out
            stage9_node = stage9_job.create_node(workflow.analysis_time,
                                                                  [stage8_out])
            workflow.add_node(stage9_node)
            stage9_out = stage9_node.output_files[0]
            stage9_outputs[inj_tag] = stage9_out
            final_outputs.append(stage9_out)

    # FIXME: Maybe contatenate and return all other outputs if needed elsewhere
    return final_outputs






    
