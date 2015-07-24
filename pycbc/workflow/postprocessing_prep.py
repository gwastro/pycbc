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
significance of candidates, following up interesting events and injections and 
making sensitivity and/or rates statements). For details of this
module and its capabilities see here:
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/workflow/postprocprep.html
"""

import logging
from pycbc.workflow.core import FileList, make_analysis_dir
from pycbc.workflow.jobsetup import select_generic_executable
from pycbc.workflow.core import get_random_label

def setup_postprocessing_preparation(workflow, triggerFiles, output_dir,
                                     tags=[], **kwargs):
    """
    This function aims to be the gateway for preparing the output of the
    coincidence and/or matched-filtering stages of the workflow for calculation 
    of the significance of triggers and any rate statements that are to made. In
    practice this normally means combining output files, performing any
    clustering and performing mapping between triggers and simulations where
    needed.

    Parameters
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
    logging.info("Entering post-processing preparation module.")
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
        postPostPreppedFiles,_,_,_ = setup_postprocprep_pipedown_workflow(
                                       workflow, triggerFiles, output_dir,
                                       tags=tags, **kwargs)
    elif postProcPrepMethod == "PIPEDOWN_REPOP":
        postPostPreppedFiles,_,_,_ = setup_postprocprep_pipedown_workflow(
                                       workflow, triggerFiles, output_dir,
                                       tags=tags, do_repop=True, **kwargs)
    elif postProcPrepMethod == "GSTLAL_POSTPROCPREP":
        postPostPreppedFiles = setup_postprocprep_gstlal_workflow(workflow,
                                 triggerFiles, output_dir, tags=tags, **kwargs)
    else:
        errMsg = "Post-processing preparation method not recognized. Must be "
        errMsg += "one of PIPEDOWN_WORKFLOW or GSTLAL_POSTPROCPREP."
        raise ValueError(errMsg)

    logging.info("Leaving post-processing preparation module.")

    return postPostPreppedFiles

def setup_postprocprep_pipedown_workflow(workflow, coincFiles, output_dir,
                                      tags=[], do_repop=False, 
                                      injectionFiles=None,
                                      vetoFiles=None, injLessTag=None,
                                      injectionTags=[], veto_cats=[]):
    """
    Parameters
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
    do_repop : Boolean
        If False, use the 'coinc_inspiral.snr' column from the coincident 
        trigger files as clustering and ranking statistic; if True, use
        a repop_coinc job before clustering to calculate a different ranking
        statistic and store in the coinc_inspiral table for later use.
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
    veto_cats : list of integers (optional, default = [])
        Decide which set of veto files should be used in the post-processing
        preparation. For example tell the workflow to only generate results
        at cumulative categories 2, 3 and 4 by supplying [2,3,4] here.

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
    if not veto_cats:
        raise ValueError("A non-empty list of veto categories is required.")

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

    if do_repop:
        repopCoincExeTag = workflow.cp.get_opt_tags("workflow-postprocprep",
                                                "postprocprep-repop-exe", tags)
        repopCoincExe = select_generic_executable(workflow, repopCoincExeTag)
        repopCoincOuts = FileList([])

    for cat in veto_cats:
        # FIXME: Some hacking is still needed while we support pipedown
        # FIXME: There are currently 3 names to say cumulative cat_3
        vetoTag = 'CUMULATIVE_CAT_%d' %(cat)
        dqSegFile = vetoFiles.find_output_with_tag(vetoTag)
        if not len(dqSegFile) == 1:
            errMsg = "Did not find exactly 1 data quality file."
            raise ValueError(errMsg)
        # Don't think this is used here, this is the tag *in* the file
        dqVetoName = 'VETO_CAT%d_CUMULATIVE' %(cat)
        # FIXME: Here we set the dqVetoName to be compatible with pipedown
        pipedownDQVetoName = 'CAT_%d_VETO' %(cat)

        sqliteCombine2Inputs = FileList([])
        # Do injection-less jobs first.

        # Choose a label for clustering the jobs
        job_label = get_random_label()

        # Combine trig files first
        currTags = tags + [injLessTag, vetoTag]
        trigVetoInpFiles = coincFiles.find_output_with_tag(pipedownDQVetoName)
        trigInpFiles = trigVetoInpFiles.find_output_with_tag(injLessTag)
        if len(trigInpFiles) == 0:
            err_msg = "No input files found. Workflow would fail."
            raise ValueError(err_msg)
        trigInpFiles.append(dqSegFile[0])
        sqliteCombine1Job = sqliteCombine1Exe(workflow.cp,
                                              sqliteCombine1ExeTag,
                                              ifo=workflow.ifo_string,
                                              out_dir=output_dir,
                                              tags=currTags)
        sqliteCombine1Node = sqliteCombine1Job.create_node(
                                          workflow.analysis_time, trigInpFiles, 
                                          workflow=workflow)
        sqliteCombine1Node.add_profile('pegasus', 'label', job_label)
        workflow.add_node(sqliteCombine1Node)
        # Node has only one output file
        sqliteCombine1Out = sqliteCombine1Node.output_files[0]
        sqliteCombine1Outs.append(sqliteCombine1Out)

        if do_repop:
            repopCoincJob = repopCoincExe(workflow.cp,
                                          repopCoincExeTag,
                                          ifo=workflow.ifo_string,
                                          out_dir=output_dir,
                                          tags=currTags)
            repopCoincNode = repopCoincJob.create_node(workflow.analysis_time,
                                                       sqliteCombine1Out)
            repopCoincNode.add_profile('pegasus', 'label', job_label)
            workflow.add_node(repopCoincNode)
            # Node has only one output file
            repopCoincOut = repopCoincNode.output_files[0]
            repopCoincOuts.append(repopCoincOut)

        # Input file plumbing allowing for possible repop_coinc job
        clusterCoincsIn = repopCoincOut if do_repop else sqliteCombine1Out
        # Cluster coincidences
        clusterCoincsJob = clusterCoincsExe(workflow.cp,
                                            clusterCoincsExeTag,
                                            ifo=workflow.ifo_string, 
                                            out_dir=output_dir, 
                                            tags=currTags)
        clusterCoincsNode = clusterCoincsJob.create_node(
                                       workflow.analysis_time, clusterCoincsIn)
        clusterCoincsNode.add_profile('pegasus', 'label', job_label)
        workflow.add_node(clusterCoincsNode)
        # Node has only one output file
        clusterCoincsOut = clusterCoincsNode.output_files[0]
        clusterCoincsOuts.append(clusterCoincsOut)
        sqliteCombine2Inputs.append(clusterCoincsOut)

        # Do injection jobs
        for injTag in injectionTags:
            # Choose a label for clustering the jobs
            job_label = get_random_label()
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
            sqliteCombine1Node = sqliteCombine1Job.create_node(
                                          workflow.analysis_time, trigInpFiles,
                                          injFile=injFile[0], injString=injTag,
                                          workflow=workflow)
            sqliteCombine1Node.add_profile('pegasus', 'label', job_label)
            workflow.add_node(sqliteCombine1Node)
            # Node has only one output file
            sqliteCombine1Out = sqliteCombine1Node.output_files[0]
            sqliteCombine1Outs.append(sqliteCombine1Out)

            if do_repop:
                repopCoincJob = repopCoincExe(workflow.cp,
                                          repopCoincExeTag,
                                          ifo=workflow.ifo_string,
                                          out_dir=output_dir,
                                          tags=currTags)
                repopCoincNode = repopCoincJob.create_node(
                                     workflow.analysis_time, sqliteCombine1Out)
                repopCoincNode.add_profile('pegasus', 'label', job_label)
                workflow.add_node(repopCoincNode)
                # Node has only one output file
                repopCoincOut = repopCoincNode.output_files[0]
                repopCoincOuts.append(repopCoincOut)

            # Input file plumbing allowing for possible repop_coinc job
            clusterCoincsIn = repopCoincOut if do_repop else sqliteCombine1Out
            # Cluster coincidences
            clusterCoincsJob = clusterCoincsExe(workflow.cp,
                                                clusterCoincsExeTag,
                                                ifo=workflow.ifo_string,
                                                out_dir=output_dir,
                                                tags=currTags)
            clusterCoincsNode = clusterCoincsJob.create_node(
                                       workflow.analysis_time, clusterCoincsIn)
            clusterCoincsNode.add_profile('pegasus', 'label', job_label)
            workflow.add_node(clusterCoincsNode)
            # Node has only one output file
            clusterCoincsOut = clusterCoincsNode.output_files[0]
            clusterCoincsOuts.append(clusterCoincsOut)
            sqliteCombine2Inputs.append(clusterCoincsOut)

        # Choose a new label for pegasus-clustering the jobs
        job_label = get_random_label()

        # Combine everything together and add veto file
        currTags = tags + [vetoTag]
        sqliteCombine2Job = sqliteCombine2Exe(workflow.cp, 
                                              sqliteCombine2ExeTag,
                                              ifo=workflow.ifo_string, 
                                              out_dir=output_dir,
                                              tags=currTags)
        sqliteCombine2Node = sqliteCombine2Job.create_node(
                                  workflow.analysis_time, sqliteCombine2Inputs)
        sqliteCombine2Node.add_profile('pegasus', 'label', job_label)
        workflow.add_node(sqliteCombine2Node)
        sqliteCombine2Out = sqliteCombine2Node.output_files[0]
        sqliteCombine2Outs.append(sqliteCombine2Out)

        # Inj finding
        injFindJob = injFindExe(workflow.cp, injFindExeTag,
                                          ifo=workflow.ifo_string,
                                          out_dir=output_dir,tags=currTags)
        injFindNode = injFindJob.create_node(workflow.analysis_time,
                                                         sqliteCombine2Out)
        injFindNode.add_profile('pegasus', 'label', job_label)
        workflow.add_node(injFindNode)
        injFindOut = injFindNode.output_files[0]
        injFindOuts.append(injFindOut)


    return injFindOuts, sqliteCombine1Outs, clusterCoincsOuts,\
           sqliteCombine2Outs


def setup_postprocprep_gstlal_workflow(workflow, coinc_files, output_dir,
                                       tags=[], injection_files=None,
                                       veto_files=None, inj_less_tag=None,
                                       injection_tags=[], veto_cat=None,
                                       summary_xml_files=None,
                                       likelihood_files=[]):
    """
    Parameters
    -----------
    workflow : workflow.Workflow
        The workflow instance that the coincidence jobs will be added to.
    coinc_files : workflow.FileList
        An FileList of the coincident trigger files that are used as
        input at this stage.
    output_dir : path
        The directory in which output files will be stored.
    tags : list of strings (optional, default = [])
        A list of the tagging strings that will be used for all jobs created
        by this call to the workflow. An example might be ['POSTPROC1'] or
        ['DENTYSNEWPOSTPROC']. This will be used in output names.
    injection_files : workflow.FileList (optional, default=None)
        The injection files to be used in this stage. An empty list (or any
        other input that evaluates as false) is valid and will imply that no
        injections are being done.
    veto_files : workflow.FileList (required)
        The data quality files to be used in this stage. This is required and
        will be used to determine the analysed times when doing post-processing.
    inj_less_tag : string (required)
        The tag that identifies files that do not have simulations in them.
        Ie. the primary search results.
    injection_tags : list of strings (optional, default = [])
        Each injection file has a unique tag. If used in the method, this
        tells the post-processing preparation code which injection tags it
        should include when creating the combined output.
    veto_cat : int (optional, default = None)
        FIXME: How does gstlal deal with veto categories?
        Hardcode to CAT1 for now.
    summary_xml_files : workflow.FileList
        An FileList of the output of the analysislogging_utils module.
        Here, this will be one file that includes the segments analysed by the
        workflow.

    Returns
    --------
    finalFiles : workflow.FileList
        A list of the single SQL database storing the clustered, injection
        found, triggers for all injections, time slid and zero lag analyses.
    initialSqlFiles : workflow.FileList
        The SQL files before clustering is applied and injection finding
        performed.
    clusteredSqlFiles : workflow.FileList
        The clustered SQL files before injection finding performed.
    combinedSqlFiles : workflow.FileList
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
    run_sqlite_exe_name = workflow.cp.get_opt_tags("workflow-postprocprep",
                                   "postprocprep-runsqlite-exe", tags)
    ligolw_sqlite_exe_name = workflow.cp.get_opt_tags("workflow-postprocprep",
                                   "postprocprep-ligolwsqlite-exe", tags) 
    inspinjfind_exe_name = workflow.cp.get_opt_tags("workflow-postprocprep",
                                   "postprocprep-inspinjfind-exe", tags)
    sql_to_xml_exe_name = workflow.cp.get_opt_tags("workflow-postprocprep",
                                   "postprocprep-sqltoxml-exe", tags)
    pycbc_picklehor_exe_name = workflow.cp.get_opt_tags("workflow-postprocprep",
                                   "postprocprep-picklehor-exe", tags)
    pycbc_combllhood_exe_name=workflow.cp.get_opt_tags("workflow-postprocprep",
                                   "postprocprep-combllhood-exe", tags)
    pycbc_genranking_exe_name=workflow.cp.get_opt_tags("workflow-postprocprep",
                                   "postprocprep-genranking-exe", tags)
    pycbc_compllhood_exe_name=workflow.cp.get_opt_tags("workflow-postprocprep",
                                   "postprocprep-compllhood-exe", tags)
    marg_likelihood_exe_name = workflow.cp.get_opt_tags("workflow-postprocprep",
                                   "postprocprep-marglikelihood-exe", tags)
    far_gstlal_exe_name = workflow.cp.get_opt_tags("workflow-postprocprep",
                                   "postprocprep-fargstlal-exe", tags)
    plot_summary_exe_name = workflow.cp.get_opt_tags("workflow-postprocprep",
                                   "postprocprep-plotsummary-exe", tags)
    plot_sensitivity_exe_name=workflow.cp.get_opt_tags("workflow-postprocprep",
                                   "postprocprep-plotsensitivity-exe", tags)
    plot_background_exe_name = workflow.cp.get_opt_tags("workflow-postprocprep",
                                   "postprocprep-plotbackground-exe", tags)
    summary_page_exe_name = workflow.cp.get_opt_tags("workflow-postprocprep",
                                   "postprocprep-summarypage-exe", tags)


    run_sqlite_exe = select_generic_executable(workflow, run_sqlite_exe_name)
    ligolw_sqlite_exe = select_generic_executable(workflow,
                                                        ligolw_sqlite_exe_name)
    inspinjfind_exe = select_generic_executable(workflow, inspinjfind_exe_name)
    sql_to_xml_exe = select_generic_executable(workflow, sql_to_xml_exe_name)
    pycbc_picklehor_exe = select_generic_executable(workflow,
                                                      pycbc_picklehor_exe_name)
    pycbc_combllhood_exe = select_generic_executable(workflow,
                                                     pycbc_combllhood_exe_name)
    pycbc_genranking_exe = select_generic_executable(workflow,
                                                     pycbc_genranking_exe_name)
    pycbc_compllhood_exe = select_generic_executable(workflow,
                                                     pycbc_compllhood_exe_name)
    marg_likelihood_exe = select_generic_executable(workflow,
                                                      marg_likelihood_exe_name)
    far_gstlal_exe = select_generic_executable(workflow, far_gstlal_exe_name)
    plot_summary_exe = select_generic_executable(workflow,
                                                         plot_summary_exe_name)
    plot_sensitivity_exe = select_generic_executable(workflow,
                                                     plot_sensitivity_exe_name)
    plot_background_exe = select_generic_executable(workflow,
                                                      plot_background_exe_name)
    summary_page_exe = select_generic_executable(workflow,
                                                         summary_page_exe_name)


    # SETUP
    # FIXME: Some hacking is still needed while we support pipedown
    # FIXME: How does gstlal deal with veto categories?
    #         Hardcode to CAT1 for now.
    veto_tag = 'CUMULATIVE_CAT_%d' %(veto_cat,)
    dq_seg_file = veto_files.find_output_with_tag(veto_tag)
    assert len(dq_seg_file) == 1
    dq_seg_file = dq_seg_file[0]
    #if not len(dqSegFile) == 1:
    #    errMsg = "Did not find exactly 1 data quality file."
    #    raise ValueError(errMsg)
    # FIXME: Here we set the dqVetoName to be compatible with pipedown
    pipedown_dq_veto_name = 'CAT_%d_VETO' %(veto_cat,)

    # First we need to covert to SQL, this is STAGE0
    # Do for all injection runs and zero lag
    stage0_outputs = {}
    for inj_tag in [inj_less_tag] + injection_tags:
        curr_tags = tags + [inj_tag, veto_tag]
        trig_veto_inp_files = \
                  coinc_files.find_output_with_tag(pipedown_dq_veto_name)
        trig_inp_files = trig_veto_inp_files.find_output_with_tag(inj_tag)
        stage0_job = ligolw_sqlite_exe(workflow.cp, ligolw_sqlite_exe_name,
                                      ifo=workflow.ifo_string,
                                      out_dir=output_dir,
                                      tags=['STAGE0'] + curr_tags)
        stage0_outputs[inj_tag] = FileList([])
        assert len(trig_inp_files) > 0
        for file in trig_inp_files:
            stage0_node = stage0_job.create_node(file.segment, [file])
            workflow.add_node(stage0_node)
            # Node has only one output file
            stage0_out = stage0_node.output_files[0]
            stage0_outputs[inj_tag].append(stage0_out)

    curr_tags = tags + [veto_tag]

    # NOW WE DO LIKELIHOOD SETUP
    pycbc_picklehor_job = pycbc_picklehor_exe(workflow.cp,
                                  pycbc_picklehor_exe_name,
                                  ifo=workflow.ifo_string,
                                  out_dir=output_dir,
                                  tags=curr_tags)
    pycbc_combllhood_job = pycbc_combllhood_exe(workflow.cp,
                                  pycbc_combllhood_exe_name,
                                  ifo=workflow.ifo_string,
                                  out_dir=output_dir,
                                  tags=curr_tags)
    pycbc_genranking_job = pycbc_genranking_exe(workflow.cp, 
                                  pycbc_genranking_exe_name,
                                  ifo=workflow.ifo_string,
                                  out_dir=output_dir,
                                  tags=curr_tags)
    marg_likelihood_job_1 = marg_likelihood_exe(workflow.cp,
                                  marg_likelihood_exe_name,
                                  ifo=workflow.ifo_string,
                                  out_dir=output_dir,
                                  tags=['MARG1']+curr_tags)
    marg_likelihood_job_2 = marg_likelihood_exe(workflow.cp,
                                  marg_likelihood_exe_name,
                                  ifo=workflow.ifo_string,
                                  out_dir=output_dir,
                                  tags=['MARG2']+curr_tags)


    # Begin with finding the horizon distances
    picklehor_inputs = stage0_outputs[inj_less_tag]
    node = pycbc_picklehor_job.create_node(workflow.analysis_time,
                                                              picklehor_inputs)
    workflow.add_node(node)
    horizon_dist_file = node.output_files[0]
    # Then combine all likelihood files
    combllhood_inputs = likelihood_files.find_output_with_tag(\
                                                         pipedown_dq_veto_name) 
    combllhood_inputs = combllhood_inputs.find_output_with_tag(inj_less_tag)
    assert len(combllhood_inputs) > 0
    node = pycbc_combllhood_job.create_node(workflow.analysis_time,
                                          combllhood_inputs, horizon_dist_file)
    workflow.add_node(node)
    likelihood_file = node.output_files[0]
    # Also compute the ranking file
    node = pycbc_genranking_job.create_node(workflow.analysis_time,
                                            likelihood_file, horizon_dist_file)
    workflow.add_node(node)
    ranking_likelihood_file = node.output_files[0]
    # And marginalize (twice for some reason!)
    node = marg_likelihood_job_1.create_node(workflow.analysis_time,
                                                       ranking_likelihood_file)
    workflow.add_node(node)
    marg_likelihood_file_1 = node.output_files[0]
    node = marg_likelihood_job_2.create_node(workflow.analysis_time,
                                                        marg_likelihood_file_1)
    workflow.add_node(node)
    marg_likelihood_file_2 = node.output_files[0]

    # Now do the sqlite conditioning. This has a few stages.
                                                  
    # STAGE 1: Populate likelihood in all input files
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
    final_outputs = FileList([])
    # Do for all injection runs and zero lag
    for inj_tag in [inj_less_tag] + injection_tags:
        curr_tags = tags + [inj_tag, veto_tag]
        trig_inp_files = stage0_outputs[inj_tag]
        stage1_job = pycbc_compllhood_exe(workflow.cp,
                                      pycbc_compllhood_exe_name,
                                      ifo=workflow.ifo_string,
                                      out_dir=output_dir,
                                      tags=['STAGE1']+curr_tags)
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
            stage9_job = ligolw_sqlite_exe(workflow.cp, ligolw_sqlite_exe_name,
                                          ifo=workflow.ifo_string,
                                          out_dir=output_dir,
                                          tags=['FINAL'] + curr_tags)

        stage1_outputs[inj_tag] = FileList([])
        stage2_outputs[inj_tag] = FileList([])
        assert len(trig_inp_files) > 0
        for file in trig_inp_files:
            stage1_node = stage1_job.create_node(file.segment, file,
                                            likelihood_file, horizon_dist_file)
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
                                    stage2_outputs[inj_tag], workflow=workflow)
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
        stage5_inputs.append(dq_seg_file)
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

    # Next we run the compute FAR from snr_chisq histograms job
    far_gstlal_outputs = {}
    for inj_tag in [inj_less_tag] + injection_tags:
        curr_tags = tags + [inj_tag, veto_tag]
        far_gstlal_job = far_gstlal_exe(workflow.cp, far_gstlal_exe_name,
                                          ifo=workflow.ifo_string,
                                          out_dir=output_dir, tags=curr_tags)
        trig_veto_inp_files = \
                  final_outputs.find_output_with_tag(veto_tag)
        trig_inp_files = trig_veto_inp_files.find_output_with_tag(inj_tag)
        assert len(trig_inp_files) == 1
        input_database = trig_inp_files[0]
        if inj_tag != inj_less_tag:
            no_inj_db = trig_veto_inp_files.find_output_with_tag(inj_less_tag)
            assert len(no_inj_db) == 1
            no_inj_db = no_inj_db[0]
            write_background = False
        else:
            # Here I don't want to provide the same file as a dependancy
            # twice. Therefore I just give non-injection DB and the code
            # assumes this is also the input-database if it is not given.
            # Also, I only want the background file once
            no_inj_db =  input_database
            input_database = None
            write_background = True
        far_gstlal_node = far_gstlal_job.create_node(workflow.analysis_time,
                                        no_inj_db, marg_likelihood_file_2,
                                        inj_database=input_database,
                                        write_background_bins=write_background)
        workflow.add_node(far_gstlal_node)
        outputs = far_gstlal_node.output_files
        if inj_tag != inj_less_tag:
            assert len(outputs) == 1
            far_gstlal_outputs[inj_tag] = outputs[0]
        else:
            assert len(outputs) == 2
            sql_out = outputs.find_output_without_tag('POSTMARG')[0]
            xml_out = outputs.find_output_with_tag('POSTMARG')[0]
            far_gstlal_outputs[inj_tag] = sql_out
            post_marginalized_file = xml_out
            

    # Finally some plotting. 
    # FIXME: These are given explicit output directories and pegasus does not
    # know about output files. Would be nice if this was done "better"  
    curr_tags = tags + [veto_tag]
    plot_summary_job = plot_summary_exe(workflow.cp, plot_summary_exe_name,
                                          ifo=workflow.ifo_string,
                                          out_dir=output_dir, tags=curr_tags)
    plot_sensitivity_job = plot_sensitivity_exe(workflow.cp,
                                          plot_sensitivity_exe_name,
                                          ifo=workflow.ifo_string,
                                          out_dir=output_dir, tags=curr_tags)
    plot_background_job = plot_background_exe(workflow.cp,
                                          plot_background_exe_name,
                                          ifo=workflow.ifo_string,
                                          out_dir=output_dir, tags=curr_tags)
    inj_dbs = []
    for inj_tag in injection_tags:
        inj_dbs.append(far_gstlal_outputs[inj_tag])
    non_inj_db = far_gstlal_outputs[inj_less_tag]
    
    plot_summary_node = plot_summary_job.create_node(non_inj_db, inj_dbs)
    plot_background_node = plot_background_job.create_node(non_inj_db,
                                                        post_marginalized_file)
    plot_sensitivity_node = plot_sensitivity_job.create_node(non_inj_db,
                                                                       inj_dbs)

    workflow.add_node(plot_summary_node)
    workflow.add_node(plot_background_node)
    workflow.add_node(plot_sensitivity_node)

    # And make the html pages
    parents = [plot_summary_node, plot_background_node, plot_sensitivity_node]
    closed_summarypage_job = summary_page_exe(workflow.cp,
                                              summary_page_exe_name,
                                              ifo=workflow.ifo_string,
                                              out_dir=output_dir,
                                              tags=['CLOSEDBOX'] + curr_tags)
    open_summarypage_job = summary_page_exe(workflow.cp, 
                                              summary_page_exe_name,
                                              ifo=workflow.ifo_string,
                                              out_dir=output_dir,
                                              tags=['OPENBOX'] + curr_tags)

    closed_summarypage_node = closed_summarypage_job.create_and_add_node(\
                                              workflow, parents)
    open_summarypage_node = open_summarypage_job.create_and_add_node(workflow,
                                              parents)

    # FIXME: Maybe contatenate and return all other outputs if needed elsewhere
    # FIXME: Move to pp utils and return the FAR files.
    return final_outputs
