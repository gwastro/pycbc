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
from glue import segments
from pycbc.ahope.ahope_utils import *
from pycbc.ahope.jobsetup_utils import *

def setup_postprocessing_preperation(workflow, triggerFiles, output_dir,
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
    Workflow : ahope.Workflow
        The ahope workflow instance that the coincidence jobs will be added to.
    triggerFiles : ahope.AhopeFileList
        An AhopeFileList of the trigger files that are used as
        input at this stage.
    output_dir : path
        The directory in which output files will be stored.
    tags : list of strings (optional, default = [])
        A list of the tagging strings that will be used for all jobs created
        by this call to the workflow. An example might be ['POSTPROC1'] or
        ['DENTYSNEWPOSTPROC']. This will be used in output names.

    Returns
    --------
    postProcPreppedFiles : ahope.AhopeFileList
        A list of files that can be used as input for the post-processing stage.
    """
    logging.info("Entering post-processing preperation stage.")
    make_analysis_dir(output_dir)

    # Parse for options in .ini file
    postProcPrepMethod = workflow.cp.get_opt_tags("ahope-postprocprep",
                                        "postprocprep-method", tags)

    # Scope here for adding different options/methods here. For now we only
    # have the single_stage ihope method which consists of converting the
    # ligolw_thinca output xml into one file, clustering, performing injection
    # finding and putting everything into one SQL database.
    if postProcPrepMethod == "PIPEDOWN_AHOPE":
        # If you want the intermediate output files, call this directly
        postPostPreppedFiles,_,_,_ = setup_postprocprep_pipedown_ahope(workflow,
                           triggerFiles, output_dir,
                           tags=tags, **kwargs) 
    else:
        errMsg = "Post-processing preperation method not recognized. Must be "
        errMsg += "one of PIPEDOWN_AHOPE (currently only one option)."
        raise ValueError(errMsg)

    logging.info("Leaving post-processing separation module.")

    return postPostPreppedFiles

def setup_postprocprep_pipedown_ahope(workflow, coincFiles, output_dir,
                                      tags=[], injectionFiles=None,
                                      vetoFiles=None, injLessTag=None,
                                      injectionTags=[], vetoCats=[]):
    """
    Properties
    -----------
    Workflow : ahope.Workflow
        The ahope workflow instance that the coincidence jobs will be added to.
    coincFiles : ahope.AhopeFileList
        An AhopeFileList of the coincident trigger files that are used as
        input at this stage.
    output_dir : path
        The directory in which output files will be stored.
    tags : list of strings (optional, default = [])
        A list of the tagging strings that will be used for all jobs created
        by this call to the workflow. An example might be ['POSTPROC1'] or
        ['DENTYSNEWPOSTPROC']. This will be used in output names.
    injectionFiles : ahope.AhopeFileList (optional, default=None)
        The injection files to be used in this stage. An empty list (or any
        other input that evaluates as false) is valid and will imply that no
        injections are being done.
    vetoFiles : ahope.AhopeFileList (required)
        The data quality files to be used in this stage. This is required and
        will be used to determine the analysed times when doing post-processing.
    injLessTag : string (required)
        The tag that identifies files that do not have simulations in them.
        Ie. the primary search results.
    injectionTags : list of strings (optional, default = [])
        Each injection file has a unique tag. If used in the method, this
        tells the post-processing preperation code which injection tags it
        should include when creating the combined output.
    vetoCats : list of integers (optional, default = [])
        Decide which set of veto files should be used in the post-processing
        preparation. This is used, for example, to tell the workflow that you
        are only interested in quoting results at categories 2, 3 and 4. In
        which case just supply [2,3,4] for those veto files here.

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
    if not vetoCats:
        raise ValueError("Veto cats is required.")

    # Setup needed exe classes
    sqliteCombine1ExeTag = workflow.cp.get_opt_tags("ahope-postprocprep",
                                   "postprocprep-combiner1-exe", tags)
    sqliteCombine1Exe = select_genericjob_instance(workflow,\
                                                   sqliteCombine1ExeTag)
    sqliteCombine2ExeTag = workflow.cp.get_opt_tags("ahope-postprocprep",
                                   "postprocprep-combiner2-exe", tags)
    sqliteCombine2Exe = select_genericjob_instance(workflow,\
                                                   sqliteCombine2ExeTag)
    clusterCoincsExeTag = workflow.cp.get_opt_tags("ahope-postprocprep",
                                   "postprocprep-cluster-exe", tags)
    clusterCoincsExe = select_genericjob_instance(workflow, clusterCoincsExeTag)
    injFindExeTag = workflow.cp.get_opt_tags("ahope-postprocprep",
                                   "postprocprep-injfind-exe", tags)
    injFindExe = select_genericjob_instance(workflow, injFindExeTag)

    sqliteCombine1Outs = AhopeFileList([])
    clusterCoincsOuts = AhopeFileList([])
    injFindOuts = AhopeFileList([])
    sqliteCombine2Outs = AhopeFileList([])
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

        sqliteCombine2Inputs = AhopeFileList([])
        # Do injection-less jobs first.

        # Combine trig files first
        currTags = tags + [injLessTag, vetoTag]
        trigVetoInpFiles = coincFiles.find_output_with_tag(pipedownDQVetoName)
        trigInpFiles = trigVetoInpFiles.find_output_with_tag(injLessTag)
        trigInpFiles.append(dqSegFile[0])
        sqliteCombine1Job = sqliteCombine1Exe.create_job(workflow.cp, 
                          workflow.ifoString, out_dir=output_dir,tags=currTags)
        sqliteCombine1Node = sqliteCombine1Job.create_node(\
                                          workflow.analysis_time, trigInpFiles)
        workflow.add_node(sqliteCombine1Node)
        # Node has only one output file
        sqliteCombine1Out = sqliteCombine1Node.output_files[0]
        sqliteCombine1Outs.append(sqliteCombine1Out)

        # Cluster coincidences
        clusterCoincsJob = clusterCoincsExe.create_job(workflow.cp, 
                         workflow.ifoString, out_dir=output_dir, tags=currTags)
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
            sqliteCombine1Job = sqliteCombine1Exe.create_job(workflow.cp,
                          workflow.ifoString, out_dir=output_dir,tags=currTags)
            sqliteCombine1Node = sqliteCombine1Job.create_node(\
                                          workflow.analysis_time, trigInpFiles,
                                          injFile=injFile[0], injString=injTag)
            workflow.add_node(sqliteCombine1Node)
            # Node has only one output file
            sqliteCombine1Out = sqliteCombine1Node.output_files[0]
            sqliteCombine1Outs.append(sqliteCombine1Out)

            # Cluster coincidences
            clusterCoincsJob = clusterCoincsExe.create_job(workflow.cp,
                         workflow.ifoString, out_dir=output_dir, tags=currTags)
            clusterCoincsNode = clusterCoincsJob.create_node(\
                                     workflow.analysis_time, sqliteCombine1Out)
            workflow.add_node(clusterCoincsNode)
            # Node has only one output file
            clusterCoincsOut = clusterCoincsNode.output_files[0]
            clusterCoincsOuts.append(clusterCoincsOut)
            sqliteCombine2Inputs.append(clusterCoincsOut)

        # Combine everything together and add veto file
        currTags = tags + [vetoTag]
        sqliteCombine2Job = sqliteCombine2Exe.create_job(workflow.cp,
                          workflow.ifoString, out_dir=output_dir,tags=currTags)
        sqliteCombine2Node = sqliteCombine2Job.create_node(\
                                  workflow.analysis_time, sqliteCombine2Inputs)
        workflow.add_node(sqliteCombine2Node)
        sqliteCombine2Out = sqliteCombine2Node.output_files[0]
        sqliteCombine2Outs.append(sqliteCombine2Out)

        # Inj finding
        injFindJob = injFindExe.create_job(workflow.cp, workflow.ifoString,
                                          out_dir=output_dir,tags=currTags)
        injFindNode = injFindJob.create_node(workflow.analysis_time,
                                                         sqliteCombine2Out)
        workflow.add_node(injFindNode)
        injFindOut = injFindNode.output_files[0]
        injFindOuts.append(injFindOut)


    return injFindOuts, sqliteCombine1Outs, clusterCoincsOuts,\
           sqliteCombine2Outs


