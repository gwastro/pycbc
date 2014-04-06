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

def setup_postprocessing(workflow, triggerFiles, output_dir, tags=[], **kwargs):
    """
    This function aims to be the gateway for running postprocessing in CBC
    offline workflows. Post-processing generally consists of calculating the
    significance of triggers and making any statements about trigger rates.
    Dedicated plotting jobs do not belong here.

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
    postProcFiles : ahope.AhopeFileList
        A list of the output from this stage.

    """
    logging.info("Entering post-processing preparation stage.")
    make_analysis_dir(output_dir)

    # Parse for options in .ini file
    postProcMethod = workflow.cp.get_opt_tags("ahope-postproc",
                                        "postproc-method", tags)

    # Scope here for adding different options/methods here. For now we only
    # have the single_stage ihope method which consists of converting the
    # ligolw_thinca output xml into one file, clustering, performing injection
    # finding and putting everything into one SQL database.
    if postProcMethod == "PIPEDOWN_AHOPE":
        # If you want the intermediate output files, call this directly
        postProcFiles = setup_postproc_pipedown_ahope(workflow,
                           triggerFiles, output_dir,
                           tags=tags, **kwargs)
    else:
        errMsg = "Post-processing method not recognized. Must be "
        errMsg += "one of PIPEDOWN_AHOPE (currently only one option)."
        raise ValueError(errMsg)

    logging.info("Leaving post-processing module.")

    return postProcFiles

def setup_postproc_pipedown_ahope(workflow, triggerFiles, output_dir, tags=[],
                                  vetoCats=[]):
    """
    This module sets up the post-processing stage in ahope, using a pipedown
    style set up. This consists of running compute_durations to determine and
    store the analaysis time (foreground and background). It then runs cfar
    jobs to determine the false alarm rate for all triggers (simulations or
    otherwise) in the input database.
    Pipedown expects to take as input (at this stage) a single database
    containing all triggers. This sub-module follows that same idea, so
    len(triggerFiles) must equal 1 (for every DQ category that we will run).

    Workflow : ahope.Workflow
        The ahope workflow instance that the coincidence jobs will be added to.
    triggerFiles : ahope.AhopeFileList
        An AhopeFileList containing the combined databases at CAT_1,2,3... that
        will be used to calculate FARs
    output_dir : path
        The directory in which output files will be stored.
    tags : list of strings (optional, default = [])
        A list of the tagging strings that will be used for all jobs created
        by this call to the workflow. An example might be ['POSTPROC1'] or
        ['DENTYSNEWPOSTPROC']. This will be used in output names.
    vetoCats : list of integers (optional, default = [])
        Decide which set of veto files should be used in the post-processing
        preparation. This is used, for example, to tell the workflow that you
        are only interested in quoting results at categories 2, 3 and 4. In
        which case just supply [2,3,4] for those veto files here.
   
    Returns
    --------
    finalFiles : ahope.AhopeFileList
        A list of the final SQL databases containing computed FARs.
    """
    if not vetoCats:
        raise ValueError("Veto cats is required.")

    # Setup needed exe classes
    computeDurationsExeTag = workflow.cp.get_opt_tags("ahope-postproc",
                                   "postproc-computedurations-exe", tags)
    computeDurationsExe = select_genericjob_instance(workflow,\
                                                     computeDurationsExeTag)
    cfarExeTag = workflow.cp.get_opt_tags("ahope-postproc", "postproc-cfar-exe",
                                       tags)
    cfarExe = select_genericjob_instance(workflow, cfarExeTag) 

    compDurationsOuts = AhopeFileList([])
    cfarOuts = AhopeFileList([])

    for vetoCat in vetoCats:
        vetoTag = 'CUMULATIVE_CAT_%d' %(vetoCat,)

        trigInpFiles = triggerFiles.find_output_with_tag(vetoTag)
        if not len(trigInpFiles) == 1:
            errMsg = "Did not find exactly 1 database input file."
            raise ValueError(errMsg)

        currTags = tags + [vetoTag]

        # Start with compute durations
        computeDurationsJob = computeDurationsExe.create_job(workflow.cp,
                         workflow.ifoString, out_dir=output_dir, tags=currTags)
        computeDurationsNode = computeDurationsJob.create_node(\
                                     workflow.analysis_time, trigInpFiles[0])
        workflow.add_node(computeDurationsNode)

        # Node has only one output file
        computeDurationsOut = computeDurationsNode.output_files[0]
        compDurationsOuts.append(computeDurationsOut)

        # Add the calculate FAR (cfar) job
        cfarJob = cfarExe.create_job(workflow.cp,
                         workflow.ifoString, out_dir=output_dir, tags=currTags)
        cfarNode = cfarJob.create_node(workflow.analysis_time,
                                       computeDurationsOut)
        workflow.add_node(cfarNode)

        # Node has only one output file
        cfarOut = cfarNode.output_files[0]
        cfarOuts.append(cfarOut)

    return cfarOuts
