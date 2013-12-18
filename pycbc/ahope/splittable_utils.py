from __future__ import division

import os
import logging
from pycbc.ahope.ahope_utils import * 
from pycbc.ahope.jobsetup_utils import *

def setup_splittable_workflow(workflow, tmpltBanks, outDir):
    '''
    Setup matched filter section of ahope workflow.
    FIXME: ADD MORE DOCUMENTATION
    '''
    logging.info("Entering split output files module.")
    # Scope here for choosing different options
    logging.info("Adding split output file jobs to workflow.")
    splitTableOuts = setup_splittable_dax_generated(workflow, tmpltBanks,
                                                    outDir)
    logging.info("Leaving split output files module.")
    
    return splitTableOuts

def setup_splittable_dax_generated(workflow, tmpltBanks, outDir):
    '''
    Setup matched-filter jobs that are generated as part of the ahope workflow.
    FIXME: ADD MORE DOCUMENTATION
    '''
    cp = workflow.cp
    splittableExe = os.path.basename(cp.get('executables', 'splittable'))
    # Select the appropriate class
    exeInstance = select_splitfilejob_instance(splittableExe, 'splittable')

    # Template banks are independent for different ifos, but might not be!
    # Begin with independent case and add after FIXME
    # FIXME: Do not hardcode value of 2
    return split_outfiles(workflow, tmpltBanks, exeInstance, outDir)

def split_outfiles(workflow, inputFileList, exeInstance, outDir):
    """
    Add documentation
    """
    # Set up output structure
    outFileGroups = AhopeFileList([])

    # Set up the condorJob class for the current executable
    currExeJob = exeInstance.create_job(workflow.cp, None, outDir)

    for input in inputFileList:
        node = currExeJob.create_node(input)
        workflow.add_node(node)
        outFileGroups += node.output_files
    return outFileGroups

