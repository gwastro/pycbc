from __future__ import division

import os
import logging
from pycbc.ahope.ahope_utils import * 
from pycbc.ahope.jobsetup_utils import *

def setup_splittable_workflow(workflow, tmplt_banks, out_dir=None):
    '''
    Setup matched filter section of ahope workflow.
    FIXME: ADD MORE DOCUMENTATION
    '''
    logging.info("Entering split output files module.")
    make_analysis_dir(out_dir)
    # Scope here for choosing different options
    logging.info("Adding split output file jobs to workflow.")
    split_table_outs = setup_splittable_dax_generated(workflow, tmplt_banks,
                                                    out_dir)
    logging.info("Leaving split output files module.")  
    return split_table_outs

def setup_splittable_dax_generated(workflow, tmplt_banks, out_dir):
    '''
    Setup matched-filter jobs that are generated as part of the ahope workflow.
    FIXME: ADD MORE DOCUMENTATION
    '''
    cp = workflow.cp
    splittable_exe = os.path.basename(cp.get('executables', 'splittable'))
    # Select the appropriate class
    exe_instance = select_splitfilejob_instance(splittable_exe, 'splittable')

    # Template banks are independent for different ifos, but might not be!
    # Begin with independent case and add after FIXME
    return split_outfiles(workflow, tmplt_banks, exe_instance, out_dir)

def split_outfiles(workflow, input_file_list, exe_instance, out_dir):
    """
    Add documentation
    """
    # Set up output structure
    out_file_groups = AhopeFileList([])

    # Set up the condorJob class for the current executable
    curr_exe_job = exe_instance.create_job(workflow.cp, None, out_dir)

    for input in input_file_list:
        node = curr_exe_job.create_node(input)
        workflow.add_node(node)
        out_file_groups += node.output_files
    return out_file_groups

