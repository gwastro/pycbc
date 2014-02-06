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
This module is responsible for setting up the splitting output files stage of ahope
workflows. For details about this module and its capabilities see here:
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/NOTYETCREATED.html
"""


from __future__ import division

import os
import logging
from pycbc.ahope.ahope_utils import * 
from pycbc.ahope.jobsetup_utils import *

def setup_splittable_workflow(workflow, tmplt_banks, out_dir=None):
    '''
    This function aims to be the gateway for code that is responsible for taking
    some input file containing some table, and splitting into multiple files
    containing different parts of that table. For now the only supported operation
    is using lalapps_splitbank to split a template bank xml file into multiple
    template bank xml files.

    Parameters
    -----------
    Workflow : ahope.Workflow
        The ahope workflow instance that the jobs will be added to.
    tmplt_banks : ahope.AhopeFileList
        The input files to be split up.
    out_dir : path
        The directory in which output will be written.

    Returns
    --------
    split_table_outs : ahope.AhopeFileList
        The list of split up files as output from this job.
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
    Function for setting up the splitting jobs as part of the workflow.

    Parameters
    -----------
    Workflow : ahope.Workflow
        The ahope workflow instance that the jobs will be added to.
    tmplt_banks : ahope.AhopeFileList
        The input files to be split up.
    out_dir : path
        The directory in which output will be written.

    Returns
    --------
    split_table_outs : ahope.AhopeFileList
        The list of split up files as output from this job.
    '''
    cp = workflow.cp
    splittable_exe = os.path.basename(cp.get('executables', 'splittable'))
    # Select the appropriate class
    exe_instance = select_splitfilejob_instance(splittable_exe, 'splittable')

    # Set up output structure
    out_file_groups = AhopeFileList([])

    # Set up the condorJob class for the current executable
    curr_exe_job = exe_instance.create_job(workflow.cp, None, out_dir)

    for input in tmplt_banks:
        node = curr_exe_job.create_node(input)
        workflow.add_node(node)
        out_file_groups += node.output_files
    return out_file_groups

