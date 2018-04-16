# Copyright (C) 2015 Larne Pekowsky
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
This module is responsible for setting up the gating files used by CBC
workflows.
"""

from __future__ import division

import os
import ConfigParser
import urlparse, urllib
import logging
from pycbc.workflow.core import File, FileList, make_analysis_dir, resolve_url

def setup_gating_workflow(workflow, output_dir=None, tags=None):
    '''
    Setup gating section of CBC workflow. At present this only supports pregenerated
    gating files, in the future these could be created within the workflow.

    Parameters
    ----------
    workflow: pycbc.workflow.core.Workflow
        An instanced class that manages the constructed workflow.
    output_dir : path string
        The directory where data products will be placed.
    tags : list of strings
        If given these tags are used to uniquely name and identify output files
        that would be produced in multiple calls to this function.

    Returns
    --------
    gate_files : pycbc.workflow.core.FileList
        The FileList holding the gate files, 0 or 1 per ifo
    '''
    if tags is None:
        tags = []
    logging.info("Entering gating module.")
    make_analysis_dir(output_dir)
    cp = workflow.cp

    # Parse for options in ini file.
    try:
        gateMethod = cp.get_opt_tags("workflow-gating", "gating-method",
                                     tags)
    except ConfigParser.Error:
        # Gating is optional, just return an empty list if not
        # provided.
        return FileList([])

    if gateMethod == "PREGENERATED_FILE":
        logging.info("Setting gating from pre-generated file(s).")
        gate_files = setup_gate_pregenerated(workflow,
                                             output_dir=output_dir, tags=tags)
    else:
        errMsg = "Gating method not recognized. Only "
        errMsg += "PREGENERATED_FILE is currently supported."
        raise ValueError(errMsg)

    logging.info("Leaving gating module.")
    return gate_files


def setup_gate_pregenerated(workflow, output_dir=None, tags=None):
    '''
    Setup CBC workflow to use pregenerated gating files.
    The file given in cp.get('workflow','gating-file-(ifo)') will
    be used as the --gating-file for all jobs for that ifo.

    Parameters
    ----------
    workflow: pycbc.workflow.core.Workflow
        An instanced class that manages the constructed workflow.
    output_dir : path string
       The directory where data products will be placed.
    tags : list of strings
        If given these tags are used to uniquely name and identify output files
        that would be produced in multiple calls to this function.

    Returns
    --------
    gate_files : pycbc.workflow.core.FileList
        The FileList holding the gating files
    '''
    if tags is None:
        tags = []
    gate_files = FileList([])

    cp = workflow.cp
    global_seg = workflow.analysis_time
    user_tag = "PREGEN_GATE"

    for ifo in workflow.ifos:
        try:
            pre_gen_file = cp.get_opt_tags('workflow-gating',
                            'gating-file-%s' % ifo.lower(),
                            tags)
            pre_gen_file = resolve_url(pre_gen_file,
                                       os.path.join(os.getcwd(),output_dir))
            file_url = urlparse.urljoin('file:',
                                         urllib.pathname2url(pre_gen_file))
            curr_file = File(ifo, user_tag, global_seg, file_url,
                                                                 tags=tags)
            curr_file.PFN(file_url, site='local')
            gate_files.append(curr_file)

            logging.info("Using gating file %s for %s", file_url, ifo)

        except ConfigParser.Error:
            logging.info("No gating file specified for %s", ifo)

    return gate_files

