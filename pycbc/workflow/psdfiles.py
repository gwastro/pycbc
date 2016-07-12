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
This module is responsible for setting up the psd files used by CBC
workflows. 
"""

from __future__ import division

import os
import ConfigParser
import urlparse, urllib
import logging
from pycbc.workflow.core import File, FileList, make_analysis_dir, resolve_url

def setup_psd_workflow(workflow, science_segs, datafind_outs,
                             output_dir=None, tags=None):
    '''
    Setup static psd section of CBC workflow. At present this only supports pregenerated
    psd files, in the future these could be created within the workflow.

    Parameters
    ----------
    workflow: pycbc.workflow.core.Workflow
        An instanced class that manages the constructed workflow.
    science_segs : Keyed dictionary of glue.segmentlist objects
        scienceSegs[ifo] holds the science segments to be analysed for each
        ifo. 
    datafind_outs : pycbc.workflow.core.FileList
        The file list containing the datafind files.
    output_dir : path string
        The directory where data products will be placed. 
    tags : list of strings
        If given these tags are used to uniquely name and identify output files
        that would be produced in multiple calls to this function.

    Returns
    --------
    psd_files : pycbc.workflow.core.FileList
        The FileList holding the psd files, 0 or 1 per ifo
    '''
    if tags is None:
        tags = []
    logging.info("Entering static psd module.")
    make_analysis_dir(output_dir)
    cp = workflow.cp
    
    # Parse for options in ini file.  
    try:
        psdMethod = cp.get_opt_tags("workflow-psd", "psd-method",
                                     tags)
    except:
        # Predefined PSD sare optional, just return an empty list if not
        # provided.
        return FileList([])

    if psdMethod == "PREGENERATED_FILE":
        logging.info("Setting psd from pre-generated file(s).")
        psd_files = setup_psd_pregenerated(workflow, tags=tags)
    else:
        errMsg = "PSD method not recognized. Only "
        errMsg += "PREGENERATED_FILE is currently supported."
        raise ValueError(errMsg)
    
    logging.info("Leaving psd module.")
    return psd_files


def setup_psd_pregenerated(workflow, tags=None):
    '''
    Setup CBC workflow to use pregenerated psd files.
    The file given in cp.get('workflow','pregenerated-psd-file-(ifo)') will 
    be used as the --psd-file argument to geom_nonspinbank, geom_aligned_bank
    and pycbc_plot_psd_file.

    Parameters
    ----------
    workflow: pycbc.workflow.core.Workflow
        An instanced class that manages the constructed workflow.
    tags : list of strings
        If given these tags are used to uniquely name and identify output files
        that would be produced in multiple calls to this function.

    Returns
    --------
    psd_files : pycbc.workflow.core.FileList
        The FileList holding the gating files
    '''
    if tags is None:
        tags = []
    psd_files = FileList([])

    cp = workflow.cp
    global_seg = workflow.analysis_time
    user_tag = "PREGEN_PSD"

    # Check for one psd for all ifos
    try:
        pre_gen_file = cp.get_opt_tags('workflow-psd',
                        'psd-pregenerated-file', tags)
        pre_gen_file = resolve_url(pre_gen_file)
        file_url = urlparse.urljoin('file:',
                                     urllib.pathname2url(pre_gen_file))
        curr_file = File(workflow.ifos, user_tag, global_seg, file_url,
                                                    tags=tags)
        curr_file.PFN(file_url, site='local')
        psd_files.append(curr_file)
    except ConfigParser.Error:
        # Check for one psd per ifo
        for ifo in workflow.ifos:
            try:
                pre_gen_file = cp.get_opt_tags('workflow-psd',
                                'psd-pregenerated-file-%s' % ifo.lower(),
                                tags)
                pre_gen_file = resolve_url(pre_gen_file)
                file_url = urlparse.urljoin('file:',
                                             urllib.pathname2url(pre_gen_file))
                curr_file = File(ifo, user_tag, global_seg, file_url,
                                                            tags=tags)
                curr_file.PFN(file_url, site='local')
                psd_files.append(curr_file)

            except ConfigParser.Error:
                # It's unlikely, but not impossible, that only some ifos
                # will have pregenerated PSDs
                logging.warn("No psd file specified for IFO %s." % (ifo,))
                pass
            
    return psd_files

