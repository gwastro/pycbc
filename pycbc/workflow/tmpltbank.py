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
This module is responsible for setting up the template bank stage of CBC
workflows. For details about this module and its capabilities see here:
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/ahope/template_bank.html
"""

from __future__ import division

import os
import ConfigParser
import urlparse, urllib
import logging
from pycbc.workflow.core import File, FileList, make_analysis_dir, resolve_url
from pycbc.workflow.jobsetup import select_tmpltbank_class, select_matchedfilter_class, sngl_ifo_job_setup

def setup_tmpltbank_workflow(workflow, science_segs, datafind_outs,
                             output_dir=None, psd_files=None, tags=[]):
    '''
    Setup template bank section of CBC workflow. This function is responsible
    for deciding which of the various template bank workflow generation
    utilities should be used.

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
    psd_files : pycbc.workflow.core.FileList
        The file list containing predefined PSDs, if provided.
    tags : list of strings
        If given these tags are used to uniquely name and identify output files
        that would be produced in multiple calls to this function.

    Returns
    --------
    tmplt_banks : pycbc.workflow.core.FileList
        The FileList holding the details of all the template bank jobs.
    '''
    logging.info("Entering template bank generation module.")
    make_analysis_dir(output_dir)
    cp = workflow.cp
    
    # Parse for options in ini file
    tmpltbankMethod = cp.get_opt_tags("workflow-tmpltbank", "tmpltbank-method",
                                      tags)

    # There can be a large number of different options here, for e.g. to set
    # up fixed bank, or maybe something else
    if tmpltbankMethod == "PREGENERATED_BANK":
        logging.info("Setting template bank from pre-generated bank(s).")
        tmplt_banks = setup_tmpltbank_pregenerated(workflow, tags=tags)
    # Else we assume template banks will be generated in the workflow
    elif tmpltbankMethod == "WORKFLOW_INDEPENDENT_IFOS":
        logging.info("Adding template bank jobs to workflow.")
        if cp.has_option_tags("workflow-tmpltbank",
                              "tmpltbank-link-to-matchedfilter", tags):
            if not cp.has_option_tags("workflow-matchedfilter",
                              "matchedfilter-link-to-tmpltbank", tags):
                errMsg = "If using tmpltbank-link-to-matchedfilter, you should "
                errMsg = "also use matchedfilter-link-to-tmpltbank."
                logging.warn(errMsg)
            linkToMatchedfltr = True
        else:
            linkToMatchedfltr = False
        if cp.has_option_tags("workflow-tmpltbank",
                              "tmpltbank-compatibility-mode", tags):
            if not linkToMatchedfltr:
                errMsg = "Compatibility mode requires that the "
                errMsg += "tmpltbank-link-to-matchedfilter option is also set."
                raise ValueError(errMsg)
            if not cp.has_option_tags("workflow-matchedfilter",
                              "matchedfilter-compatibility-mode", tags):
                errMsg = "If using compatibility mode it must be set both in "
                errMsg += "the template bank and matched-filtering stages."
                raise ValueError(errMsg)
            compatibility_mode = True
        else:
            compatibility_mode = False
        tmplt_banks = setup_tmpltbank_dax_generated(workflow, science_segs,
                                         datafind_outs, output_dir, tags=tags,
                                         link_to_matchedfltr=linkToMatchedfltr,
                                         compatibility_mode=compatibility_mode,
                                         psd_files=psd_files)
    elif tmpltbankMethod == "WORKFLOW_INDEPENDENT_IFOS_NODATA":
        logging.info("Adding template bank jobs to workflow.")
        tmplt_banks = setup_tmpltbank_without_frames(workflow, output_dir,
                                         tags=tags, independent_ifos=True,
                                         psd_files=psd_files)
    elif tmpltbankMethod == "WORKFLOW_NO_IFO_VARIATION_NODATA":
        logging.info("Adding template bank jobs to workflow.")
        tmplt_banks = setup_tmpltbank_without_frames(workflow, output_dir,
                                         tags=tags, independent_ifos=False,
                                         psd_files=psd_files)
    else:
        errMsg = "Template bank method not recognized. Must be either "
        errMsg += "PREGENERATED_BANK, WORKFLOW_INDEPENDENT_IFOS "
        errMsg += "or WORKFLOW_INDEPENDENT_IFOS_NODATA."
        raise ValueError(errMsg)
    
    logging.info("Leaving template bank generation module.")
    return tmplt_banks

def setup_tmpltbank_dax_generated(workflow, science_segs, datafind_outs,
                                  output_dir, tags=[],
                                  link_to_matchedfltr=True,
                                  compatibility_mode=False,
                                  psd_files=None):
    '''
    Setup template bank jobs that are generated as part of the CBC workflow.
    This function will add numerous jobs to the CBC workflow using
    configuration options from the .ini file. The following executables are
    currently supported:

    * lalapps_tmpltbank
    * pycbc_geom_nonspin_bank

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
    link_to_matchedfltr : boolean, optional (default=True)
        If this option is given, the job valid_times will be altered so that
        there will be one inspiral file for every template bank and they will
        cover the same time span. Note that this option must also be given
        during matched-filter generation to be meaningful.
    psd_file : pycbc.workflow.core.FileList
        The file list containing predefined PSDs, if provided.

    Returns
    --------
    tmplt_banks : pycbc.workflow.core.FileList
        The FileList holding the details of all the template bank jobs.
    '''
    cp = workflow.cp
    # Need to get the exe to figure out what sections are analysed, what is
    # discarded etc. This should *not* be hardcoded, so using a new executable
    # will require a bit of effort here .... 

    ifos = science_segs.keys()
    tmplt_bank_exe = os.path.basename(cp.get('executables', 'tmpltbank'))
    # Select the appropriate class
    exe_class = select_tmpltbank_class(tmplt_bank_exe)

    # The exe instance needs to know what data segments are analysed, what is
    # discarded etc. This should *not* be hardcoded, so using a new executable
    # will require a bit of effort here .... 

    if link_to_matchedfltr:
        # Use this to ensure that inspiral and tmpltbank jobs overlap. This
        # means that there will be 1 inspiral job for every 1 tmpltbank and
        # the data read in by both will overlap as much as possible. (If you
        # ask the template bank jobs to use 2000s of data for PSD estimation
        # and the matched-filter jobs to use 4000s, you will end up with
        # twice as many matched-filter jobs that still use 4000s to estimate a
        # PSD but then only generate triggers in the 2000s of data that the
        # template bank jobs ran on.
        tmpltbank_exe = os.path.basename(cp.get('executables', 'inspiral'))
        link_exe_instance = select_matchedfilter_class(tmpltbank_exe)
    else:
        link_exe_instance = None

    # Set up class for holding the banks
    tmplt_banks = FileList([])


    # Template banks are independent for different ifos, but might not be!
    # Begin with independent case and add after FIXME
    for ifo in ifos:
        job_instance = exe_class(workflow.cp, 'tmpltbank', ifo=ifo, 
                                               out_dir=output_dir,
                                               tags=tags)
        # Check for the write_psd flag
        if cp.has_option_tags("workflow-tmpltbank", "tmpltbank-write-psd-file", tags):
            job_instance.write_psd = True
        else:
            job_instance.write_psd = False

        if link_exe_instance:
            link_job_instance = link_exe_instance(cp, 'inspiral', ifo=ifo,
                        out_dir=output_dir, tags=tags)
        else:
            link_job_instance = None
        sngl_ifo_job_setup(workflow, ifo, tmplt_banks, job_instance, 
                           science_segs[ifo], datafind_outs,
                           link_job_instance=link_job_instance, 
                           allow_overlap=True,
                           compatibility_mode=compatibility_mode)
    return tmplt_banks

def setup_tmpltbank_without_frames(workflow, output_dir,
                                   tags=[], independent_ifos=False,
                                   psd_files=None):
    '''
    Setup CBC workflow to use a template bank (or banks) that are generated in
    the workflow, but do not use the data to estimate a PSD, and therefore do
    not vary over the duration of the workflow. This can either generate one
    bank that is valid for all ifos at all times, or multiple banks that are
    valid only for a single ifo at all times (one bank per ifo).

    Parameters
    ----------
    workflow: pycbc.workflow.core.Workflow
        An instanced class that manages the constructed workflow.
    output_dir : path string
        The directory where the template bank outputs will be placed. 
    tags : list of strings
        If given these tags are used to uniquely name and identify output files
        that would be produced in multiple calls to this function.
    independent_ifos : Boolean, optional (default=False)
        If given this will produce one template bank per ifo. If not given
        there will be on template bank to cover all ifos.
    psd_file : pycbc.workflow.core.FileList
        The file list containing predefined PSDs, if provided.

    Returns
    --------
    tmplt_banks : pycbc.workflow.core.FileList
        The FileList holding the details of the template bank(s).
    '''
    cp = workflow.cp
    # Need to get the exe to figure out what sections are analysed, what is
    # discarded etc. This should *not* be hardcoded, so using a new executable
    # will require a bit of effort here .... 

    ifos = workflow.ifos
    fullSegment = workflow.analysis_time

    tmplt_bank_exe = os.path.basename(cp.get('executables','tmpltbank'))
    # Can not use lalapps_template bank with this
    if tmplt_bank_exe == 'lalapps_tmpltbank':
        errMsg = "Lalapps_tmpltbank cannot be used to generate template banks "
        errMsg += "without using frames. Try another code."
        raise ValueError(errMsg)

    # Select the appropriate class
    exe_instance = select_tmpltbank_class(tmplt_bank_exe)

    tmplt_banks = FileList([])

    # Make the distinction between one bank for all ifos and one bank per ifo
    if independent_ifos:
        ifoList = [ifo for ifo in ifos]
    else:
        ifoList = [[ifo for ifo in ifos]]

    # Check for the write_psd flag
    if cp.has_option_tags("workflow-tmpltbank", "tmpltbank-write-psd-file", tags):
        exe_instance.write_psd = True
    else:
        exe_instance.write_psd = False

    for ifo in ifoList:
        job_instance = exe_instance(workflow.cp, 'tmpltbank', ifo=ifo, 
                                               out_dir=output_dir,
                                               tags=tags,
                                               psd_files=psd_files)
        node = job_instance.create_nodata_node(fullSegment)
        workflow.add_node(node)
        tmplt_banks += node.output_files
        
    return tmplt_banks

def setup_tmpltbank_pregenerated(workflow, tags=[]):
    '''
    Setup CBC workflow to use a pregenerated template bank.
    The bank given in cp.get('workflow','pregenerated-template-bank') will be used
    as the input file for all matched-filtering jobs. If this option is
    present, workflow will assume that it should be used and not generate
    template banks within the workflow.

    Parameters
    ----------
    workflow: pycbc.workflow.core.Workflow
        An instanced class that manages the constructed workflow.
    tags : list of strings
        If given these tags are used to uniquely name and identify output files
        that would be produced in multiple calls to this function.

    Returns
    --------
    tmplt_banks : pycbc.workflow.core.FileList
        The FileList holding the details of the template bank.
    '''
    # Currently this uses the *same* fixed bank for all ifos.
    # Maybe we want to add capability to analyse separate banks in all ifos?
    
    # Set up class for holding the banks
    tmplt_banks = FileList([])

    cp = workflow.cp
    global_seg = workflow.analysis_time
    user_tag = "PREGEN_TMPLTBANK"
    try:
        # First check if we have a bank for all ifos
        pre_gen_bank = cp.get_opt_tags('workflow-tmpltbank',
                                           'tmpltbank-pregenerated-bank', tags)
        pre_gen_bank = resolve_url(pre_gen_bank)
        file_url = urlparse.urljoin('file:', urllib.pathname2url(pre_gen_bank))
        curr_file = File(workflow.ifos, user_tag, global_seg, file_url,
                                                                     tags=tags)
        curr_file.PFN(file_url, site='local')
        tmplt_banks.append(curr_file)
    except ConfigParser.Error:
        # Okay then I must have banks for each ifo
        for ifo in workflow.ifos:
            try:
                pre_gen_bank = cp.get_opt_tags('workflow-tmpltbank',
                                'tmpltbank-pregenerated-bank-%s' % ifo.lower(),
                                tags)
                pre_gen_bank = resolve_url(pre_gen_bank)
                file_url = urlparse.urljoin('file:',
                                             urllib.pathname2url(pre_gen_bank))
                curr_file = File(ifo, user_tag, global_seg, file_url,
                                                                     tags=tags)
                curr_file.PFN(file_url, site='local')
                tmplt_banks.append(curr_file)

            except ConfigParser.Error:
                err_msg = "Cannot find pregerated template bank in section "
                err_msg += "[workflow-tmpltbank] or any tagged sections. "
                if tags:
                    tagged_secs = " ".join("[workflow-tmpltbank-%s]" \
                                           %(ifo,) for ifo in workflow.ifos)
                    err_msg += "Tagged sections are %s. " %(tagged_secs,)
                err_msg += "I looked for 'tmpltbank-pregenerated-bank' option "
                err_msg += "and 'tmpltbank-pregenerated-bank-%s'." %(ifo,)
                raise ConfigParser.Error(err_msg)
            
    return tmplt_banks

