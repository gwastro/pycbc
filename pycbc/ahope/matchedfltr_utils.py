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
This module is responsible for setting up the matched-filtering stage of ahope
workflows. For details about this module and its capabilities see here:
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/NOTYETCREATED.html
"""

from __future__ import division

import os
from pycbc.ahope.ahope_utils import * 
from pycbc.ahope.jobsetup_utils import *

def setup_matchedfltr_workflow(workflow, science_segs, datafind_outs,
                               tmplt_banks, output_dir=None, injection_file=None,
                               tags=[]):
    '''
    Setup matched filter section of ahope workflow.
    FIXME: ADD MORE DOCUMENTATION
    '''
    logging.info("Entering matched-filtering setup module.")
    make_analysis_dir(output_dir)
    # Scope here for choosing different options
    logging.info("Adding matched-filtering jobs to workflow.")

    # There should be a number of different options here, for e.g. to set
    # up fixed bank, or maybe something else
    inspiral_outs = setup_matchedfltr_dax_generated(workflow, science_segs, 
                                     datafind_outs, tmplt_banks, output_dir,
                                     injection_file, tags)
    logging.info("Leaving matched-filtering setup module.")    
    return inspiral_outs

def setup_matchedfltr_dax_generated(workflow, science_segs, datafind_outs,
                                    tmplt_banks, output_dir, injection_file=None,
                                    tags=[], link_to_tmpltbank=True):
    '''
    Setup matched-filter jobs that are generated as part of the ahope workflow.
    FIXME: ADD MORE DOCUMENTATION
    '''

    # Need to get the exe to figure out what sections are analysed, what is
    # discarded etc. This should *not* be hardcoded, so using a new executable
    # will require a bit of effort here .... 
    # There is also a stub for a default class using values given in the .ini
    # file.

    cp = workflow.cp
    ifos = science_segs.keys()
    match_fltr_exe = os.path.basename(cp.get('executables','inspiral'))
    # Select the appropriate class
    exe_instance = select_matchedfilterjob_instance(match_fltr_exe, 'inspiral')

    if link_to_tmpltbank:
        # Use this to ensure that inspiral and tmpltbank jobs overlap. This
        # means that there will be 1 inspiral job for every 1 tmpltbank and
        # the data read in by both will overlap as much as possible. (If you
        # ask the template bank jobs to use 2000s of data for PSD estimation
        # and the matched-filter jobs to use 4000s, you will end up with
        # twice as many matched-filter jobs that still use 4000s to estimate a
        # PSD but then only generate triggers in the 2000s of data that the
        # template bank jobs ran on.
        tmpltbank_exe = os.path.basename(cp.get('executables', 'tmpltbank'))
        link_exe_instance = select_tmpltbankjob_instance(tmpltbank_exe, 
                                                        'tmpltbank')
    else:
        link_exe_instance = None

    # Set up class for holding the banks
    inspiral_outs = AhopeFileList([])

    # Template banks are independent for different ifos, but might not be!
    # Begin with independent case and add after FIXME
    for ifo in ifos:
        job_instance = exe_instance.create_job(workflow.cp, ifo, output_dir, 
                                               injection_file=injection_file, 
                                               tags=tags)
        sngl_ifo_job_setup(workflow, ifo, inspiral_outs, job_instance, 
                           science_segs[ifo], datafind_outs, output_dir,
                           parents=tmplt_banks, 
                           link_exe_instance=link_exe_instance,
                           allow_overlap=False)
    return inspiral_outs
