from __future__ import division

import os
from pycbc.ahope import select_job_instance,sngl_ifo_job_setup,AhopeOutFileList

def setup_matchedfltr_workflow(cp,scienceSegs,ahopeDax,tmpltBanks):
    '''
    Setup matched filter section of ahope workflow.
    FIXME: ADD MORE DOCUMENTATION
    '''

    # There should be a number of different options here, for e.g. to set
    # up fixed bank, or maybe something else
    inspiralOuts = setup_matchedfltr_dax_generated(cp,scienceSegs,ahopeDax,\
                                                 tmpltBanks)
    
    return inspiralOuts

def setup_matchedfltr_dax_generated(cp,scienceSegs,ahopeDax,tmpltBanks):
    '''
    Setup matched-filter jobs that are generated as part of the ahope workflow.
    FIXME: ADD MORE DOCUMENTATION
    '''

    # Need to get the exe to figure out what sections are analysed, what is
    # discarded etc. This should *not* be hardcoded, so using a new executable
    # will require a bit of effort here .... 
    # There is also a stub for a default class using values given in the .ini
    # file.

    ifos = scienceSegs.keys()
    matchFltrExe = os.path.basename(cp.get('executables','inspiral'))
    # Select the appropriate class
    exeInstance = select_job_instance(matchFltrExe,'inspiral')

    # Set up class for holding the banks
    inspiralOuts = AhopeOutFileList([])

    # Template banks are independent for different ifos, but might not be!
    # Begin with independent case and add after FIXME
    for ifo in ifos:
        sngl_ifo_job_setup(cp,ifo,inspiralOuts,exeInstance,scienceSegs[ifo],\
                           ahopeDax,parents=tmpltBanks)

    return inspiralOuts
