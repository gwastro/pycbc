from __future__ import division

import os
from glue import segments
from pycbc.ahope import select_job_instance,sngl_ifo_job_setup,AhopeOutFileList
from pycbc.ahope import AhopeOutFile

def setup_tmpltbank_workflow(cp,scienceSegs,ahopeDax):
    '''
    Setup template bank section of ahope workflow.
    FIXME: ADD MORE DOCUMENTATION
    '''

    # There should be a number of different options here, for e.g. to set
    # up fixed bank, or maybe something else
    
    # First thing is to check if a pre-generated bank is supplied:
    if (cp.has_option('ahope','pregenerated-template-bank')):
        tmpltBanks = setup_tmpltbank_pregenerated(cp,scienceSegs)
    # Else we assume template banks will be generated in the workflow
    else:
        tmpltBanks = setup_tmpltbank_dax_generated(cp,scienceSegs,ahopeDax)
    
    return tmpltBanks

def setup_tmpltbank_dax_generated(cp,scienceSegs,ahopeDax):
    '''
    Setup template bank jobs that are generated as part of the ahope workflow.
    FIXME: ADD MORE DOCUMENTATION
    '''

    # Need to get the exe to figure out what sections are analysed, what is
    # discarded etc. This should *not* be hardcoded, so using a new executable
    # will require a bit of effort here .... 
    # There is also a stub for a default class using values given in the .ini
    # file.

    ifos = scienceSegs.keys()
    tmpltBankExe = os.path.basename(cp.get('executables','tmpltbank'))
    # Select the appropriate class
    exeInstance = select_job_instance(tmpltBankExe,'tmpltbank')

    # Set up class for holding the banks
    tmpltBanks = AhopeOutFileList([])

    # Template banks are independent for different ifos, but might not be!
    # Begin with independent case and add after FIXME
    for ifo in ifos:
        sngl_ifo_job_setup(cp,ifo,tmpltBanks,exeInstance,scienceSegs[ifo],\
                           ahopeDax)

    return tmpltBanks

def setup_tmpltbank_pregenerated(cp,scienceSegs):
    '''
    Setup ahope workflow to use a pregenerated template bank.
    FIXME: ADD MORE DOCUMENTATION
    '''
    # Currently this uses the *same* fixed bank for all ifos.
    # Maybe we want to add capability to analyse separate banks in all ifos?
    ifos = scienceSegs.keys()
    
    # Set up class for holding the banks
    tmpltBanks = AhopeOutFileList([])

    preGenBank = cp.get('ahope','pregenerated-template-bank')
    globalSeg = segments.segment([0,9999999999])

    for ifo in ifos:
        # Add bank for that ifo
        currFile = AhopeOutFile()
        currFile.set_output(preGenBank)
        currFile.set_ifo(ifo)
        currFile.set_time(segments.segmentlist([globalSeg]))
        tmpltBanks.append(currFile)

    return tmpltBanks
