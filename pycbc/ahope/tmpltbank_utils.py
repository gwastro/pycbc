from __future__ import division

import os
from glue import segments
from pycbc.ahope import select_job_instance,sngl_ifo_job_setup,AhopeOutFileList
from pycbc.ahope import AhopeOutFile

def setup_tmpltbank_workflow(cp,scienceSegs,ahopeDax):
    '''
    Setup template bank section of ahope workflow. This function is responsible
    for deciding which of the various template bank workflow generation
    utilities should be used. Currently the logic is as follows:
    
    * If pregenerated-template-bank=BANK option in [ahope] section, use
      setup_tmpltbank_pregenerated. This will assume a pregenerated bank called
      BANK is to be used for **all** ahope matched-filtering.

    * Otherwise use setup_tmpltbank_dax_generated. This will assume that the
      template banks are to be generated within the workflow. This will
      generate numerous banks for every ifo, dependent on configuration
      options.

    Parameters
    ----------
    cp : ConfigParser
        The ConfigParser holding all the options used by the ahope workflow.
    scienceSegs : Keyed dictionary of glue.segmentlist objects
        scienceSegs[ifo] holds the science segments to be analysed for each
        ifo. 
    ahopeDax : Instanced CondorDag class
        The CondorDag class that will hold all the jobs that the ahope workflow
        needs to run.

    Returns
    --------
    AhopeOutFileList
        The AhopeOutFileList holding the details of all the template bank jobs.
    '''

    # There should be a number of different options here, for e.g. to set
    # up fixed bank, or maybe something else
    
    # First thing is to check if a pre-generated bank is supplied:
    if (cp.has_option('ahope','pregenerated-template-bank')):
        tmpltBanks = setup_tmpltbank_pregenerated(cp,scienceSegs.keys())
    # Else we assume template banks will be generated in the workflow
    else:
        tmpltBanks = setup_tmpltbank_dax_generated(cp,scienceSegs,ahopeDax)
    
    return tmpltBanks

def setup_tmpltbank_dax_generated(cp,scienceSegs,ahopeDax):
    '''
    Setup template bank jobs that are generated as part of the ahope workflow.
    This function will add numerous jobs to the ahope workflow using
    configuration options from the .ini file. The following executables are
    currently supported:

    * lalapps_tmpltbank

    Parameters
    ----------
    cp : ConfigParser
        The ConfigParser holding all the options used by the ahope workflow.
    scienceSegs : Keyed dictionary of glue.segmentlist objects
        scienceSegs[ifo] holds the science segments to be analysed for each
        ifo. 
    ahopeDax : Instanced CondorDag class
        The CondorDag class that will hold all the jobs that the ahope workflow
        needs to run.

    Returns
    --------
    AhopeOutFileList
        The AhopeOutFileList holding the details of all the template bank jobs.
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

def setup_tmpltbank_pregenerated(cp,ifos):
    '''
    Setup ahope workflow to use a pregenerated template bank.
    The bank given in cp.get('ahope','pregenerated-template-bank') will be used
    as the input file for all matched-filtering jobs. If this option is
    present, ahope will assume that it should be used and not generate
    template banks within the workflow.

    Parameters
    ----------
    cp : ConfigParser
        The ConfigParser holding all the options used by the ahope workflow.
    ifos : list of strings
        The list of ifos that are used in this analysis.

    Returns
    --------
    AhopeOutFileList
        The AhopeOutFileList holding the details of the template bank.
    '''
    # Currently this uses the *same* fixed bank for all ifos.
    # Maybe we want to add capability to analyse separate banks in all ifos?
    
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
