from __future__ import division

import os
import urlparse, urllib
from glue import segments
from pycbc.ahope import sngl_ifo_job_setup, AhopeOutFileList
from pycbc.ahope import AhopeOutFile, select_tmpltbankjob_instance
from pycbc.ahope import select_matchedfilterjob_instance

def setup_tmpltbank_workflow(cp, scienceSegs, datafindOuts, ahopeDax,\
                             outputDir):
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
        tmpltBanks = setup_tmpltbank_pregenerated(cp, scienceSegs.keys())
    # Else we assume template banks will be generated in the workflow
    else:
        tmpltBanks = setup_tmpltbank_dax_generated(cp, scienceSegs, \
                       datafindOuts, ahopeDax, outputDir)
    
    return tmpltBanks

def setup_tmpltbank_dax_generated(cp, scienceSegs, datafindOuts, ahopeDax,\
                                  outputDir, link_to_matchedfltr=True):
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
    exeInstance = select_tmpltbankjob_instance(tmpltBankExe,'tmpltbank')

    if link_to_matchedfltr:
        # Use this to ensure that inspiral and tmpltbank jobs overlap. This
        # means that there will be 1 inspiral job for every 1 tmpltbank and
        # the data read in by both will overlap as much as possible. (If you
        # ask the template bank jobs to use 2000s of data for PSD estimation
        # and the matched-filter jobs to use 4000s, you will end up with
        # twice as many matched-filter jobs that still use 4000s to estimate a
        # PSD but then only generate triggers in the 2000s of data that the
        # template bank jobs ran on.
        tmpltbankExe = os.path.basename(cp.get('executables', 'inspiral'))
        linkExeInstance = select_matchedfilterjob_instance(tmpltbankExe, \
                                                           'inspiral')
    else:
        linkExeInstance = None

    # Set up class for holding the banks
    tmpltBanks = AhopeOutFileList([])

    # Template banks are independent for different ifos, but might not be!
    # Begin with independent case and add after FIXME
    for ifo in ifos:
        sngl_ifo_job_setup(cp, ifo, tmpltBanks, exeInstance, scienceSegs[ifo],\
                           datafindOuts, ahopeDax, outputDir, \
                           linkExeInstance=linkExeInstance, allowOverlap=True)

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
        userTag = "PREGEN_TMPLTBANK"
        fileUrl = urlparse.urljoin('file:', urllib.pathname2url(preGenBank))
        currFile = AhopeOutFile(ifo, userTag, globalSeg, fileUrl)
        tmpltBanks.append(currFile)

    return tmpltBanks
