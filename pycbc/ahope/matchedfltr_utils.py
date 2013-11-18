from __future__ import division

import os
from pycbc.ahope import sngl_ifo_job_setup, AhopeOutFileList
from pycbc.ahope import select_matchedfilterjob_instance


def setup_matchedfltr_workflow(cp, scienceSegs, datafindOuts, ahopeDax,\
                               tmpltBanks, outputDir):
    '''
    Setup matched filter section of ahope workflow.
    FIXME: ADD MORE DOCUMENTATION
    '''

    # There should be a number of different options here, for e.g. to set
    # up fixed bank, or maybe something else
    inspiralOuts = setup_matchedfltr_dax_generated(cp, scienceSegs, \
                       datafindOuts, ahopeDax, tmpltBanks, outputDir)
    
    return inspiralOuts

def setup_matchedfltr_dax_generated(cp, scienceSegs, datafindOuts, ahopeDax,\
                                    tmpltBanks, outputDir,\
                                    link_to_tmpltbank=True):
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
    exeInstance = select_matchedfilterjob_instance(matchFltrExe,'inspiral')

    if link_to_tmpltbank:
        # Use this to ensure that inspiral and tmpltbank jobs overlap. This
        # means that there will be 1 inspiral job for every 1 tmpltbank and
        # the data read in by both will overlap as much as possible. (If you
        # ask the template bank jobs to use 2000s of data for PSD estimation
        # and the matched-filter jobs to use 4000s, you will end up with
        # twice as many matched-filter jobs that still use 4000s to estimate a
        # PSD but then only generate triggers in the 2000s of data that the
        # template bank jobs ran on.
        tmpltbankExe = os.path.basename(cp.get('executables', 'tmpltbank'))
        linkExeInstance = select_tmpltbankjob_instance(tmpltbankExe, \
                                                       'tmpltbank')
    else:
        linkExeInstance = None

    # Set up class for holding the banks
    inspiralOuts = AhopeOutFileList([])

    # Template banks are independent for different ifos, but might not be!
    # Begin with independent case and add after FIXME
    for ifo in ifos:
        sngl_ifo_job_setup(cp, ifo, inspiralOuts, exeInstance, \
                           scienceSegs[ifo], datafindOuts, ahopeDax, outputDir,\
                           parents=tmpltBanks, linkExeInstance=linkExeInstance,\
                           allowOverlap=False)

    return inspiralOuts

def select_matchedfilterjob_instance(currExe,currSection):
    """This function returns an instance of the class that is appropriate for
    matched-filtering within ahope.
    
    Parameters
    ----------
    currExe : string
        The name of the executable that is being used.
    currSection : string
        The name of the section storing options for this executble

    Returns
    --------
    Instanced class : exeClass
        An instance of the class that holds the utility functions appropriate
        for the given executable. This class **must** contain
        * exeClass.get_valid_times()
        * exeClass.create_condorjob()
        * exeClass.create_condornode()
    """

    # This is basically a list of if statements
    if currExe == 'lalapps_inspiral':
        exeClass = legacy_sngl_job_utils(currSection)
    # Some elif statements
    else:
        # Should we try some sort of default class??
        errString = "No class exists for executable %s" %(currExe,)
        raise NotImplementedError(errString)

    return exeClass
