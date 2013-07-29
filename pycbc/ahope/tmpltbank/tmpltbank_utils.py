from __future__ import division

import os
import math
import ConfigParser
from glue import segments
from glue import pipeline
from pycbc.ahope import AhopeOutFileList,AhopeOutFile
from pycbc.ahope.tmpltbank import select_tmpltbank_instance

def setup_tmplbank_workflow(cp,scienceSegs,ahopeDax):
    '''
    Setup template bank section of ahope workflow.
    FIXME: ADD MORE DOCUMENTATION
    '''

    # There should be a number of different options here, for e.g. to set
    # up fixed bank, or maybe something else
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

    tmpltBankExe = os.path.basename(cp.get('executables','tmpltbank'))
    # Select the appropriate class
    exeInstance = select_tmpltbank_instance(tmpltBankExe)

    # Set up class for holding the banks
    tmpltBanks = AhopeOutFileList([])

    # Template banks are independent for different ifos, but might not be!
    # Begin with independent case and add after FIXME
    ifos = scienceSegs.keys()

    for ifo in ifos:
        # Begin by getting analysis start and end, and start and end of time
        # that the bank is valid for
        dataLength,validChunk = exeInstance.get_tmpltbank_valid_times(cp,ifo)

        # Set up the tmpltbankJob class
        tmpltBankJob = exeInstance.create_tmpltbank_condorjob(cp,ifo)
        
        dataLoss = dataLength - abs(validChunk)
        if dataLoss < 0:
            raise ValueError("Ahope needs fixing! Please contact a developer")
        # Loop over science segments and set up banks
        for currSeg in scienceSegs[ifo]:
            # How many banks do we need
            currSegLength = abs(currSeg)
            numBanks = int( math.ceil( \
                       (currSegLength - dataLoss) / abs(validChunk) ))
            # What is the incremental shift between banks
            timeShift = (abs(currSeg) - dataLength) / (numBanks - 1)
            for bankNum in range(numBanks):
                # Get the science segment for this bank
                shiftDur = currSeg[0] + int(timeShift * bankNum)
                dataChunk = segments.segment([0,dataLength])
                bankDataSeg = dataChunk.shift(shiftDur)
                tmpltBankNode,bankFile = \
                        exeInstance.create_tmpltbank_condornode(ahopeDax,\
                                                      tmpltBankJob,bankDataSeg)
                # Make the TmpltBankList instance
                currBank = AhopeOutFile()
                currBank.set_bank(bankFile)
                currBank.set_ifo(ifo)
                currBank.set_time(segments.segmentlist([validChunk]))
                currBank.set_job(tmpltBankJob)
                tmpltBanks.append(currBank)

    return tmpltBanks
