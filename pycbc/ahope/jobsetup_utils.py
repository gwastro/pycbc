import math
from glue import segments
from pycbc.ahope import AhopeFile
from pycbc.ahope.legacy_ihope import *

def select_tmpltbankjob_instance(currExe, currSection):
    """This function returns an instance of the class that is appropriate for
    creating a template bank within ihope.
    
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
    if currExe == 'lalapps_tmpltbank':
        exeClass = LegacyTmpltbankExec(currSection)
    elif currExe == 'pycbc_geom_nonspin':
        exeClass = PyCBCTmpltbankExec(currSection)
    # Some elif statements
    else:
        # Should we try some sort of default class??
        errString = "No class exists for executable %s" %(currExe,)
        raise NotImplementedError(errString)

    return exeClass

def select_matchedfilterjob_instance(currExe, currSection):
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
        exeClass = LegacyInspiralExec(currSection)
    elif currExe == 'pycbc_inspiral':
        exeClass = PyCBCInspiralExec(currSection)
    # Some elif statements
    else:
        # Should we try some sort of default class??
        errString = "No class exists for executable %s" %(currExe,)
        raise NotImplementedError(errString)

    return exeClass

def select_splitfilejob_instance(currExe, currSection):
    """This function returns an instance of the class that is appropriate for
    splitting an output file up within ahope (for e.g. splitbank).
    
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
        * exeClass.create_condorjob()
        * exeClass.create_condornode()
    """

    # This is basically a list of if statements
    if currExe == 'lalapps_splitbank':
        exeClass = LegacySplitBankExec(currSection)
    # Some elif statements
    else:
        # Should we try some sort of default class??
        errString = "No class exists for executable %s" %(currExe,)
        raise NotImplementedError(errString)

    return exeClass

def sngl_ifo_job_setup(workflow, ifo, outFiles, exeInstance, scienceSegs, 
                       datafindOuts, outputDir,
                       parents=None, linkExeInstance=False, allowOverlap=True):
    """
    This function sets up a set of single ifo jobs.

    Parameters
    -----------
    cp : ConfigParser
        The ConfigParser object holding the parameters of the ahope workflow.
    ifo : string
        The name of the ifo to set up the jobs for
    outFiles : AhopeOutFileList or AhopeOutGroupList
        The AhopeOutFileList containing the list of jobs. Jobs will be appended
        to this list, and it does not need to be empty when supplied.
    exeInstance : Instanced class
        An instanced class that contains the functions needed to set up things
        that are specific to the executable being run.
    scienceSegs : segments.segmentlist
        The list of times that the jobs should cover
    ahopeDax : CondorDAG object
        The condorDAG object holding the ahope workflow being constructed.
    parents : AhopeOutFileList (optional, kwarg, default=None)
        The AhopeOutFileList containing the list of jobs that are parents to
        the one being set up.
    allowOverlap : boolean (optional, kwarg, default = True)
        If this is set the times that jobs are valid for will be allowed to
        overlap. This may be desired for template banks which may have some
        overlap in the times they cover. This may not be desired for inspiral
        jobs, where you probably want triggers recorded by jobs to not overlap
        at all.
    """
    cp = workflow.cp
    
    # Begin by getting analysis start and end, and start and end of time
    # that the output file is valid for
    dataLength, validChunk = exeInstance.get_valid_times(cp, ifo)

    dataChunk = segments.segment([0, dataLength])
    jobTag = exeInstance.exe_name.upper()
    
    if linkExeInstance:
        _, linkValidChunk = linkExeInstance.get_valid_times(cp, ifo)
        validChunkStart = max(validChunk[0], linkValidChunk[0])
        validChunkEnd = min(validChunk[1], linkValidChunk[1])
        validChunk = segments.segment([validChunkStart, validChunkEnd])


    # Set up the condorJob class for the current executable
    currExeJob = exeInstance.create_job(cp, ifo, outputDir)
    
    dataLoss = dataLength - abs(validChunk)
    
    if dataLoss < 0:
        raise ValueError("Ahope needs fixing! Please contact a developer")
    # Loop over science segments and set up jobs
    for currSeg in scienceSegs:
        # Is there enough data to analyse?
        currSegLength = abs(currSeg)
        if currSegLength < dataLength:
            continue
        # How many jobs do we need
        currSegLength = float(abs(currSeg))
        numJobs = int( math.ceil( \
                 (currSegLength - dataLoss) / float(abs(validChunk)) ))
        # What is the incremental shift between jobs
        timeShift = (currSegLength - dataLength) / float(numJobs - 1)
        for jobNum in range(numJobs):
            # Get the science segment for this job
            shiftDur = currSeg[0] + int(timeShift * jobNum)
            jobDataSeg = dataChunk.shift(shiftDur)
            jobValidSeg = validChunk.shift(shiftDur)
            # If we need to recalculate the valid times to avoid overlap
            if not allowOverlap:
                dataPerJob = (currSegLength - dataLoss) / float(numJobs)
                lowerBoundary = int(jobNum*dataPerJob \
                                    + validChunk[0] + currSeg[0])
                upperBoundary = int(dataPerJob + lowerBoundary)
                if lowerBoundary < jobValidSeg[0] or \
                        upperBoundary > jobValidSeg[1]:
                    errMsg = "Ahope is attempting to generate output "
                    errMsg += "from a job at times where it is not valid."
                    raise ValueError(errMsg)
                jobValidSeg = segments.segment([lowerBoundary, upperBoundary])
                
            if parents:
                currParent = parents.find_output(ifo, jobValidSeg)
                if not currParent:
                    errString = "No parent jobs found overlapping %d to %d." \
                                %(jobValidSeg[0],jobValidSeg[1])
                    errString += "\nThis is a bad error! Contact a developer."
                    raise ValueError(errString)
            else:
                currParent = None

            if datafindOuts:
                currDfOuts = datafindOuts.find_all_output_in_range(ifo, 
                                                                   jobDataSeg)
                if not currDfOuts:
                    errString = "No datafind jobs found overlapping %d to %d."\
                                %(jobDataSeg[0],jobDataSeg[1])
                    errString += "\nThis shouldn't happen. Contact a developer."
                    raise ValueError(errString)

            currExeNode = currExeJob.create_node(jobDataSeg,
                                     jobValidSeg, parent=currParent,
                                     dfParents=currDfOuts)
            workflow.add_node(currExeNode)
            outFiles += currExeNode.output_files
    return outFiles

