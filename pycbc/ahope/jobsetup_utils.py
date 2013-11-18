from glue import segments
from pycbc.ahope import AhopeOutGroup, AhopeOutFile

def sngl_ifo_job_setup(cp, ifo, outFiles, exeInstance, scienceSegs, \
                       datafindOuts, ahopeDax, outputDir,\
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
    # Begin by getting analysis start and end, and start and end of time
    # that the output file is valid for
    dataLength,validChunk = exeInstance.get_valid_times(cp, ifo)
    dataChunk = segments.segment([0, dataLength])
    jobTag = exeInstance.exeName.upper()

    if linkExeInstance:
        _, linkValidChunk = linkExeInstance.get_valid_times(cp, ifo)
        validChunkStart = max(validChunk[0], linkValidChunk[0])
        validChunkEnd = min(validChunk[1], linkValidChunk[1])
        validChunk = segments.segment([validChunkStart, validChunkEnd])


    # Set up the condorJob class for the current executable
    currExeJob = exeInstance.create_condorjob(cp, ifo, outputDir)

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
        currSegLength = abs(currSeg)
        numJobs = int( math.ceil( \
                 (currSegLength - dataLoss) / float(abs(validChunk)) ))
        # What is the incremental shift between jobs
        timeShift = (abs(currSeg) - dataLength) / float(numJobs - 1)
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
            # Get the parent job if necessary
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
                currDfOuts = datafindOuts.find_all_output_in_range(ifo, \
                                                                   jobDataSeg)
                if not currDfOuts:
                    errString = "No datafind jobs found overlapping %d to %d."\
                                %(jobDataSeg[0],jobDataSeg[1])
                    errString += "\nThis shouldn't happen. Contact a developer."
                    raise ValueError(errString)

            # If the parent produces a group of output files, such as
            # lalapps_splitbank, a number of condor jobs are needed
            if currParent.__class__.__name__ == 'AhopeOutGroup':
                # Set up the global outputs
                currFiles = AhopeOutGroup(ifo, jobTag, jobValidSeg)
                nodeList = []
                urlList = []
                for parentJob in currParent.get_output():
                    currExeNode, fileUrl = exeInstance.create_condornode(\
                                     ahopeDax, currExeJob, jobDataSeg,\
                                     jobValidSeg, parent=parentJob,\
                                     dfParents=currDfOuts)
                    nodeList.append(currExeNode)
                    urlList.append(fileUrl)
                currFiles.set_output(urlList, nodeList)
                outFiles.append(currFiles)
            else:
                currExeNode, fileUrl = exeInstance.create_condornode(\
                                         ahopeDax, currExeJob, jobDataSeg,\
                                         jobValidSeg, parent=currParent,\
                                         dfParents=currDfOuts)
                # Make the AhopeOutFile instance
                currFile = AhopeOutFile(ifo, jobTag, jobValidSeg, fileUrl,\
                                        job=currExeNode )
                outFiles.append(currFile)
    return outFiles

