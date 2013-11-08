import math
import numpy
from glue import lal
from glue import segments

class AhopeOutFile(lal.CacheEntry):
    '''This class holds the details of an individual output file in the ahope
    workflow. This file may be pre-supplied, generated from within the ahope
    command line script, or generated within the workflow. This class inherits
    from the glue.lal.CacheEntry class and has all attributes/methods of that
    class. It also adds some additional stuff for ahope. The important stuff
    from both is:

    * The location of the output file (which may not yet exist)
    * The ifo that the AhopeOutFile is valid for
    * The time span that the AhopeOutFile is valid for
    * A short description of what the file is
    * The dax job that will generate the output file (if appropriate). If the
      file is generated within the workflow the dax job object will hold all
      the job-specific information that may be relevant for later stages.

    An example of initiating this class:
    
    c = AhopeOutFile("H1", "INSPIRAL_S6LOWMASS", segments.segment(815901601, 815902177.5), "file://localhost/home/kipp/tmp/1/H1-815901601-576.xml", job=CondorDagNodeInstance)
    '''

    def __init__(self, ifo, description, timeSeg, fileUrl, job=None, **kwargs):
        lal.CacheEntry.__init__(self, ifo, description, timeSeg, fileUrl,\
                                **kwargs)
        self.job = job

class AhopeOutFileList(lal.Cache):
    '''This class holds a list of AhopeOutFile objects. It inherits from the
    built-in list class, but also allows a number of features. ONLY
    AhopeOutFile instances should be within a AhopeOutFileList instance.
    '''
    entry_class = AhopeOutFile

    def find_output(self,ifo,time):
        '''Return AhopeOutFile that covers the given time, or is most
        appropriate for the supplied time range.

        Parameters
        -----------
        ifo : string
           Name of the ifo that the 
        time : int/float/LIGOGPStime or tuple containing two values
           If int/float/LIGOGPStime (or similar may of specifying one time) is
           given, return the AhopeOutFile corresponding to the time. This calls
           self.find_output_at_time(ifo,time).
           If a tuple of two values is given, return the AhopeOutFile that is
           **most appropriate** for the time range given. This calls
           self.find_output_in_range

        Returns
        --------
        AhopeOutFile class
           The AhopeOutFile that corresponds to the time/time range
        '''
        # Determine whether I have a specific time, or a range of times
        try:
            lenTime = len(time)
        except TypeError:
            # This is if I have a single time
            outFile = self.find_output_at_time(ifo,time)                
        else:
            # This is if I have a range of times
            if lenTime == 2:
                outFile = self.find_output_in_range(ifo,time[0],time[1])
            # This is if I got a list that had more (or less) than 2 entries
            if len(time) != 2:
                raise TypeError("I do not understand the input variable time")
        return outFile

    def find_output_at_time(self,ifo,time):
       '''Return AhopeOutFile that covers the given time.

        Parameters
        -----------
        ifo : string
           Name of the ifo that the AhopeOutFile should correspond to
        time : int/float/LIGOGPStime
           Return the AhopeOutFiles that covers the supplied time. If no
           AhopeOutFile covers the time this will return None.

        Returns
        --------
        list of AhopeOutFile classes
           The AhopeOutFiles that corresponds to the time.
        '''
       # Get list of AhopeOutFiles that overlap time, for given ifo
       outFiles = [i for i in self if ifo == i.observatory \
                                   and time in i.segment] 
       if len(outFiles) == 0:
           # No AhopeOutFile at this time
           return None
       elif len(outFiles) == 1:
           # 1 AhopeOutFile at this time (good!)
           return outFiles
       else:
           # Multiple output files. Currently this is valid, but we may want
           # to demand exclusivity later, or in certain cases. Hence the
           # separation.
           return outFiles

    def find_output_in_range(self,ifo,start,end):
        '''Return the AhopeOutFile that is most appropriate for the supplied
        time range. That is, the AhopeOutFile whose coverage time has the
        largest overlap with the supplied time range. If no AhopeOutFiles
        overlap the supplied time window, will return None.

        Parameters
        -----------
        ifo : string
           Name of the ifo that the AhopeOutFile should correspond to
        start : int/float/LIGOGPStime 
           The start of the time range of interest.
        end : int/float/LIGOGPStime
           The end of the time range of interest

        Returns
        --------
        AhopeOutFile class
           The AhopeOutFile that is most appropriate for the time range
        '''
        # First filter AhopeOutFiles corresponding to ifo
        outFiles = [i for i in self if ifo == i.observatory] 
        if len(outFiles) == 0:
            # No AhopeOutFiles correspond to that ifo
            return None
        # Filter AhopeOutFiles to those overlapping the given window
        currSeg = segments.segment([start,end])
        outFiles = [i for i in outFiles if \
                               i.segment.intersects(currSeg)]
        if len(outFiles) == 0:
            # No AhopeOutFile overlap that time period
            return None
        elif len(outFiles) == 1:
            # One AhopeOutFile overlaps that period
            return outFiles[0]
        else:
            # More than one AhopeOutFile overlaps period. Find lengths of
            # overlap between time window and AhopeOutFile window
            overlapWindows = [abs(i.segment & currSeg) \
                                  for i in outFiles]
            # Return the AhopeOutFile with the biggest overlap
            # Note if two AhopeOutFile have identical overlap, this will return
            # the first AhopeOutFile in the list
            overlapWindows = numpy.array(overlapWindows,dtype = int)
            return outFiles[overlapWindows.argmax()]

    def find_all_output_in_range(self, ifo, currSeg):
        """
        Return all files that overlap the specified segment.
        """
        outFiles = [i for i in self if ifo == i.observatory]
        outFiles = [i for i in outFiles if \
                               i.segment.intersects(currSeg)]
        return self.__class__(outFiles)

class AhopeOutGroupList(AhopeOutFileList):
    '''
    This class holds a list of AhopeOutGroup objects. It inerits from the
    built-in list class, but also allows a number of additional features. ONLY
    AhopeOutGroup instances should be within a AhopeOutGroupList instance.

    NOTE: This class should not use some of the stuff in the underlying
    lal.Cache class
    '''
    # For now this is simply an AhopeOutFileList, but we may want to diverge
    # the two classes in the future.
    pass


class AhopeOutGroup:
    """
    This calls holds the details of a group of files that collectively cover
    a specified stretch of time in the ahope workflow. An example of where
    this might be used is if the template bank has been split up into a 
    number of components and analysed separately there will be a number of
    files (the split template banks and the split matched-filter output)
    that will cover that time region. This allows these files to be grouped
    together. This group may also have a single file that holds information
    about the group. An example of this is datafind, where a group of frame
    files covers the desired region, but codes generally read in a frame cache
    file which points to the locations of these. This class holds

    * The ifo that the AhopeOutGroup is valid for 
    * The time span that the AhopeOutGroup is valid for
    * A AhopeOutFileList storing the files and associated jobs in the group
    * A URL to the group summary file (for e.g. a frame-cache for datafind)
      if appropriate for the group (OPTIONAL)
    * A job that generated the URL group summary file, if that file exists
      and if it was generated from within the pegasus workflow
    """
    def __init__(self, ifo, description, time, summaryUrl=None, \
                 summaryJob=None):
        # Set everything to None to start
        self.__outFileList = None
        self.observatory = ifo
        self.segment = time
        self.description = description
        self.summaryUrl = summaryUrl
        self.summaryJob = summaryJob
        
    def get_output(self):
        '''Return self.__outFileList if it has been set, fail if not.

        Parameters
        ----------
        None

        Returns
        ----------
        AhopeOutFileList
            The AhopeOutFileList containing the list of files and jobs that
            will run them.
        '''
        if self.__outFileList:
            return self.__outFileList
        else:
            raise ValueError("Output file list has not been set.")

    def set_output(self, outFiles, outFileJobs):
        '''Set self.__outFile to outFile.

        Parameters
        ----------
        outFiles : list of strings
            URL paths to the list of output files
        outFileJobs : pipeline.CondorDagNode(s)
            This is the list of jobs used to generate the outFiles.
            If len(outFiles) == len(outFileJobs) then assume that these
            correspond so that outFileJobs[i] will generate outFiles[i].
            If len(outFiles) == 1 assume that one object generated all output.
            If outFileJobs == None then assume these files already exist.

        Returns
        ----------
        None
        '''
        if (not self.observatory) or (not self.segment):
            errMsg = "Ifo and time must be set before setting the output for " 
            errMsg += "AhopeOutGroup instances."
            raise ValueError(errMsg)

        # First some input conversion if needed
        if isinstance(outFiles, basestring):
            # If we got a string, make it a list of one string
            outFiles = [outFiles]

        outputList = AhopeOutFileList([])
        for i, fileUrl in enumerate(outFiles):
            if not outFileJobs:
                currJob = None
            elif len(outFileJobs) == 1:
                currJob = outFileJobs[0]
            elif len(outFileJobs) == len(outFiles):
                currJob = outFileJobs[i]
            else:
                errMsg = "The number of jobs given to .set_output must be "
                errMsg += "equal to the length of files or equal to 1."
                raise ValueError(errMsg)
            # Add the number to discriminate each job. This will be used later
            currDesc = self.description + "_%d" %(i)
            currFile = AhopeOutFile(self.observatory, self.description, \
                                    self.segment, fileUrl, job=currJob)
            outputList.append(currFile)
        self.__outFileList = outputList

    @property
    def summaryUrl(self):
        """
        The summary file's URL.  The URL is constructed from the
        values of the scheme, host, and path attributes.  Assigning
        a value to the URL attribute causes the value to be parsed
        and the scheme, host and path attributes updated.
        Stolen shamelessly from Kipp's glue.lal module.
        """
        return urlparse.urlunparse((self.summaryScheme, self.summaryHost, \
                                    self.summaryPath, None, None, None))

    @summaryUrl.setter
    def summaryUrl(self, url):
        self.scheme, self.host, self.path = urlparse.urlparse(url)[:3]

        
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
                                         jobValidSeg, parent=currParent)
                # Make the AhopeOutFile instance
                currFile = AhopeOutFile(ifo, jobTag, jobValidSeg, fileUrl,\
                                        job=currExeNode )
                outFiles.append(currFile)
    return outFiles

def split_outfiles(cp, inputFileList, exeInstance, numBanks, ahopeDax):
    """
    Add documentation
    """
    # Set up output structure
    outFileGroups = AhopeOutGroupList([])

    # Set up the condorJob class for the current executable
    currExeJob = exeInstance.create_condorjob(cp, None)

    for input in inputFileList:
        jobTag = input.description + "_" + exeInstance.exeName.upper()
        currExeNode, outUrlList = exeInstance.create_condornode(\
                                      ahopeDax, currExeJob, numBanks, input)
        outFileGroup = AhopeOutGroup(input.observatory, jobTag, input.segment)
        outFileGroup.set_output(outUrlList, [currExeNode])
        outFileGroups.append(outFileGroup)
    return outFileGroups
