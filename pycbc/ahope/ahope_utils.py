import math
import numpy
from glue import segments

class AhopeOutFileList(list):
    '''This class holds a list of AhopeOutFile objects. It inherits from the
    built-in list class, but also allows a number of features. ONLY
    AhopeOutFile instances should be within a AhopeOutFileList instance.
    '''

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
            # This is if I have a range of times
            if len(time) == 2:
                outFile = self.find_output_in_range(ifo,time[0],time[1])
        except TypeError:
            # This is if I have a single time
            outFile = self.find_output_at_time(ifo,time)                
        else:
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
       outFiles = [i for i in self if ifo == i.get_ifo() \
                                   and time in i.get_time()] 
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
        outFiles = [i for i in self if ifo == i.get_ifo()] 
        if len(outFiles) == 0:
            # No AhopeOutFiles correspond to that ifo
            return None
        # Filter AhopeOutFiles to those overlapping the given window
        currSeg = segments.segment([start,end])
        currSegList = segments.segmentlist([currSeg])
        outFiles = [i for i in outFiles if \
                               i.get_time().intersects(currSegList)]
        if len(outFiles) == 0:
            # No AhopeOutFile overlap that time period
            return None
        elif len(outFiles) == 1:
            # One AhopeOutFile overlaps that period
            return outFiles[0]
        else:
            # More than one AhopeOutFile overlaps period. Find lengths of
            # overlap between time window and AhopeOutFile window
            overlapWindows = [abs(i.get_time() & currSegList) \
                                  for i in outFiles]
            # Return the AhopeOutFile with the biggest overlap
            # Note if two AhopeOutFile have identical overlap, this will return
            # the first AhopeOutFile in the list
            overlapWindows = numpy.array(overlapWindows,dtype = int)
            return outFiles[overlapWindows.argmax()]

class AhopeOutGroupList(AhopeOutFileList):
    '''
    This class holds a list of AhopeOutGroup objects. It inerits from the
    built-in list class, but also allows a number of additional features. ONLY
    AhopeOutFile instances should be within a AhopeOutFileList instance.
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
    together. This class holds

    * The ifo that the AhopeOutGroup is valid for 
    * The time span that the AhopeOutGroup is valid for
    * A AhopeOutFileList storing the files and associated jobs in the group
    """
    def __init__(self, outFileList=None, ifo=None, time=None):
        # Set everything to None to start
        self.__outFileList = None
        self.__ifo = None
        self.__time = None
        
        # Parse and set kwargs if given
        if ifo:
            self.set_ifo(ifo)
        if time:
            self.set_time(time)
        if outFileList:
            self.set_output(outFileList)

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
            Path to the list of output files
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
        if (not self.__ifo) or (not self.__time):
            errMsg = "Ifo and time must be set before setting the output for " 
            errMsg += "AhopeOutGroup instances."
            raise ValueError(errMsg)

        # First some input conversion if needed
        if isinstance(outFiles, basestring):
            # If we got a string, make it a list of one string
            outFiles = [outFiles]

        outputList = AhopeOutFileList([])
        for i, filePath in enumerate(outFiles):
            if len(outFileJobs) == 1:
                currJob = outFileJobs[0]
            elif len(outFileJobs) == len(outFiles):
                currJob = outFileJobs[i]
            else:
                errMsg = "The number of jobs given to .set_output must be "
                errMsg += "equal to the length of files or equal to 1."
                raise ValueError(errMsg)
            currFile = AhopeOutFile(outFile=filePath ,ifo=self.get_ifo(),\
                       time=self.get_time(), job=currJob)
            outputList.append(currFile)
        self.__outFileList = outputList
        

    def get_ifo(self):
        '''Return self.__ifo if it has been set, Fail if not.

        Parameters
        ----------
        None

        Returns
        ----------
        ifo : string
           The ifo that the AhopeOutGroup is intended to be used
           for. If an output (such as a pregenerated file) is intended for
           **all** ifos then you should create multiple
           AhopeOutGroup classes, each with a different
           ifo. Will fail if the ifo has not yet been set.
        '''
        if self.__ifo:
            return self.__ifo
        else:
            raise ValueError("ifo has not been set for this instance.")

    def set_ifo(self,ifo):
        '''Set self.__ifo to ifo.

        Parameters
        ----------
        ifo : string
           Sets the ifo used for this AhopeOutGroup

        Returns
        ----------
        None
        '''
        self.__ifo = ifo

    def get_time(self):
        '''Return self.__time if it has been set, fail if not.

        Parameters
        ----------
        None

        Returns
        ----------
        time : glue.ligolw segmentlist
           This is the time for
           which the AhopeOutGroup should be valid.
           The output files may be generated using data that covers
           part of this time, or even a completely
           different time. This simply says
           that during the times given in the segmentlist,
           this output file is the preferred one to use.
           Will fail if the time has not yet been set.
        '''
        if self.__time:
            return self.__time
        else:
            raise ValueError("ifo has not been set for this instance.")

    def set_time(self,time):
        '''Set self.__time to time.

        Parameters
        ----------
        time : glue.ligolw segmentlist
           Sets the time for which the output file is valid.

        Returns
        ----------
        None
        '''
        if type(time) != segments.segmentlist:
            raise TypeError("Variable supplied to AhopeOutGroup.set_time() "+\
                            "must be a glue.segment.segmentlist class.")
        self.__time = time


class AhopeOutFile:
    '''This class holds the details of an individual output file in the ahope
    workflow. This file may be pre-supplied, generated from within the ahope
    command line script, or generated within the workflow. This stores 4 pieces
    of information, which will be used in later stages of the workflow. The 4
    values this class holds are:

    * The location of the output file (which may not yet exist)
    * The ifo that the AhopeOutFile is valid for
    * The time span that the AhopeOutFile is valid for
    * The dax job that will generate the output file (if appropriate). If the
      file is generated within the workflow the dax job object will hold all
      the job-specific information that may be relevant for later stages.
    '''

    def __init__(self,outFile=None,ifo=None,time=None,job=None):
        # Set everything to None to start
        self.__outFile = None
        self.__ifo = None
        self.__time = None
        self.__job = None
        
        # Parse and set kwargs if given
        if outFile:
            self.set_output(outFile)
        if ifo:
            self.set_ifo(ifo)
        if time:
            self.set_time(time)
        if job:
            self.set_job(job)

    def get_output(self):
        '''Return self.__outFile if it has been set, fail if not.

        Parameters
        ----------
        None

        Returns
        ----------
        outFile : string
          The location of the output file, this may not
          exist yet if the output file will be
          generated by the workflow. In that case
          this gives the location that the output file will be written to.
          Will fail if the output file has not yet been set.
        '''
        if self.__outFile:
            return self.__outFile
        else:
            raise ValueError("Bank has not been set for this instance.")

    def set_output(self,outFile):
        '''Set self.__outFile to outFile.

        Parameters
        ----------
        outFile : string
           The location of the output file, this may not
           exist yet if the outFile will be generated by the workflow.

        Returns
        ----------
        None
        '''
        self.__outFile = outFile

    def get_ifo(self):
        '''Return self.__ifo if it has been set, Fail if not.

        Parameters
        ----------
        None

        Returns
        ----------
        ifo : string
           The ifo that the outFile is intended to be used
           for. If an output (such as a pregenerated file) is intended for
           **all** ifos then you should create multiple
           AhopeOutFile classes, each with a different
           ifo. Will fail if the ifo has not yet been set.
        '''
        if self.__ifo:
            return self.__ifo
        else:
            raise ValueError("ifo has not been set for this instance.")

    def set_ifo(self,ifo):
        '''Set self.__ifo to ifo.

        Parameters
        ----------
        ifo : string
           Sets the ifo used for this output file

        Returns
        ----------
        None
        '''
        self.__ifo = ifo

    def get_time(self):
        '''Return self.__time if it has been set, fail if not.

        Parameters
        ----------
        None

        Returns
        ----------
        time : glue.ligolw segmentlist
           This is the time for
           which the output file should be valid.
           The output file may be generated using data that covers
           part of this time, or even a completely
           different time. This simply says
           that during the times given in the segmentlist,
           this output file is the preferred one to use.
           Will fail if the time has not yet been set.
        '''
        if self.__time:
            return self.__time
        else:
            raise ValueError("ifo has not been set for this instance.")

    def set_time(self,time):
        '''Set self.__time to time.

        Parameters
        ----------
        time : glue.ligolw segmentlist
           Sets the time for which the output file is valid.

        Returns
        ----------
        None
        '''
        if type(time) != segments.segmentlist:
            raise TypeError("Variable supplied to AhopeOutFile.set_time() "+\
                            "must be a glue.segment.segmentlist class.")
        self.__time = time

    def get_job(self):
        '''Return self.__job if it has been set, Fail if not.

        Parameters
        ----------
        None

        Returns
        ----------
        job : Condor job
           This is the the condor job that will generate the output file.
           If no job has been set, this
           will return None. For the case where
           the output file will not be generated from
           within the workflow, this should not be set.
        '''
        return self.__job

    def set_job(self,job):
        '''Set self.__job to job.

        Parameters
        ----------
        job : Condor job
           Sets the condor job that will generate the output file.

        Returns
        ----------
        None

        '''
        # FIXME: What sanity check is appropriate here?
        self.__job = job

def sngl_ifo_job_setup(cp, ifo, outFiles, exeInstance, scienceSegs, ahopeDax,\
                       parents=None, linkExeInstance=False):
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
    """
    # Begin by getting analysis start and end, and start and end of time
    # that the output file is valid for
    dataLength,validChunk = exeInstance.get_valid_times(cp, ifo)
    dataChunk = segments.segment([0, dataLength])
    print validChunk, dataLength

    if linkExeInstance:
        _, linkValidChunk = linkExeInstance.get_valid_times(cp, ifo)
        validChunkStart = max(validChunk[0], linkValidChunk[0])
        validChunkEnd = min(validChunk[1], linkValidChunk[1])
        validChunk = segments.segment([validChunkStart, validChunkEnd])
    

    # Set up the condorJob class for the current executable
    currExeJob = exeInstance.create_condorjob(cp,ifo)

    dataLoss = dataLength - abs(validChunk)
    print validChunk, dataLength
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
            if 0:
                dataPerJob = (currSegLength - dataLoss) / float(numJobs)
                lowerBoundary = jobNum*dataPerJob + validChunk[0]
                upperBoundary = dataPerJob + lowerBoundary
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

            # If the parent produces a group of output files, such as
            # lalapps_splitbank, a number of condor jobs are needed
            if currParent.__class__.__name__ == 'AhopeOutGroup':
                # Set up the global outputs
                currFiles = AhopeOutGroup()
                currFiles.set_ifo(ifo)
                currFiles.set_time(segments.segmentlist([jobValidSeg]))
                nodeList = []
                fileList = []
                for parentJob in currParent.get_output():
                    currExeNode, outFile = exeInstance.create_condornode(\
                                     ahopeDax, currExeJob, jobDataSeg,\
                                     parent=parentJob)
                    nodeList.append(currExeNode)
                    fileList.append(outFile)
                currFiles.set_output(fileList, nodeList) 
                outFiles.append(currFiles)
            else:
                currExeNode,outFile = exeInstance.create_condornode(\
                                         ahopeDax, currExeJob, jobDataSeg,\
                                         parent=currParent)
                # Make the AhopeOutFile instance
                currFile = AhopeOutFile()
                currFile.set_output(outFile)
                currFile.set_ifo(ifo)
                currFile.set_time(segments.segmentlist([jobValidSeg]))
                currFile.set_job(currExeNode)
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
        currExeNode, outFileList = exeInstance.create_condornode(\
                                      ahopeDax, currExeJob, numBanks, input)
        outFileGroup = AhopeOutGroup()
        outFileGroup.set_time(input.get_time())
        outFileGroup.set_ifo(input.get_ifo())
        outFileGroup.set_output(outFileList, [currExeNode])
        outFileGroups.append(outFileGroup)
    return outFileGroups
