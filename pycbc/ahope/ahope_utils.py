import os, sys
import subprocess
import logging
import math
import numpy
import urlparse
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

class AhopeOutSegFile(AhopeOutFile):
    '''
    This class inherits from the AhopeOutFile class, and is designed to store
    ahope output files containing a segment list. This is identical in
    usage to AhopeOutFile except for an additional kwarg for holding the
    segment list, if it is known at ahope run time.
    '''
    def __init__(self, ifo, description, timeSeg, fileUrl,\
                 segList=None, **kwargs):
        """
        ADD DOCUMENTATION
        """
        AhopeOutFile.__init__(self, ifo, description, timeSeg, fileUrl,\
                              **kwargs)
        self.segmentList = segList

class AhopeOutFileList(lal.Cache):
    '''This class holds a list of AhopeOutFile objects. It inherits from the
    built-in list class, but also allows a number of features. ONLY
    AhopeOutFile instances should be within a AhopeOutFileList instance.
    '''
    entry_class = AhopeOutFile

    def find_output(self,ifo,time):
        '''
        Return AhopeOutFile that covers the given time, or is most
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


class AhopeOutGroup(object):
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
        if self.__outFileList is not None:
            return self.__outFileList
        else:
            raise ValueError("Output file list has not been set.")

    def set_output(self, outFiles, outFileJobs, outSegs=None):
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
        outSegs : List of glue.segments.segment objects
            If given len(outSegs) must equal len(outFiles). This is the
            segment that each individual outFile covers. If not given then
            assume that the segment for each job is self.segment.

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
            if outSegs:
                currSeg = outSegs[i]
            else:
                currSeg = self.segment
            # Add the number to discriminate each job. This will be used later
            currDesc = self.description + "_%d" %(i)
            currFile = AhopeOutFile(self.observatory, self.description, \
                                    currSeg, fileUrl, job=currJob)
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
        if self.summaryScheme and self.summaryHost and self.summaryPath:
            return urlparse.urlunparse((self.summaryScheme, self.summaryHost, \
                                    self.summaryPath, None, None, None))
        else:
            return None

    @summaryUrl.setter
    def summaryUrl(self, summaryUrl):
        """
        WRITE THIS
        """
        if summaryUrl is None:
            self.summaryScheme = self.summaryHost = self.summaryPath = None
        else:
            self.summaryScheme, self.summaryHost, self.summaryPath = \
                                             urlparse.urlparse(summaryUrl)[:3]


def make_external_call(cmdList, outDir=None, outBaseName='external_call',\
                       shell=False, fail_on_error=True):
    """
    Use this to make an external call using the python subprocess module.
    See the subprocess documentation for more details of how this works.
    http://docs.python.org/2/library/subprocess.html

    Parameters
    -----------
    cmdList : list of strings
        This list of strings contains the command to be run. See the subprocess
        documentation for more details.
    outDir : string
        If given the stdout and stderr will be redirected to
        os.path.join(outDir,outBaseName+[".err",".out])
        If not given the stdout and stderr will not be recorded
    outBaseName : string
        The value of outBaseName used to construct the file names used to
        store stderr and stdout. See outDir for more information.
    shell : boolean, default=False
        This value will be given as the shell kwarg to the subprocess call.
        **WARNING** See the subprocess documentation for details on this
        Kwarg including a warning about a serious security exploit. Do not
        use this unless you are sure it is necessary **and** safe.
    fail_on_error : boolean, default=True
        If set to true an exception will be raised if the external command does
        not return a code of 0. If set to false such failures will be ignored.
        Stderr and Stdout can be stored in either case using the outDir
        and outBaseName options.

    Returns
    --------
    exitCode : int
        The code returned by the process.
    """
    if outDir:
        outBase = os.path.join(outDir,outBaseName)
        errFile = outBase + '.err'
        errFP = open(errFile, 'w')
        outFile = outBase + '.out'
        outFP = open(outFile, 'w')
        cmdFile = outBase + '.sh'
        cmdFP = open(cmdFile, 'w')
        cmdFP.write(' '.join(cmdList))
        cmdFP.close()
    else:
        errFile = None
        outFile = None
        cmdFile = None
        errFP = None
        outFP = None

    msg = "Making external call %s" %(' '.join(cmdList))
    logging.debug(msg)
    errCode = subprocess.call(cmdList, stderr=errFP, stdout=outFP,\
                              shell=shell)
    if errFP:
        errFP.close()
    if outFP:
        outFP.close()

    if errCode and fail_on_error:
        raise CalledProcessErrorMod(errCode, ' '.join(cmdList), \
                errFile=errFile, outFile=outFile, cmdFile=cmdFile)
    logging.debug("Call successful, or error checking disabled.")

class CalledProcessErrorMod(Exception):
    """
    This exception is raised when subprocess.call returns a non-zero exit code
    and checking has been requested
    """
    def __init__(self, returncode, cmd, errFile=None, outFile=None, \
                 cmdFile=None):
        self.returncode = returncode
        self.cmd = cmd
        self.errFile = errFile
        self.outFile = outFile
        self.cmdFile = cmdFile
    def __str__(self):
        msg = "Command '%s' returned non-zero exit status %d.\n" \
              %(self.cmd, self.returncode)
        if self.errFile:
            msg += "Stderr can be found in %s.\n" %(self.errFile)
        if self.outFile:
            msg += "Stdout can be found in %s.\n" %(self.outFile)
        if self.cmdFile:
            msg += "The failed command has been printed in %s." %(self.cmdFile)
        return msg
              
