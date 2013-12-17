import os, sys
import subprocess
import logging
import math
import numpy
import urlparse
from os.path import splitext, basename
from glue import lal
from glue import segments, pipeline
import pycbc.ahope as ahope
import pylal.dq.dqSegmentUtils as dqUtils

class Executable(object):
    """
    """
    def __init__(self):
        pass

class Job(object):
    def __init__(self):
        pass
    def create_node(self):
        pass

class Node(pipeline.CondorDADNode):
    def __init__(self, job):
        pipeline.CondorDAGNode.__init__(self, job)
        self.input_files = []
        self.output_files = []
        
    def add_input(self, files, opts=None):
        """files can be an AhopeFile or an AhopeFileGroup
        """
        for file in files:
            self.input_files.append(file)
            self.add_input_file(file.filename)
            self.add_parent(file.node)
        if opt:
            for file, opt in zip(files, opts):
                self.add_var_opt(opt, file.filename)
            
    def add_output(self, files, opts=None):
        for file in files:
            self.output_files.append(file)
            self.add_input_file(file.filename)
            self.add_parent(file.node)
        if opt:
            for file, opt in zip(files, opts):
                self.add_var_opt(opt, file.filename)
           
    def get_output(self):
        pass
        
    def get_input(self):
        pass
    

class Workflow(object):
    """This class manages an aHOPE style workflow. It provides convenience 
    functions for finding input files using time and keywords. It can also
    generate cache files from the inputs.
    """
    def __init__(self, config):
        """Create an aHOPE workflow
        """
        # Parse ini file
        self.cp = ahope.parse_ahope_ini_file(config)
        self.basename = basename(splitext(config)[0])
        
        # Initialize the dag
        logfile = self.basename + '.log'
        fh = open( logfile, "w" )
        fh.close()
        self.dag = pipeline.CondorDAG(logfile, dax=False)
        self.dag.set_dax_file(self.basename)
        self.dag.set_dag_file(self.basename)
        
    def add_files(self, files):
        """Add files to the workflow file collection. These can later be 
        used for input to nodes and may be queried for. File names must satisfy
        the lal cache name standard. 
        """
        pass
        
    def find_files(self, desc, time=None, ifo=None, **kwds):
        pass
        
    def add_node(self, node):
        self.dag.add_node(node)
        
    def write_plans():
        self.dag.write_sub_files()
        self.dag.write_dag()
        #self.dag.write_abstract_dag()
        self.dag.write_script()

class AhopeFile(lal.CacheEntry):
    '''This class holds the details of an individual output file in the ahope
    workflow. This file may be pre-supplied, generated from within the ahope
    command line script, or generated within the workflow. This class inherits
    from the glue.lal.CacheEntry class and has all attributes/methods of that
    class. It also adds some additional stuff for ahope. The important stuff
    from both is:

    * The location of the output file (which may not yet exist)
    * The ifo that the AhopeFile is valid for
    * The time span that the AhopeOutFile is valid for
    * A short description of what the file is
    * The dax node that will generate the output file (if appropriate). If the
      file is generated within the workflow the dax job object will hold all
      the job-specific information that may be relevant for later stages.

    An example of initiating this class:
    
    c = AhopeFile("H1", "INSPIRAL_S6LOWMASS", segments.segment(815901601, 815902177.5), "file://localhost/home/kipp/tmp/1/H1-815901601-576.xml", job=CondorDagNodeInstance)
    '''
    def __init__(self, description, extension, directory, ifo='N', 
                 time_seg=segments.segment(0, 99999999999), **kwargs):
        
        self.node=None
        self.kwargs = kwargs         
        self.filename = self._filename(ifo, description, extension, time_seg)
        self.path = os.path.join(directory, filename)
        file_url = urlparse.urlunparse(['file', 'localhost', self.path, None, None, None])

        lal.CacheEntry.__init__(self, ifo, description, time_seg, file_url)

        
    @classmethod
    def from_url(self, ifo, description, time_seg, file_url, node=None, **kwargs):
        pass
        
    def _filename(self, ifo, description, extension, time_seg, part=None):
        """ Construct the standard output filename
        """
        if part:
            description += '_' + str(part)
         
        duration = int(time_seg[1] - time_seg[0])
        start = time_seg[0]
        
        filename = ifo + '-' + description.upper() + '-' + start + '-' + duration + extension

class AhopeFileGroup(list):
    pass        


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

    def removeShortSciSegs(self, minSegLength):
        """
        Function to remove all science segments
        shorter than a specific length. Also updates the file on disk to remove
        these segments.

        Parameters
        -----------
        minSegLength : int
            Maximum length of science segments. Segments shorter than this will
            be removed.
        """
        newSegList = segments.segmentlist()
        for seg in self.segmentList:
            if abs(seg) > minSegLength:
                newSegList.append(seg)
        newSegList.coalesce()
        self.segmentlist = newSegList
        self.toSegmentXml()

    def toSegmentXml(self):
        """
        Write the segment list in self.segmentList to the url in self.url.
        """
        filePointer = open(self.path, 'w')
        dqUtils.tosegmentxml(filePointer, self.segmentList)
        filePointer.close()

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
              
