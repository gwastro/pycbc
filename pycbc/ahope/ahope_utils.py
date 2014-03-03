# Copyright (C) 2013  Ian Harry, Alex Nitz
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""
This module provides the worker functions and classes that are used when
creating an ahope workflow. For details about ahope see here:
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/ahope.html
"""

import os, sys, subprocess, logging, math
import numpy
import urlparse
from itertools import combinations
from os.path import splitext, basename, isfile
from glue import lal
from glue import segments, pipeline
from configparserutils import AhopeConfigParser
import pylal.dq.dqSegmentUtils as dqUtils
import copy

#REMOVE THESE FUNCTIONS  FOR PYTHON >= 2.7 ####################################
def check_output(*popenargs, **kwargs):
    """
    This function is used to obtain the stdout of a command. It is only used
    internally, recommend using the make_external_call command if you want
    to call external executables.
    """
    if 'stdout' in kwargs:
        raise ValueError('stdout argument not allowed, it will be overridden.')
    process = subprocess.Popen(stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               *popenargs, **kwargs)
    output, unused_err = process.communicate()
    retcode = process.poll()
    return output

###############################################################################

def make_analysis_dir(path):
    """
    Make the analysis directory path, any parent directories that don't already
    exist, and the 'logs' subdirectory of path.
    """
    if path is not None:
        makedir(os.path.join(path, 'logs'))

def makedir(path):
    """
    Make the analysis directory path and any parent directories that don't
    already exist. Will do nothing if path already exists.
    """
    if not os.path.exists(path) and path is not None:
        os.makedirs(path)

def is_condor_exec(exe_path):
    """
    Determine if an executable is condor-compiled

    Parameters
    ----------
    exe_path : str
          The executable path

    Returns
    -------
    truth_value  : boolean
        Return True if the exe is condor compiled, False otherwise.
    """
    if check_output(['nm', '-a', exe_path]).find('condor') != -1:
        return True
    else:
        return False

class Job(pipeline.AnalysisJob, pipeline.CondorDAGJob):
    """
    The pyCBC Job class is an extension of the pipeline.AnalysisJob class.
    From the Dagman perspective a single AnalysisJob instance corresponds to
    a single condor submit file, which may be called by several condor nodes.
    In the DAX perspective it can be used to add several nodes to a workflow
    which have a number of common arguments.

    The specific pyCBC implementation adds a number of helper functions as
    documented below and allows for some extra functionality when initializing
    the instance, including checking if the executable is condor_compiled at
    run time and setting the stderr, stdout and log files to ahope standards.
    """

    def __init__(self, cp, exe_name, universe=None, ifo=None, 
                                     out_dir=None, tags=[]):
        """
        Initialize the pycbc.ahope.Job class.
   
        Parameters
        -----------
        cp : ConfigParser object
            The ConfigParser object holding the ahope configuration settings
        exec_name : string
            Executable name
        universe : string, optional
            Condor universe to run the job in
        ifo: string, optional
        out_dir: path, optional
            The folder to store output files of this job. 
        tags : list of strings
            A list of strings that is used to identify this job.
        """
        self.exe_name = exe_name
        self.cp = cp
        self.ifo = ifo
        
        tags = [tag.upper() for tag in tags]
        self.tags = tags
        
        # Determine the sub file name      
        if self.ifo:
            tags = tags + [self.ifo.upper()]
        self.tag_desc = '_'.join(tags)
        
        if len(tags) != 0:
            self.tagged_name = "%s-%s" % (exe_name, self.tag_desc)
        else:
            self.tagged_name = self.exe_name
        
        # Determine the output directory
        if out_dir is not None:
            self.out_dir = out_dir
        elif len(tags) == 0:
            self.out_dir = self.exe_name
        else:
            self.out_dir = self.tagged_name
            
        if not os.path.isabs(self.out_dir):
            self.out_dir = os.path.join(os.getcwd(), self.out_dir) 
            
        self.base_dir = os.path.basename(self.out_dir)
        
        exe_path = cp.get('executables', exe_name)
        
        # Check that the executable actually exists
        if os.path.isfile(exe_path):
	    logging.debug("Using %s executable "
                          "at %s" % (exe_name, exe_path))
        else:
            raise TypeError("Failed to find %s executable " 
                            "at %s" % (exe_name, exe_path))

        # Determine the condor universe if we aren't given one 
        # Default is vanilla unless the executable is condor-compiled
        # in which case we detect this and use standard universe
        if universe is None:
            if is_condor_exec(exe_path):
                universe = 'standard'
            else:
                universe = 'vanilla'

        logging.debug("%s executable will run as %s universe"
                     % (exe_name, universe))
        
        pipeline.CondorDAGJob.__init__(self, universe, exe_path)
        pipeline.AnalysisJob.__init__(self, cp, dax=True)       
        
        if not (universe == 'standard'):
            self.add_condor_cmd('getenv', 'True')
        self.add_condor_cmd('copy_to_spool','False')
        
        sections = [self.exe_name]
        for tag in tags:
             section = '%s-%s' %(self.exe_name, tag.lower())
             if cp.has_section(section):
                sections.append(section)
             
        # Do some basic sanity checking on the options
        for sec1, sec2 in combinations(sections, 2):
            cp.check_duplicate_options(sec1, sec2, raise_error=True)
             
        for sec in sections:
            if cp.has_section(sec):
                self.add_ini_opts(cp, sec)
            else:
                warnString = "warning: config file is missing section [%s]"\
                             %(sec,)
                logging.warn(warnString)

        # What would be a better more general logname ?
        logBaseNam = 'logs/%s-$(macrogpsstarttime)' %(exe_name,)
        logBaseNam += '-$(macrogpsendtime)-$(cluster)-$(process)'
        
        self.add_condor_cmd("initialdir", self.out_dir)
        self.set_stdout_file('%s.out' % (logBaseNam,) )
        self.set_stderr_file('%s.err' % (logBaseNam,) )
        self.set_sub_file('%s.sub' % (self.tagged_name) )
          
        # Set default requirements for a JOB, these can be changed
        # through methods of this class
        self.set_memory(1000)
        self.set_storage(100)
        
        # Make sure that the directories for this jobs output exist
        makedir(self.out_dir)
        makedir(os.path.join(self.out_dir, 'logs'))
    
    def set_memory(self, ram_value):
        """
        Set the amount of the RAM that this job requires.
        
        Parameters
        ----------
        ram_value: int
              The amount of ram that this job requires in MB.
        """
        self.add_condor_cmd('request_memory', '%d' %(ram_value))
        
    def set_storage(self, storage_value):
        """
        Set the amount of harddrive storage that this job requires.
        
        Parameters
        ----------
        ram_value: int
              The amount of storage that this job requires in MB.
        """
        self.add_condor_cmd('request_disk', '%dM' %(storage_value))
    
    def needs_gpu(self):
        """
        Call this function to indicate that this job wants to run on a GPU.
        FIXME: THIS IS NOT PROPERLY SUPPORTED YET!
        """
        # Satisfy the requirements to use GPUs on cit, sugar
        # FIXME add in the ATLAS support
        self.add_condor_cmd('+WantsGPU', 'true')
        self.add_condor_cmd('+WantGPU', 'true')
        self.add_condor_cmd('+HAS_GPU', 'true')
        self.add_condor_cmd('Requirements', '( GPU_PRESENT =?= true) || (HasGPU =?= "gtx580") || (Target.HAS_GPU =?= True)')

    def create_node(self):
        """
        Create a condor node from this job. This provides a basic interface to
        the Node class. Most jobs in an ahope workflow will subclass the 
        pycbc.ahope.Job class and overwrite this to give more details when
        initializing the node.
        """
        return Node(self)
        

class Node(pipeline.CondorDAGNode):
    """
    The pyCBC Node class is an extension of the pipeline.CondorDAGNode class.
    A single Node instance corresponds to a single command that will be run
    within the workflow. *Every* command/node in the workflow will have a Node
    instance.

    The specific pyCBC implementation adds a number of helper functions as
    documented below and allows for some extra functionality when initializing
    the instance, including allowing us to deal with cases where a set of
    partitioned input files are given. For e.g. if I have run lalapps_splitbank
    I now have a set of files corresponding to *one* template bank. With this
    function I can just feed the Node that AhopeFile pointing to the complete
    template bank (ie. the list of files making up that template bank) and this
    will automatically create a job for each file in the template bank.
    """
    def __init__(self, job):
        """
        Initialize the pycbc.ahope.Node class. This is often overridden in
        subclassed instances.
        
        Parameters
        -----------
        Job : pycbc.ahope.Job instance
            The pycbc.ahope.Job instance to create a Node from
        """

        pipeline.CondorDAGNode.__init__(self, job)
        self.input_files = AhopeFileList([])
        self.output_files = AhopeFileList([])
        self.set_category(job.exe_name)
        self.executed = False
        
    def add_input(self, file, opt=None, argument=False, recombine=False):
        """
        Add a file, or partitioned file, as input to this node. 
        
        Parameters
        ----------
        file : AhopeFile
            The AhopeFile that this node needs to run
        opt : string, optional
            The command line option that the executable needs
            in order to set the associated file as input.
        argument : Boolean, optional
            If present this indicates that the file should be supplied as an
            argument to the job (ie. file name with no option at the end of
            the command line call). Using argument=True and opt != None will
            result in a failure.
        recombine : Boolean, optional
            If present this indicates that the partitioned input file will be
            recombined in this Node. So only one Node is needed and all the
            partitioned input files will be parents. It is possible to
            sequentially add partitioned files to be recombined and other 
            partitioned files that are *not* to be recombined. 
        """
        if argument and opt != None:
            errMsg = "You cannot supply an option and tell the code that this "
            errMsg += "is an argument. Choose one or the other."
            raise ValueError(errMsg)
        
        # Add the files to the nodes internal lists of input
        self.input_files.append(file) 
        self.add_input_file(file.path)  
                                
        # If the file was created by another node, then make that
        # node a parent of this one
        if file.node and file.node.executed is False:
            self.add_parent(file.node)

        if opt:
            self.add_var_opt(opt, file.path)

        if argument:
            self.add_var_arg(file.path)

    def add_input_list(self, fileList, opt=None, argument=False, delimiter=' '):
        """
        Add a list of files as input to this node, under a single option, ie:
        --option-name file1,file2,file3, or as an argument. Use delimiter to
        specify how to separate the files.

        Parameters
        -----------
        fileList: AhopeFileList
            The list of AhopeFiles that this node needs to run.
        opt : string, optional
            The command line option that the executable needs
            in order to set the associated file list as input.
        argument : Boolean, optional
            If present this indicates that the files should be supplied as an
            argument to the job (ie. file names with no option at the end of
            the command line call). Using argument=True and opt != None will
            result in a failure.
        delimiter : string, optional (default = ' ')
            Set the delimiter that is used to separate the file names when
            supplied as an option or argument.
        """
        if argument and opt != None:
            errMsg = "You cannot supply an option and tell the code that this "
            errMsg += "is an argument. Choose one or the other."
            raise ValueError(errMsg)

        # Add the files to the nodes internal lists of input
        for file in fileList:
            self.input_files.append(file)
            self.add_input_file(file.path)
            if file.node and not hasattr(file.node, 'executed'):
                self.add_parent(file.node)

        fileListString = delimiter.join([file.path for file in fileList])
     
        if opt:
            self.add_var_opt(opt, fileListString)

        if argument:
            self.add_var_arg(fileListString)
        
            
    def add_output(self, file, opt=None, argument=False): 
        """
        Add a file, or partitioned file, as an output to this node. 
        
        Parameters
        ----------
        file : AhopeFile
            A file object that this node generates
        opt : string, optional
            The command line options that the executable needs
            in order to set the names of this file.
        argument : Boolean, optional
            If present this indicates that the file should be supplied as an
            argument to the job (ie. file name with no option at the end of
            the command line call). Using argument=True and opt != None will
            result in a failure.
        """  
        if argument and opt != None:
            errMsg = "You cannot supply an option and tell the code that this "
            errMsg += "is an argument. Choose one or the other."
            raise ValueError(errMsg)

        # make sure the files know which node creates them   
        self.output_files.append(file)
        self.add_output_file(file.path)
        file.node = self
        if opt:
            file.opt = opt
            self.add_var_opt(opt, file.path)
            
        if argument:
            if (len(file.paths) == 1):
                self.add_var_arg(file.path)
            else:
                errMsg = "Do not yet have support for taking partitioned "
                errMsg += "output files as arguments. Ask for this feature "
                errMsg += "to be added."
                
    def make_and_add_output(self, valid_seg, extension, option_name, 
                                 tags=[]):
        """
        This function will create a AhopeFile corresponding to the given
        information and then add that file as output of this node.

        Parameters
        -----------
        valid_seg : glue.segments.segment
            The time span over which the job is valid for.
        extension : string
            The extension to be used at the end of the filename. 
            E.g. '.xml' or '.sqlite'.
        option_name : string
            The option that is used when setting this job as output. For e.g.
            'output-name' or 'output-file', whatever is appropriate for the
            current executable.
        tags : list of strings, (optional, default=[])
            These tags will be added to the list of tags already associated with
            the job. They can be used to uniquely identify this output file.
        """
        job = self.job()

        # Changing this from set(tags) to enforce order. It might make sense
        # for all jobs to have file names with tags in the same order.
        all_tags = job.tags
        if extra_tags:
            for tag in tags:
                if tag not in all_tags:
                    all_tags.append(tag)

        insp = AhopeFile(job.ifo, job.exe_name, extension=extension,
                         segment=valid_seg,
                         directory=job.out_dir,
                         tags=all_tags)    

        self.add_output(insp, opt=option_name)

                
    def is_unreliable(script):
        """
        Make this job run two instances of itself and check the results using
        the given script. This is primarily used for GPU jobs where GPUs can
        produce unreliable results some of the time.
        FIXME: This has not been fully implemented yet. This mode breaks script
        output. Does pegasus know how to deal with this? Two jobs writing the
        same output file sounds dangerous.
        """
        raise NotImplementedError
        
        #FIXME this mode breaks script output
        # FIXME: Does pegasus know how to deal with this? Two jobs writing the
        # same output file sounds dangerous to me.
        
        # Make two instances 
        self.job().__queue = 2
        
        # change the output files so that they are unique
        self.unreliable = True
        
        # Call a script on the output of both 
        self.set_post_script(script)
        
    def finalize(self):
        """
        This is a stub, it does nothing at the moment, do not call it.
        FIXME: This is part of the GPU support, needs implementing.
        """
        # If the job needs to be run twice and checked, change the output 
        # name format accordingly (the postscript will produce a file of the
        # original name
        if hasattr(self, 'unreliable'):
            pass

class Executable(object):
    """
    This class is a reprentation of an executable and its capabilities.
    It can be used to create condor job(s) for this executable and is normally
    sub-classed, hence the lack of stuff in the stock instance.
    """
    def __init__(self, exe_name, universe=None):
        """
        Initialize the Executable class.

        Parameters
        -----------
        exe_name: string
            The string corresponding to the executable. We demand that this tag
            points to the executable path in the [executables] section of the
            config file and that any options in the [exe_name] section will be
            sent to the nodes resulting in the workflow.
        universe=None: string, (optional, default=None)
            The condor universe to run the job in. If not given the condor
            universe will be automatically checked by determining if condor
            libraries are present in the executable. We recommend that the
            auto-detect feature is used.
        """

        self.exe_name = exe_name
        self.condor_universe = universe

    def create_job(self, cp, ifo=None, out_dir=None, tags=[]):
        """
        Create an pycbc.ahope.Job instance for this Executable.

        Parameters
        -----------
        cp : ConfigParser.ConfigParser instance
            An in-memory representation of the configuration file
        ifo : string, (optional, default=None)
            The ifo string appropriate for the job. This is used in naming the
            output files for some executables and is necessary in some cases.
            This is normally reflected in the sub-classes of this class.
        out_dir : string, (optional, default=None)
            The path to the output directory for this job. This is used in a lot
            of cases to determine the location of the output files. Again this
            is normally reflected in the sub-classes of this class.
        Returns
        --------
        pycbc.ahope.Job instance
            The pycbc.ahope.Job instance requested.
        """
        return Job(cp, self.exe_name, self.condor_universe, ifo=ifo, 
                   out_dir=out_dir, tags=tags)

class Workflow(object):
    """
    This class manages an aHOPE style workflow. It provides convenience 
    functions for finding input files using time and keywords. It can also
    generate cache files from the inputs. It makes heavy use of the
    pipeline.CondorDAG class, which is instantiated under self.dag.
    """
    def __init__(self, args):
        """
        Create an aHOPE workflow
        
        Parameters
        ----------
        args : argparse.ArgumentParser
            The command line options to initialize an ahope workflow.
        """
        # Parse ini file
        self.cp = AhopeConfigParser.from_args(args)

        # Set global values
        start_time = int(self.cp.get("ahope","start-time"))
        end_time = int(self.cp.get("ahope", "end-time"))
        self.analysis_time = segments.segment([start_time,end_time])

        # Set the ifos to analyse
        ifos = []
        for ifo in self.cp.options('ahope-ifos'):
            ifos.append(ifo.upper())
        self.ifos = ifos
        self.ifos.sort(key=str.lower)
        self.ifoString = ''.join(self.ifos)
        
        # FIXME: Not sure if this is this is really what we should be doing to
        # get the basename!
        self.basename = basename(splitext(args.config_files[0])[0])
        
        # Initialize the dag
        logfile = self.basename + '.log'
        fh = open( logfile, "w" )
        fh.close()
        self.dag = pipeline.CondorDAG(logfile, dax=False)
        self.dag.set_dax_file(self.basename)
        self.dag.set_dag_file(self.basename)
        
        # Set up input and output file lists for workflow
        self.input_files = AhopeFileList([])
        self.output_files = AhopeFileList([])

        # Dump the parsed config file
        iniFile = os.path.abspath(self.basename + 'PARSED.ini')
        if not os.path.isfile(iniFile):
            fp = open(iniFile, 'w')
            self.cp.write(fp)
            fp.close()
        else:
            logging.warn("Cowardly refusing to overwrite %s." %(iniFile))
                 
    def add_node(self, node):
        """
        Use this function to add a pycbc.ahope.Node instance to the workflow.

        Parameters
        -----------
        node : pycbc.ahope.Node instance
            The pyCBC Node instance to be added.
        """
        self.dag.add_node(node)
        
        if node.input_files:
            self.input_files += [f for f in node.input_files 
                                if (not f.node or (f.node and f.node.executed)) 
                                and (f not in self.input_files)]
        
        if node.output_files:
            self.output_files += [f for f in node.output_files 
                                    if f not in self.output_files]
            
    def execute_node(self, node):
        """ Execute this node immediately.
        """
        node.executed = True  
        
        cmd_list = [node.job().get_executable()]
        cmd_tuples = node.get_cmd_tuple_list()   
        for cmd in cmd_tuples:
            cmd_list += list(cmd)
        cmd_list = [c for cmd in cmd_list for c in cmd.split(' ')]
        cmd_list = filter(None, cmd_list)
        job_dir = node.job().out_dir
        
        if len(node.output_files) > 0:
            base_name = node.output_files[0].filename
        else:
            base_name = node.job().exe_name
        
        # Must execute in output directory.
        currDir = os.getcwd()
        os.chdir(job_dir)
        # Make call
        make_external_call(cmd_list, outDir=os.path.join(job_dir, 'logs'),
                                     outBaseName=base_name) 
        # Change back
        os.chdir(currDir)
        
    def write_plans(self):
        """
        This will create the workflow and write it out to disk, only call this
        after the workflow has been completely created.
        """
        self.dag.write_script()
        self.dag.write_abstract_dag()
       
        # FIXME this should be in a lower level code
        f = open('output_map.dat', 'w')
        for fil in self.output_files:
            f.write(fil.map_str() + '\n')
            
        f = open('input_map.dat', 'w')
        for fil in self.input_files:
            f.write(fil.map_str() + '\n')

class AhopeFile(object):
    '''
    This class holds the details of an individual output file *or* a group
    of partitioned output files in the output workflow.
    An example of partitioned output is a template bank file split up with
    splitbank, or the matched-filter outputs from running on each of these
    split template banks in turn.
    This file(s) may be pre-supplied, generated from within the ahope
    command line script, or generated within the workflow. The important stuff
    is:

    * The ifo that the AhopeFile is valid for
    * The time span that the AhopeOutFile is valid for
    * A short description of what the file is
    * The extension that the file should have
    * The url where the file should be located

    An example of initiating this class:
    
    c = AhopeFile("H1", "INSPIRAL_S6LOWMASS", segments.segment(815901601, 815902001),file_url="file://localhost/home/spxiwh/H1-INSPIRAL_S6LOWMASS-815901601-400.xml.gz" )

    another where the file url is generated from the inputs:

    c = AhopeFile("H1", "INSPIRAL_S6LOWMASS", segments.segment(815901601, 815902001), directory="/home/spxiwh", extension="xml.gz" )
    '''
    def __init__(self, ifos, description, segs, file_url=None, 
                 extension=None, directory=None, tags=None, **kwargs):       
        """
        Create an AhopeFile instance
        
        Parameters
        ----------
        ifos : string or list
            The ifo(s) that the AhopeFile is valid for. If the file is
            independently valid for multiple ifos it can be provided as a list.
            Ie. ['H1',L1','V1'], if the file is only valid for the combination
            of ifos (for e.g. ligolw_thinca output) then this can be supplied
            as, for e.g. "H1L1V1".
        description: string
            A short description of what the file is, normally used in naming of
            the output files.
            FIXME: I think that this is now executable description, tagging
            only the program that ran this job. Can we change the name
            accordingly?
        segs : glue.segment or glue.segmentlist
            The time span that the AhopeOutFile is valid for. Note that this is
            *not* the same as the data that the job that made the file reads in.
            Lalapps_inspiral jobs do not analyse the first an last 72s of the
            data that is read, and are therefore not valid at those times. If
            the time is not continuous a segmentlist can be supplied.
        file_url : url (optional, default=None)
            If this is *not* supplied, extension and directory must be given.
            If specified this explicitly points to the url of the file, or the
            url where the file will be generated when made in the workflow.
        extension : string (optional, default=None)
            Either supply this *and* directory *or* supply only file_url.
            If given this gives the extension at the end of the file name. The
            full file name will be inferred from the other arguments
            following the ahope standard.
        directory : string (optional, default=None)
            Either supply this *and* extension *or* supply only file_url.
            If given this gives the directory in which the file exists, or will
            exists. The file name will be inferred from the other arguments
            following the ahope standard.
        tags : list of strings (optional, default=None)
            This is a list of descriptors describing what this file is. For
            e.g. this might be ["BNSINJECTIONS" ,"LOWMASS","CAT_2_VETO"].
            These are used in file naming.
        """
        # Set the science metadata on the file
        self.node=None
        if isinstance(ifos, (str, unicode)):
            self.ifoList = [ifos]
        else:
            self.ifoList = ifos
        self.ifoString = ''.join(self.ifoList)
        self.description = description
        if isinstance(segs, (segments.segment)):
            self.segList = segments.segmentlist([segs])
        elif isinstance(segs, (segments.segmentlist)):
            self.segList = segs
        else:
            errMsg = "segs input must be either glue.segments.segment or "
            errMsg += "segments.segmentlist. Got %s." %(str(type(segs)),)
            raise ValueError(errMsg)
        self.tags = tags 
        if tags is not None:
            self.tag_str = '_'.join(tags)
            tagged_description = '_'.join([description] + tags)
        else:
            tagged_description = description
        self.kwargs = kwargs
            
        # Follow the capitals-for-naming convention
        self.ifoString = self.ifoString.upper()
        self.tagged_description = tagged_description.upper()
      
        if not file_url:
            if not extension:
                raise TypeError("a file extension required if a file_url "
                                "is not provided")
            if not directory:
                raise TypeError("a directory is required if a file_url is "
                                "not provided")
            
            filename = self._filename(self.ifoString, self.tagged_description,
                                      extension, self.segList.extent())
            path = os.path.join(directory, filename)
            if not os.path.isabs(path):
                path = os.path.join(os.getcwd(), path) 
            file_url = urlparse.urlunparse(['file', 'localhost', path, None,
                                            None, None])
       
        if not isinstance(file_url, list):
            file_url = [file_url]
       
        self.cache_entries = []
        for url in file_url:
            cache_entry = lal.CacheEntry(self.ifoString,
                       self.tagged_description, self.segList.extent(), url)
            self.cache_entries.append(cache_entry)   
    
    @property
    def url(self):
        return self.cache_entries[0].url
       
    @property
    def path(self):
        """
        If only one file is contained in this instance this will be that path.
        Otherwise a TypeError is raised.
        """
        return self.cache_entries[0].path

    @property
    def ifo(self):
        """
        If only one ifo in the ifoList this will be that ifo. Otherwise an
        error is raised.
        """
        if len(self.ifoList) == 1:
            return self.ifoList[0]
        else:
            errMsg = "self.ifoList must contain only one ifo to access the "
            errMsg += "ifo property. %s." %(str(self.ifoList),)
            raise TypeError(errMsg)

    @property
    def segment(self):
        """
        If only one segment in the segmentlist this will be that segment.
        Otherwise an error is raised.
        """
        if len(self.segList) == 1:
            return self.segList[0]
        else:
            errMsg = "self.segList must only contain one segment to access"
            errMsg += " the segment property. %s." %(str(self.segList),)
            raise TypeError(errMsg)
        
    @property
    def filename(self):
        """
        If only one file is contained in this instance this will be that
        file's name. Otherwise a TypeError is raised.
        """
        return basename(self.cache_entries[0].path)
        
    def _filename(self, ifo, description, extension, segment):
        """
        Construct the standard output filename. Should only be used internally
        of the AhopeFile class.
        """        
        if extension.startswith('.'):
            extension = extension[1:]
        # Follow the frame convention of using integer filenames, but stretching
        # to cover partially covered seconds.
        start = int(segment[0])
        end = int(math.ceil(segment[1]))
        duration = str(end-start)
        start = str(start)
        
        return "%s-%s-%s-%s.%s" % (ifo, description.upper(), start, duration, extension)    
       
    @property
    def lfn(self):
        return self.filename
        
    @property
    def pfn(self):
        return self.path

    def map_str(self, site='local'):
        return '%s %s pool="%s"' % (self.lfn, self.pfn, site) 
    
class AhopeFileList(list):
    '''
    This class holds a list of AhopeFile objects. It inherits from the
    built-in list class, but also allows a number of features. ONLY
    AhopeFile instances should be within an AhopeFileList instance.
    '''
    entry_class = AhopeFile

    def find_output(self, ifo, time):
        '''
        Return one AhopeFile that covers the given time, or is most
        appropriate for the supplied time range.

        Parameters
        -----------
        ifo : string
           Name of the ifo (or ifos) that the file should be valid for.
        time : int/float/LIGOGPStime or tuple containing two values
           If int/float/LIGOGPStime (or similar may of specifying one time) is
           given, return the AhopeFile corresponding to the time. This calls
           self.find_output_at_time(ifo,time).
           If a tuple of two values is given, return the AhopeFile that is
           **most appropriate** for the time range given. This calls
           self.find_output_in_range

        Returns
        --------
        AhopeFile class
           The AhopeFile that corresponds to the time/time range
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

    def find_output_at_time(self, ifo, time):
       '''
       Return AhopeFile that covers the given time.

        Parameters
        -----------
        ifo : string
           Name of the ifo (or ifos) that the AhopeFile should correspond to
        time : int/float/LIGOGPStime
           Return the AhopeFiles that covers the supplied time. If no
           AhopeFile covers the time this will return None.

        Returns
        --------
        list of AhopeFile classes
           The AhopeFiles that corresponds to the time.
        '''
       # Get list of AhopeFiles that overlap time, for given ifo
       outFiles = [i for i in self if\
                                    ifo in i.ifoList and time in i.segList] 
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

    def find_outputs_in_range(self, ifo, current_segment):
        """
        Return the list of AhopeFiles that is most appropriate for the supplied
        time range. That is, the AhopeFiles whose coverage time has the
        largest overlap with the supplied time range.

        Parameters
        -----------
        ifo : string
           Name of the ifo (or ifos) that the AhopeFile should correspond to
        current_segment : glue.segment.segment
           The segment of time that files must intersect.

        Returns
        --------
        AhopeFileList class
           The list of AhopeFiles that are most appropriate for the time range
        """
        currSegList = segments.segmentlist([current_segment])

        # First filter AhopeFiles corresponding to ifo
        ifo_files = [i for i in self if ifo in i.ifoList] 

        # Filter AhopeOutFiles to those overlapping the given window
        overlap_files = [i for i in ifo_files \
                          if i.segList.intersects_segment(current_segment)]
        overlap_windows = [abs(i.segList & currSegList) \
                                                        for i in overlap_files]

        # FIXME: Error handling for the overlap_files == [] case?

        # Return the AhopeFile with the biggest overlap
        # Note if two AhopeFile have identical overlap, the first is used
        # to define the valid segment
        overlap_windows = numpy.array(overlap_windows, dtype = int)
        segmentLst = overlap_files[overlap_windows.argmax()].segList
        
        # Get all output files with the exact same segment definition
        output_files = [f for f in overlap_files if f.segList==segmentLst]
        return output_files

    def find_output_in_range(self, ifo, start, end):
        '''
        Return the AhopeFile that is most appropriate for the supplied
        time range. That is, the AhopeFile whose coverage time has the
        largest overlap with the supplied time range. If no AhopeFiles
        overlap the supplied time window, will return None. 

        Parameters
        -----------
        ifo : string
           Name of the ifo (or ifos) that the AhopeFile should correspond to
        start : int/float/LIGOGPStime 
           The start of the time range of interest.
        end : int/float/LIGOGPStime
           The end of the time range of interest

        Returns
        --------
        AhopeFile class
           The AhopeFile that is most appropriate for the time range
        '''
        currSegList = segments.segmentlist([current_segment])

        # First filter AhopeFiles corresponding to ifo
        outFiles = [i for i in self if ifo in i.ifoList]

        if len(outFiles) == 0:
            # No AhopeOutFiles correspond to that ifo
            return None
        # Filter AhopeOutFiles to those overlapping the given window
        currSeg = segments.segment([start,end])
        outFiles = [i for i in outFiles \
                                  if i.segList.intersects_segment(currSeg)]

        if len(outFiles) == 0:
            # No AhopeOutFile overlap that time period
            return None
        elif len(outFiles) == 1:
            # One AhopeOutFile overlaps that period
            return outFiles[0]
        else:
            overlap_windows = [abs(i.segList & currSegList) \
                                                        for i in outFiles]
            # Return the AhopeFile with the biggest overlap
            # Note if two AhopeFile have identical overlap, this will return
            # the first AhopeFile in the list
            overlap_windows = numpy.array(overlap_windows, dtype = int)
            return outFiles[overlap_windows.argmax()]

    def find_all_output_in_range(self, ifo, currSeg):
        """
        Return all files that overlap the specified segment.
        """
        outFiles = [i for i in self if ifo in i.ifoList]
        outFiles = [i for i in outFiles \
                                      if i.segList.intersects_segment(currSeg)]
        return self.__class__(outFiles)

    def find_output_with_tag(self, tag):
        """
        Find all files who have tag in self.tags
        """
        return AhopeFileList([i for i in self if tag in i.tags])

    def find_output_with_ifo(self, ifo):
        """
        Find all files who have ifo = ifo
        """
        return AhopeFileList([i for i in self if ifo in i.ifoList])

    def get_times_covered_by_files(self):
        """
        Find the coalesced intersection of the segments of all files in the
        list.
        """
        times = segments.segmentlist([])
        for entry in self:
            times.extend(entry.segList)
        times.coalesce()
        return times

    def convert_to_lal_cache(self):
        """
        Return all files in this object as a lal.Cache object
        """
        lalCache = lal.Cache([])
        for entry in self:
            lalCache.extend(entry.cache_entries)
        return lalCache


class AhopeOutSegFile(AhopeFile):
    '''
    This class inherits from the AhopeFile class, and is designed to store
    ahope output files containing a segment list. This is identical in
    usage to AhopeFile except for an additional kwarg for holding the
    segment list, if it is known at ahope run time.
    '''
    def __init__(self, ifo, description, segment, fileUrl,
                 segList=None, **kwargs):
        """
        See AhopeFile.__init__ for a full set of documentation for how to
        call this class. The only thing unique and added to this class is
        the required option timeSeg, as described below:

        Parameters:
        ------------
        ifo : string or list (required)
            See AhopeFile.__init__
        description : string (required)
            See AhopeFile.__init__
        segment : glue.segments.segment or glue.segments.segmentlist
            See AhopeFile.__init__
        fileUrl : string (required)
            See AhopeFile.__init__
            FIXME: This is a kwarg in AhopeFile and should be here as well,
            if this is removed from the explicit arguments it would allow for
            either fileUrls or constructed file names to be used in AhopeFile.
        segList : glue.segments.segmentlist (optional, default=None)
            A glue.segments.segmentlist covering the times covered by the
            segmentlist associated with this file. If this is the science time
            or CAT_1 file this will be used to determine analysis time. Can
            be added by setting self.segList after initializing an instance of
            the class.

        """
        AhopeFile.__init__(self, ifo, description, segment, fileUrl,
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
        self.segmentList = newSegList
        self.toSegmentXml()

    def toSegmentXml(self):
        """
        Write the segment list in self.segmentList to the url in self.url.
        """
        filePointer = open(self.path, 'w')
        dqUtils.tosegmentxml(filePointer, self.segmentList)
        filePointer.close()

def make_external_call(cmdList, outDir=None, outBaseName='external_call',
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
        raise CalledProcessErrorMod(errCode, ' '.join(cmdList), 
                errFile=errFile, outFile=outFile, cmdFile=cmdFile)
    logging.debug("Call successful, or error checking disabled.")

class CalledProcessErrorMod(Exception):
    """
    This exception is raised when subprocess.call returns a non-zero exit code
    and checking has been requested. This should not be accessed by the user
    it is used only within make_external_call.
    """
    def __init__(self, returncode, cmd, errFile=None, outFile=None, 
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
            msg += "Stderr can be found in %s .\n" %(self.errFile)
        if self.outFile:
            msg += "Stdout can be found in %s .\n" %(self.outFile)
        if self.cmdFile:
            msg += "The failed command has been printed in %s ." %(self.cmdFile)
        return msg
              
def get_full_analysis_chunk(science_segs):
    """
    Function to find the first and last time point contained in the science segments
    and return a single segment spanning that full time.

    Parameters
    -----------
    science_segs : ifo-keyed dictionary of glue.segments.segmentlist instances
        The list of times that are being analysed in this workflow.
    Returns
    --------
    fullSegment : glue.segments.segment
        The segment spanning the first and last time point contained in science_segs.
    """
    extents = [science_segs[ifo].extent() for ifo in science_segs.keys()]
    min, max = extents[0]
    for lo, hi in extents:
        if min > lo:
            min = lo
        if max < hi:
            max = hi
    fullSegment = segments.segment(min, max)
    return fullSegment
        
