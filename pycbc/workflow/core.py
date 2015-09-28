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
creating a workflow. For details about the workflow module see here:
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/ahope.html
"""
import os, stat, subprocess, logging, math, string, urllib2, urlparse, ConfigParser, copy, time
import numpy, cPickle, random
from itertools import combinations, groupby
from operator import attrgetter
import lal as lalswig
from glue import lal, segments
from pycbc.workflow.configuration import WorkflowConfigParser
from pycbc.workflow import pegasus_workflow

# workflow should never be using the glue LIGOTimeGPS class, override this with
# the nice SWIG-wrapped class in lal
lal.LIGOTimeGPS = lalswig.LIGOTimeGPS

#REMOVE THESE FUNCTIONS  FOR PYTHON >= 2.7 ####################################
def check_output_error_and_retcode(*popenargs, **kwargs):
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
    output, error = process.communicate()
    retcode = process.poll()
    return output, error, retcode

def check_output(*popenargs, **kwargs):
    output, unused_error, unused_retcode = \
                           check_output_error_and_retcode(*popenargs, **kwargs)
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
        
class Executable(pegasus_workflow.Executable):
    # These are the file retention levels
    INTERMEDIATE_PRODUCT = 1
    NON_CRITICAL = 2
    CRITICAL = 3
    FINAL_RESULT = 4
    
    # This is the default value. It will give a warning if a class is
    # used where the retention level is not set. The file will still be stored
    KEEP_BUT_RAISE_WARNING = 5
    _warned_classes_list = ['Executable']

    # Sub classes, or instances, should override this. If not overriden the
    # file will be retained, but a warning given
    current_retention_level = KEEP_BUT_RAISE_WARNING
    def __init__(self, cp, name, 
                       universe=None, ifos=None, out_dir=None, tags=[]):
        """
        Initialize the Executable class.

        Parameters
        -----------
        cp : ConfigParser object
            The ConfigParser object holding the workflow configuration settings
        exec_name : string
            Executable name
        universe : string, optional
            Condor universe to run the job in
        ifos : string or list, optional
            The ifo(s) that the Job is valid for. If the job is
            independently valid for multiple ifos it can be provided as a list.
            Ie. ['H1',L1','V1'], if the job is only valid for the combination
            of ifos (for e.g. ligolw_thinca) then this can be supplied
            as, for e.g. "H1L1V1".
        out_dir: path, optional
            The folder to store output files of this job. 
        tags : list of strings
            A list of strings that is used to identify this job.
        """
        tags = [tag.upper() for tag in tags]
        self.tags = tags
        if isinstance(ifos, (str, unicode)):
            self.ifo_list = [ifos]
        else:
            self.ifo_list = ifos
        if self.ifo_list is not None:
            self.ifo_string = ''.join(self.ifo_list)
        else:
            self.ifo_string = None
        self.cp = cp
        self.universe=universe
        
        if len(tags) != 0:
            self.tagged_name = "%s-%s" % (name, '_'.join(tags))
        else:
            self.tagged_name = name
        if self.ifo_string is not None:
            self.tagged_name = "%s-%s" % (self.tagged_name, self.ifo_string)

        try:
            self.installed = cp.getboolean('pegasus_profile-%s' % name, 'pycbc|installed')
        except:
            self.installed = True

        super(Executable, self).__init__(self.tagged_name, installed=self.installed)
        
        self.name=name
        
        # Determine the output directory
        if out_dir is not None:
            self.out_dir = out_dir
        elif len(tags) == 0:
            self.out_dir = name
        else:
            self.out_dir = self.tagged_name            
        if not os.path.isabs(self.out_dir):
            self.out_dir = os.path.join(os.getcwd(), self.out_dir) 
              
        # Check that the executable actually exists locally or 
        # looks like a URL, in which case trust Pegasus to be
        # able to fetch it.
        exe_path = cp.get('executables', name)
        valid_path = False

        if exe_path.find('://') > 0:
            if exe_path.startswith('file://'):
                valid_path = os.path.isfile(exe_path[7:])
            else:
                valid_path = True
        else:
            valid_path = os.path.isfile(exe_path)

        if valid_path:
            logging.debug("Using %s executable "
                          "at %s" % (name, exe_path))
        else:
            raise TypeError("Failed to find %s executable " 
                            "at %s" % (name, exe_path))
        
        self.add_pfn(exe_path)

        # Determine the condor universe if we aren't given one 
        if self.universe is None:
            if is_condor_exec(exe_path):
                self.universe = 'standard'
            else:
                self.universe = 'vanilla'
                
        logging.debug("%s executable will run as %s universe"
                     % (name, self.universe))  
    
        self.set_universe(self.universe)

        # Determine the sections from the ini file that will configure
        # this executable
        sections = [name]
        if self.ifo_string:
            sec_tags = tags + [self.ifo_string]
        else:
            sec_tags = tags
        for tag in sec_tags:
             section = '%s-%s' %(name, tag.lower())
             if cp.has_section(section):
                sections.append(section)
        self.sections = sections   
        # Do some basic sanity checking on the options      
        for sec1, sec2 in combinations(sections, 2):
            cp.check_duplicate_options(sec1, sec2, raise_error=True)
             
        # collect the options and profile information 
        # from the ini file section(s)
        self.common_options = []        
        for sec in sections:
            if cp.has_section(sec):
                self.add_ini_opts(cp, sec)
            else:
                warnString = "warning: config file is missing section [%s]"\
                             %(sec,)
                logging.warn(warnString)
            if cp.has_section('pegasus_profile-%s' % sec):
                self.add_ini_profile(cp, 'pegasus_profile-%s' % sec)
                
        # Add executable non-specific profile information
        if cp.has_section('pegasus_profile'):
            self.add_ini_profile(cp, 'pegasus_profile')

        # Determine the level at which output files should be kept
        try:
            global_retention_level = \
                cp.get_opt_tags("workflow", "file-retention-level",
                                   tags+[name])
        except:
            msg="Cannot find file-retention-level in [workflow] section "
            msg+="of the configuration file. Setting a default value of "
            msg+="retain all files."
            logging.warn(msg)
            self.retain_files = True
            self.global_retention_threshold = 1
            cp.set("workflow", "file-retention-level", "all_files")
        else:
            # FIXME: Are these names suitably descriptive?
            if global_retention_level == 'all_files':
                self.global_retention_threshold = 1
            elif global_retention_level == 'no_intermediates':
                self.global_retention_threshold = 2
            elif global_retention_level == 'all_triggers':
                self.global_retention_threshold = 3
            elif global_retention_level == 'results_only':
                self.global_retention_threshold = 4
            else:
                err_msg = "Cannot recognize the file-retention-level in the "
                err_msg += "[workflow] section of the ini file. "
                err_msg += "Got : %s." %(global_retention_level,)
                err_msg += "Valid options are: 'all_files', 'no_intermediates',"
                err_msg += "'all_triggers' or 'results_only' "
                raise ValueError(err_msg)
            if self.current_retention_level == 5:
                self.retain_files = True
                if type(self).__name__ in Executable._warned_classes_list:
                    pass
                else:
                    warn_msg = "Attribute current_retention_level has not "
                    warn_msg += "been set in class %s. " %(type(self),)
                    warn_msg += "This value should be set explicitly. "
                    warn_msg += "All output from this class will be stored."
                    logging.warn(warn_msg)
                    Executable._warned_classes_list.append(type(self).__name__)
            elif self.global_retention_threshold > self.current_retention_level:
                self.retain_files = False
            else:
                self.retain_files = True
                
        if hasattr(self, "group_jobs"):
            self.add_profile('pegasus', 'clusters.size', self.group_jobs)        

    @property
    def ifo(self):
        """
        If only one ifo in the ifo list this will be that ifo. Otherwise an
        error is raised.
        """
        if self.ifo_list and len(self.ifo_list) == 1:
            return self.ifo_list[0]
        else:
            errMsg = "self.ifoList must contain only one ifo to access the "
            errMsg += "ifo property. %s." %(str(self.ifo_list),)
            raise TypeError(errMsg)

    def add_ini_profile(self, cp, sec):
        for opt in cp.options(sec):
            namespace = opt.split('|')[0]
            if namespace == 'pycbc':
                continue

            value = string.strip(cp.get(sec, opt))
            key = opt.split('|')[1]
            self.add_profile(namespace, key, value)

            # Remove if Pegasus can apply this hint in the TC
            if namespace == 'hints' and key == 'execution.site':
                self.execution_site = value

    def add_ini_opts(self, cp, sec):
        for opt in cp.options(sec):
            value = string.strip(cp.get(sec, opt))
            self.common_options += ['--%s' % opt, value]
            
    def add_opt(self, opt, value=None):
        if value is None:
            self.common_options += [opt]
        else:
            self.common_options += [opt, value]
            
    def get_opt(self, opt):
        for sec in self.sections:
            try:
                key = self.cp.get(sec, opt)
                if key:
                    return key
            except ConfigParser.NoOptionError:
                pass

        return None

    def has_opt(self, opt):
        for sec in self.sections:
            val = self.cp.has_option(sec, opt)
            if val:
                return val

        return False

    def create_node(self):
        """ Default node constructor. This is usually overridden by subclasses
        of Executable.
        """
        return Node(self)

class Workflow(pegasus_workflow.Workflow):
    """
    This class manages a pycbc workflow. It provides convenience 
    functions for finding input files using time and keywords. It can also
    generate cache files from the inputs.
    """
    def __init__(self, args, name):
        """
        Create a pycbc workflow
        
        Parameters
        ----------
        args : argparse.ArgumentParser
            The command line options to initialize a CBC workflow.
        """
        super(Workflow, self).__init__(name)
        
        # Parse ini file
        self.cp = WorkflowConfigParser.from_args(args)
        
        # Dump the parsed config file
        symlink = os.path.abspath(self.name + '_parsed.ini')

        if os.path.isfile(symlink):
            os.remove(symlink)

        ini_file = os.path.abspath(self.name + '_parsed_%d.ini' % time.time())
        # This shouldn't already exist, but just in case
        if os.path.isfile(ini_file):
            os.remove(ini_file)

        fp = open(ini_file, 'w')
        self.cp.write(fp)
        fp.close()

        os.symlink(ini_file, symlink)

        # Set global values
        start_time = int(self.cp.get("workflow", "start-time"))
        end_time = int(self.cp.get("workflow", "end-time"))
        self.analysis_time = segments.segment([start_time, end_time])

        # Set the ifos to analyse
        ifos = []
        for ifo in self.cp.options('workflow-ifos'):
            ifos.append(ifo.upper())
        self.ifos = ifos
        self.ifos.sort(key=str.lower)
        self.ifo_string = ''.join(self.ifos)
        
        # Set up input and output file lists for workflow
        self._inputs = FileList([])
        self._outputs = FileList([])
 
    @property
    def output_map(self):  
        if self.in_workflow != False:
            name = self.name + '.map'
        else:
            name = 'output.map'
        path =  os.path.join(os.getcwd(), name)
        return path
        

    def execute_node(self, node, verbatim_exe = False):
        """ Execute this node immediately on the local machine
        """
        node.executed = True
        cmd_list = node.get_command_line(verbatim_exe=verbatim_exe)
        
        # Must execute in output directory.
        curr_dir = os.getcwd()
        out_dir = node.executable.out_dir
        os.chdir(out_dir)
        
        # Make call
        make_external_call(cmd_list, out_dir=os.path.join(out_dir, 'logs'),
                                     out_basename=node.executable.name) 
        # Change back
        os.chdir(curr_dir)

        for fil in node._outputs:
            fil.node = None
            fil.PFN(fil.storage_path, site='local')
            
    def save(self):
        self.as_job.addArguments('-Dpegasus.dir.storage.mapper.replica.file=%s' % self.output_map) 
        self.as_job.addArguments('-Dpegasus.dir.storage.mapper.replica=File') 
        self.as_job.addArguments('--cache %s' % os.path.join(os.getcwd(), '_reuse.cache')) 
        self.as_job.addArguments('--output-site local')     
        self.as_job.addArguments('--cleanup inplace')
        self.as_job.addArguments('--cluster label,horizontal')

        # add executable pfns for local site to dax
        for exe in self._executables:
            exe.insert_into_dax(self._adag)
            
        # add workflow input files pfns for local site to dax
        for fil in self._inputs:
            fil.insert_into_dax(self._adag)
            
        # save the dax file
        super(Workflow, self).save()
        
        # add workflow storage locations to the output mapper
        f = open(self.output_map, 'w')
        for out in self._outputs:
            try:
                f.write(out.output_map_str() + '\n')
            except ValueError:
                # There was no storage path
                pass
    
class Node(pegasus_workflow.Node):
    def __init__(self, executable):
        super(Node, self).__init__(executable)
        self.executed = False
        self.set_category(executable.name)
        
        if executable.universe == 'vanilla' and executable.installed:
            self.add_profile('condor', 'getenv', 'True')
        
        if hasattr(executable, 'execution_site'):
            self.add_profile('hints', 'execution.site', executable.execution_site)
            self.add_profile('hints', 'executionPool', executable.execution_site)
            
        self._options += self.executable.common_options
    
    def get_command_line(self, verbatim_exe=False):
        self._finalize()
        arglist = self._dax_node.arguments
        
        tmpargs = []
        for a in arglist:
            if not isinstance(a, File):
                tmpargs += a.split(' ')
            else:
                tmpargs.append(a)
        arglist = tmpargs
        
        arglist = [a for a in arglist if a != '']
        
        arglist = [a.storage_path if isinstance(a, File) else a for a in arglist]
       
        # This allows the pfn to be an http(s) URL, which will be
        # downloaded by resolve_url
        if verbatim_exe:
            exe_path = self.executable.get_pfn()
        else:
            exe_path = urlparse.urlsplit(self.executable.get_pfn()).path

        return [exe_path] + arglist
        
    def new_output_file_opt(self, valid_seg, extension, option_name, tags=[],
                            store_file=None, use_tmp_subdirs=False):
        """
        This function will create a workflow.File object corresponding to the given
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
        store_file : Boolean, (optional, default=True)
            This file is to be added to the output mapper and will be stored
            in the specified output location if True. If false file will be
            removed when no longer needed in the workflow.
        """
        
        # Changing this from set(tags) to enforce order. It might make sense
        # for all jobs to have file names with tags in the same order.
        all_tags = copy.deepcopy(self.executable.tags)
        for tag in tags:
            if tag not in all_tags:
                all_tags.append(tag)

        store_file = store_file if store_file is not None else self.executable.retain_files

        fil = File(self.executable.ifo_list, self.executable.name,
                   valid_seg, extension=extension, store_file=store_file, 
                   directory=self.executable.out_dir, tags=all_tags,
                   use_tmp_subdirs=use_tmp_subdirs)
        self.add_output_opt(option_name, fil)
        return fil
        
    @property    
    def output_files(self):
        return FileList(self._outputs)

    @property
    def output_file(self):
        """
        If only one output file return it. Otherwise raise an exception.
        """
        out_files = self.output_files
        if len(out_files) != 1:
            err_msg = "output_file property is only valid if there is a single"
            err_msg += " output file. Here there are "
            err_msg += "%d output files." %(len(out_files))
            raise ValueError(err_msg)
        return out_files[0]
    
class File(pegasus_workflow.File):
    '''
    This class holds the details of an individual output file 
    This file(s) may be pre-supplied, generated from within the workflow
    command line script, or generated within the workflow. The important stuff
    is:

    * The ifo that the File is valid for
    * The time span that the OutFile is valid for
    * A short description of what the file is
    * The extension that the file should have
    * The url where the file should be located

    An example of initiating this class:
    
    c = File("H1", "INSPIRAL_S6LOWMASS", segments.segment(815901601, 815902001), file_url="file://localhost/home/spxiwh/H1-INSPIRAL_S6LOWMASS-815901601-400.xml.gz" )

    another where the file url is generated from the inputs:

    c = File("H1", "INSPIRAL_S6LOWMASS", segments.segment(815901601, 815902001), directory="/home/spxiwh", extension="xml.gz" )
    '''
    def __init__(self, ifos, exe_name, segs, file_url=None, 
                 extension=None, directory=None, tags=None, 
                 store_file=True, use_tmp_subdirs=False):
        """
        Create a File instance
        
        Parameters
        ----------
        ifos : string or list
            The ifo(s) that the File is valid for. If the file is
            independently valid for multiple ifos it can be provided as a list.
            Ie. ['H1',L1','V1'], if the file is only valid for the combination
            of ifos (for e.g. ligolw_thinca output) then this can be supplied
            as, for e.g. "H1L1V1".
        exe_name: string
            A short description of the executable description, tagging
            only the program that ran this job.
        segs : glue.segment or glue.segmentlist
            The time span that the OutFile is valid for. Note that this is
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
            following the workflow standard.
        directory : string (optional, default=None)
            Either supply this *and* extension *or* supply only file_url.
            If given this gives the directory in which the file exists, or will
            exists. The file name will be inferred from the other arguments
            following the workflow standard.
        tags : list of strings (optional, default=None)
            This is a list of descriptors describing what this file is. For
            e.g. this might be ["BNSINJECTIONS" ,"LOWMASS","CAT_2_VETO"].
            These are used in file naming.
        """
        self.metadata = {}
        
        # Set the science metadata on the file
        if isinstance(ifos, (str, unicode)):
            self.ifo_list = [ifos]
        else:
            self.ifo_list = ifos
        self.ifo_string = ''.join(self.ifo_list)
        self.description = exe_name
        
        if isinstance(segs, (segments.segment)):
            self.segment_list = segments.segmentlist([segs])
        elif isinstance(segs, (segments.segmentlist)):
            self.segment_list = segs
        else:
            err = "segs input must be either glue.segments.segment or "
            err += "segments.segmentlist. Got %s." %(str(type(segs)),)
            raise ValueError(err)
            
        self.tags = tags 
        if tags is not None:
            self.tag_str = '_'.join(tags)
            tagged_description = '_'.join([self.description] + tags)
        else:
            tagged_description = self.description
            
        # Follow the capitals-for-naming convention
        self.ifo_string = self.ifo_string.upper()
        self.tagged_description = tagged_description.upper()
      
        if not file_url:
            if not extension:
                raise TypeError("a file extension required if a file_url "
                                "is not provided")
            if not directory:
                raise TypeError("a directory is required if a file_url is "
                                "not provided")
            
            filename = self._filename(self.ifo_string, self.tagged_description,
                                      extension, self.segment_list.extent())
            path = os.path.join(directory, filename)
            if not os.path.isabs(path):
                path = os.path.join(os.getcwd(), path) 
            file_url = urlparse.urlunparse(['file', 'localhost', path, None,
                                            None, None])

        # Let's do a test here
        if use_tmp_subdirs and len(self.segment_list):
            pegasus_lfn = str(int(self.segment_list.extent()[0]))[:-4]
            pegasus_lfn = pegasus_lfn + '/' + os.path.basename(file_url)
        else:
            pegasus_lfn = os.path.basename(file_url)
        super(File, self).__init__(pegasus_lfn)
        
        if store_file:
            self.storage_path = urlparse.urlsplit(file_url).path
        else:
            self.storage_path = None

    def __getstate__(self):
        """ Allow the ahope file to be picklable. This disables the usage of
        the internal cache entry.
        """
        for i, seg in enumerate(self.segment_list):
            self.segment_list[i] = segments.segment(float(seg[0]), float(seg[1]))
        self.cache_entry = None
        safe_dict = copy.copy(self.__dict__)
        safe_dict['cache_entry'] = None
        return safe_dict   

    def add_metadata(self, key, value):
        """ Add arbitrary metadata to this file """
        self.metadata[key] = value

    @property
    def ifo(self):
        """
        If only one ifo in the ifo_list this will be that ifo. Otherwise an
        error is raised.
        """
        if len(self.ifo_list) == 1:
            return self.ifo_list[0]
        else:
            err = "self.ifo_list must contain only one ifo to access the "
            err += "ifo property. %s." %(str(self.ifo_list),)
            raise TypeError(err)

    @property
    def segment(self):
        """
        If only one segment in the segmentlist this will be that segment.
        Otherwise an error is raised.
        """
        if len(self.segment_list) == 1:
            return self.segment_list[0]
        else:
            err = "self.segment_list must only contain one segment to access"
            err += " the segment property. %s." %(str(self.segment_list),)
            raise TypeError(err)

    @property
    def cache_entry(self):
        """
        Returns a CacheEntry instance for File.
        """
        if self.storage_path is None:
            raise ValueError('This file is temporary and so a lal '
                             'cache entry cannot be made')
            
        file_url = urlparse.urlunparse(['file', 'localhost', self.storage_path, None,
                                            None, None])
        cache_entry = lal.CacheEntry(self.ifo_string,
                   self.tagged_description, self.segment_list.extent(), file_url)
        cache_entry.workflow_file = self
        return cache_entry

    def _filename(self, ifo, description, extension, segment):
        """
        Construct the standard output filename. Should only be used internally
        of the File class.
        """        
        if extension.startswith('.'):
            extension = extension[1:]
            
        # Follow the frame convention of using integer filenames,
        # but stretching to cover partially covered seconds.
        start = int(segment[0])
        end = int(math.ceil(segment[1]))
        duration = str(end-start)
        start = str(start)
        
        return "%s-%s-%s-%s.%s" % (ifo, description.upper(), start, 
                                   duration, extension)  
    
class FileList(list):
    '''
    This class holds a list of File objects. It inherits from the
    built-in list class, but also allows a number of features. ONLY
    pycbc.workflow.File instances should be within a FileList instance.
    '''
    entry_class = File

    def categorize_by_attr(self, attribute):
        '''
        Function to categorize a FileList by a File object
        attribute (eg. 'segment', 'ifo', 'description').

        Parameters
        -----------
        attribute : string
           File object attribute to categorize FileList

        Returns
        --------
        keys : list
           A list of values for an attribute
        groups : list
           A list of FileLists
        '''

        # need to sort FileList otherwise using groupby without sorting does
        # 'AAABBBCCDDAABB' -> ['AAA','BBB','CC','DD','AA','BB']
        # and using groupby with sorting does
        # 'AAABBBCCDDAABB' -> ['AAAAA','BBBBB','CC','DD']
        flist = sorted(self, key=attrgetter(attribute), reverse=True)

        # use groupby to create lists
        groups = []
        keys = []
        for k, g in groupby(flist, attrgetter(attribute)):
            groups.append(FileList(g))
            keys.append(k)

        return keys, groups

    def find_output(self, ifo, time):
        '''
        Return one File that covers the given time, or is most
        appropriate for the supplied time range.

        Parameters
        -----------
        ifo : string
           Name of the ifo (or ifos) that the file should be valid for.
        time : int/float/LIGOGPStime or tuple containing two values
           If int/float/LIGOGPStime (or similar may of specifying one time) is
           given, return the File corresponding to the time. This calls
           self.find_output_at_time(ifo,time).
           If a tuple of two values is given, return the File that is
           **most appropriate** for the time range given. This calls
           self.find_output_in_range

        Returns
        --------
        File class
           The File that corresponds to the time/time range
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
       Return File that covers the given time.

        Parameters
        -----------
        ifo : string
           Name of the ifo (or ifos) that the File should correspond to
        time : int/float/LIGOGPStime
           Return the Files that covers the supplied time. If no
           File covers the time this will return None.

        Returns
        --------
        list of File classes
           The Files that corresponds to the time.
        '''
       # Get list of Files that overlap time, for given ifo
       outFiles = [i for i in self if ifo in i.ifo_list and time in i.segment_list] 
       if len(outFiles) == 0:
           # No OutFile at this time
           return None
       elif len(outFiles) == 1:
           # 1 OutFile at this time (good!)
           return outFiles
       else:
           # Multiple output files. Currently this is valid, but we may want
           # to demand exclusivity later, or in certain cases. Hence the
           # separation.
           return outFiles

    def find_outputs_in_range(self, ifo, current_segment, useSplitLists=False):
        """
        Return the list of Files that is most appropriate for the supplied
        time range. That is, the Files whose coverage time has the
        largest overlap with the supplied time range.

        Parameters
        -----------
        ifo : string
           Name of the ifo (or ifos) that the File should correspond to
        current_segment : glue.segment.segment
           The segment of time that files must intersect.

        Returns
        --------
        FileList class
           The list of Files that are most appropriate for the time range
        """
        currsegment_list = segments.segmentlist([current_segment])

        # Get all files overlapping the window
        overlap_files = self.find_all_output_in_range(ifo, current_segment,
                                                    useSplitLists=useSplitLists)

        # By how much do they overlap?
        overlap_windows = [abs(i.segment_list & currsegment_list) for i in overlap_files]

        if not overlap_windows:
            return []

        # Return the File with the biggest overlap
        # Note if two File have identical overlap, the first is used
        # to define the valid segment
        overlap_windows = numpy.array(overlap_windows, dtype = int)
        segmentLst = overlap_files[overlap_windows.argmax()].segment_list
        
        # Get all output files with the exact same segment definition
        output_files = [f for f in overlap_files if f.segment_list==segmentLst]
        return output_files

    def find_output_in_range(self, ifo, start, end):
        '''
        Return the File that is most appropriate for the supplied
        time range. That is, the File whose coverage time has the
        largest overlap with the supplied time range. If no Files
        overlap the supplied time window, will return None. 

        Parameters
        -----------
        ifo : string
           Name of the ifo (or ifos) that the File should correspond to
        start : int/float/LIGOGPStime 
           The start of the time range of interest.
        end : int/float/LIGOGPStime
           The end of the time range of interest

        Returns
        --------
        File class
           The File that is most appropriate for the time range
        '''
        currsegment_list = segments.segmentlist([current_segment])

        # First filter Files corresponding to ifo
        outFiles = [i for i in self if ifo in i.ifo_list]

        if len(outFiles) == 0:
            # No OutFiles correspond to that ifo
            return None
        # Filter OutFiles to those overlapping the given window
        currSeg = segments.segment([start,end])
        outFiles = [i for i in outFiles \
                                  if i.segment_list.intersects_segment(currSeg)]

        if len(outFiles) == 0:
            # No OutFile overlap that time period
            return None
        elif len(outFiles) == 1:
            # One OutFile overlaps that period
            return outFiles[0]
        else:
            overlap_windows = [abs(i.segment_list & currsegment_list) \
                                                        for i in outFiles]
            # Return the File with the biggest overlap
            # Note if two File have identical overlap, this will return
            # the first File in the list
            overlap_windows = numpy.array(overlap_windows, dtype = int)
            return outFiles[overlap_windows.argmax()]

    def find_all_output_in_range(self, ifo, currSeg, useSplitLists=False):
        """
        Return all files that overlap the specified segment.
        """
        if not useSplitLists:
            # Slower, but simpler method
            outFiles = [i for i in self if ifo in i.ifo_list]
            outFiles = [i for i in outFiles \
                                      if i.segment_list.intersects_segment(currSeg)]
        else:
            # Faster, but more complicated
            # Basically only check if a subset of files intersects_segment by
            # using a presorted list. Sorting only happens once.
            if not self._check_split_list_validity():
                # FIXME: DO NOT hard code this.
                self._temporal_split_list(100)
            startIdx = int( (currSeg[0] - self._splitListsStart) / \
                                                          self._splitListsStep )
            # Add some small rounding here
            endIdx = (currSeg[1] - self._splitListsStart) / self._splitListsStep
            endIdx = int(endIdx - 0.000001)

            outFiles = []
            for idx in range(startIdx, endIdx + 1):
                outFilesTemp = [i for i in self._splitLists[idx] \
                                                            if ifo in i.ifo_list]
                outFiles.extend([i for i in outFilesTemp \
                                      if i.segment_list.intersects_segment(currSeg)])
                # Remove duplicates
                outFiles = list(set(outFiles))

        return self.__class__(outFiles)

    def find_output_with_tag(self, tag):
        """
        Find all files who have tag in self.tags
        """
        # Enforce upper case
        tag = tag.upper()
        return FileList([i for i in self if tag in i.tags])

    def find_output_without_tag(self, tag):
        """
        Find all files who do not have tag in self.tags
        """
        # Enforce upper case
        tag = tag.upper()
        return FileList([i for i in self if not tag in i.tags])

    def find_output_with_ifo(self, ifo):
        """
        Find all files who have ifo = ifo
        """
        # Enforce upper case
        ifo = ifo.upper()
        return FileList([i for i in self if ifo in i.ifo_list])

    def get_times_covered_by_files(self):
        """
        Find the coalesced intersection of the segments of all files in the
        list.
        """
        times = segments.segmentlist([])
        for entry in self:
            times.extend(entry.segment_list)
        times.coalesce()
        return times

    def convert_to_lal_cache(self):
        """
        Return all files in this object as a lal.Cache object
        """
        lal_cache = lal.Cache([])
        for entry in self:
            try:
                lal_cache.append(entry.cache_entry)
            except ValueError:
                pass
        return lal_cache

    def _temporal_split_list(self,numSubLists):
        """
        This internal function is used to speed the code up in cases where a
        number of operations are being made to determine if files overlap a
        specific time. Normally such operations are done on *all* entries with
        *every* call. However, if we predetermine which files are at which
        times, we can avoid testing *every* file every time.
  
        We therefore create numSubLists distinct and equal length time windows
        equally spaced from the first time entry in the list until the last.
        A list is made for each window and files are added to lists which they
        overlap.
 
        If the list changes it should be captured and these split lists become
        invalid. Currently the testing for this is pretty basic
        """
        # Assume segment lists are coalesced!
        startTime = float( min([i.segment_list[0][0] for i in self]))
        endTime = float( max([i.segment_list[-1][-1] for i in self]))
        step = (endTime - startTime) / float(numSubLists)

        # Set up storage
        self._splitLists = []
        for idx in range(numSubLists):
            self._splitLists.append(FileList([]))
        
        # Sort the files

        for ix, currFile in enumerate(self):
            segExtent = currFile.segment_list.extent()
            segExtStart = float(segExtent[0])
            segExtEnd = float(segExtent[1])
            startIdx = (segExtent[0] - startTime) / step
            endIdx = (segExtent[1] - startTime) / step
            # Add some small rounding here
            startIdx = int(startIdx - 0.001) 
            endIdx = int(endIdx + 0.001)

            if startIdx < 0:
                startIdx = 0
            if endIdx >= numSubLists:
                endIdx = numSubLists - 1

            for idx in range(startIdx, endIdx + 1):
                self._splitLists[idx].append(currFile)

        # Set information needed to detect changes and to be used elsewhere
        self._splitListsLength = len(self)
        self._splitListsStart = startTime
        self._splitListsEnd = endTime
        self._splitListsStep = step
        self._splitListsSet = True

    def _check_split_list_validity(self):
        """
        See _temporal_split_list above. This function checks if the current
        split lists are still valid.
        """
        # FIXME: Currently very primitive, but needs to be fast
        if not (hasattr(self,"_splitListsSet") and (self._splitListsSet)):
            return False
        elif len(self) != self._splitListsLength:
            return False
        else:
            return True

    @classmethod
    def load(self, filename):
        """
        Load an AhopeFileList from a pickle file
        """
        f = open(filename, 'r')
        return cPickle.load(f)
    
    def dump(self, filename):
        """
        Output this AhopeFileList to a pickle file
        """
        f = open(filename, 'w')
        cPickle.dump(self, f)
        
    def to_file_object(self, name, out_dir):
        """Dump to a pickle file and return an File object reference of this list
        
        Parameters
        ----------
        name : str
            An identifier of this file. Needs to be unique.
        out_dir : path 
            path to place this file
            
        Returns
        -------
        file : AhopeFile
        """
        make_analysis_dir(out_dir)
        
        file_ref = File('ALL', name, self.get_times_covered_by_files(),
                             extension='.pkl', directory=out_dir)
        self.dump(file_ref.storage_path)
        return file_ref

class OutSegFile(File):
    '''
    This class inherits from the File class, and is designed to store
    workflow output files containing a segment list. This is identical in
    usage to File except for an additional kwarg for holding the
    segment list, if it is known at workflow run time.
    '''
    def __init__(self, ifo, description, segment, fileUrl,
                 segment_list=None, **kwargs):
        """
        See File.__init__ for a full set of documentation for how to
        call this class. The only thing unique and added to this class is
        the required option timeSeg, as described below:

        Parameters:
        ------------
        ifo : string or list (required)
            See File.__init__
        description : string (required)
            See File.__init__
        segment : glue.segments.segment or glue.segments.segmentlist
            See File.__init__
        fileUrl : string (required)
            See File.__init__
            FIXME: This is a kwarg in File and should be here as well,
            if this is removed from the explicit arguments it would allow for
            either fileUrls or constructed file names to be used in File.
        segment_list : glue.segments.segmentlist (optional, default=None)
            A glue.segments.segmentlist covering the times covered by the
            segmentlist associated with this file. If this is the science time
            or CAT_1 file this will be used to determine analysis time. Can
            be added by setting self.segment_list after initializing an instance of
            the class.

        """
        super(OutSegFile, self).__init__(ifo, description, segment, fileUrl,
                              **kwargs)
        self.segmentList = segment_list

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
        newsegment_list = segments.segmentlist()
        for seg in self.segmentList:
            if abs(seg) > minSegLength:
                newsegment_list.append(seg)
        newsegment_list.coalesce()
        self.segmentList = newsegment_list
        self.toSegmentXml()

    def toSegmentXml(self):
        """
        Write the segment list in self.segmentList to the url in self.url.
        """
        from pycbc.events import segments_to_file
        segments_to_file(self.segmentList, self.storage_path, 
                             self.tagged_description,  ifo=self.ifo_string)

def make_external_call(cmdList, out_dir=None, out_basename='external_call',
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
    out_dir : string
        If given the stdout and stderr will be redirected to
        os.path.join(out_dir,out_basename+[".err",".out])
        If not given the stdout and stderr will not be recorded
    out_basename : string
        The value of out_basename used to construct the file names used to
        store stderr and stdout. See out_dir for more information.
    shell : boolean, default=False
        This value will be given as the shell kwarg to the subprocess call.
        **WARNING** See the subprocess documentation for details on this
        Kwarg including a warning about a serious security exploit. Do not
        use this unless you are sure it is necessary **and** safe.
    fail_on_error : boolean, default=True
        If set to true an exception will be raised if the external command does
        not return a code of 0. If set to false such failures will be ignored.
        Stderr and Stdout can be stored in either case using the out_dir
        and out_basename options.

    Returns
    --------
    exitCode : int
        The code returned by the process.
    """
    resolvedExe = resolve_url(cmdList[0])
    if resolvedExe != cmdList[0]:
        os.chmod(resolvedExe, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)
        cmdList[0] = resolvedExe

    if out_dir:
        outBase = os.path.join(out_dir,out_basename)
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
        
def get_random_label():
    """
    Get a random label string to use when clustering jobs.
    """
    return ''.join(random.choice(string.ascii_uppercase + string.digits) \
                   for _ in range(15))


def resolve_url(url):
    """
    Resolves a URL to a local file, and returns the path to
    that file.
    """
    if url.startswith('http://') or url.startswith('https://'):
        filename = url.split('/')[-1]
        filename = os.path.join(os.getcwd(), filename)
        succeeded = False
        num_tries = 5
        t_sleep   = 10

        while not succeeded and num_tries > 0: 
            try:
                response = urllib2.urlopen(url)
                result   = response.read()
                out_file = open(filename, 'w')
                out_file.write(result)
                out_file.close()
                succeeded = True
            except:
                logging.warn("Unable to download %s, retrying" % url)
                time.sleep(t_sleep)
                num_tries -= 1
                t_sleep   *= 2
                
        if not succeeded:
            errMsg  = "Unable to download %s " % (url)
            raise ValueError(errMsg)

    elif url.startswith('file://'):
        filename = url[7:]
    elif url.find('://') != -1:
        # TODO: We could support other schemes such as gsiftp by
        # calling out to globus-url-copy
        errMsg  = "%s: Only supported URL schemes are\n" % (url)
        errMsg += "   file: http: https:" 
        raise ValueError(errMsg)
    else:
        filename = url

    if not os.path.isfile(filename):
        errMsg = "File %s does not exist." %(url)
        raise ValueError(errMsg)

    return filename
   
