# Copyright (C) 2013, 2017  Ian Harry, Alex Nitz, Duncan Brown
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
import sys, os, stat, subprocess, logging, math, string
from six.moves import configparser as ConfigParser
from six.moves import urllib
from six.moves.urllib.request import pathname2url
from six.moves.urllib.parse import urljoin
from six.moves import cPickle
import copy
import numpy, random
from itertools import combinations, groupby, permutations
from operator import attrgetter
from six import string_types
import lal
import lal.utils
import Pegasus.DAX3
from glue import lal as gluelal
from ligo import segments
from glue.ligolw import table, lsctables, ligolw
from glue.ligolw import utils as ligolw_utils
from glue.ligolw.utils import segments as ligolw_segments
from glue.ligolw.utils import process as ligolw_process
from pycbc import makedir
from pycbc.workflow.configuration import WorkflowConfigParser, resolve_url
from pycbc.workflow import pegasus_workflow

class ContentHandler(ligolw.LIGOLWContentHandler):
    pass

lsctables.use_in(ContentHandler)

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
    output, _, _ = check_output_error_and_retcode(*popenargs, **kwargs)
    return output

###############################################################################

def make_analysis_dir(path):
    """
    Make the analysis directory path, any parent directories that don't already
    exist, and the 'logs' subdirectory of path.
    """
    if path is not None:
        makedir(os.path.join(path, 'logs'))

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
    if str(check_output(['nm', '-a', exe_path])).find('condor') != -1:
        return True
    else:
        return False

file_input_from_config_dict = {}

class Executable(pegasus_workflow.Executable):
    # These are the file retention levels
    INTERMEDIATE_PRODUCT = 1
    ALL_TRIGGERS = 2
    MERGED_TRIGGERS = 3
    FINAL_RESULT = 4

    # Set this parameter to indicate that this option is used to specify a
    # file and is *not* handled explicitly in the create_node or __init__
    # methods of the sub-class. Usually that is to say that this option is a
    # file and is normally specified in an file, e.g. a PSD file. As files
    # need to be identified as such to pegasus, this attempts to catch this
    # case.
    # file_input_options = ['--psd-file, '--bank-file'] (as an example)
    file_input_options = []

    # This is the default value. It will give a warning if a class is
    # used where the retention level is not set. The file will still be stored
    KEEP_BUT_RAISE_WARNING = 5
    _warned_classes_list = ['Executable']

    # Sub classes, or instances, should override this. If not overriden the
    # file will be retained, but a warning given
    current_retention_level = KEEP_BUT_RAISE_WARNING
    def __init__(self, cp, name,
                 universe=None, ifos=None, out_dir=None, tags=None):
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
        if isinstance(ifos, string_types):
            self.ifo_list = [ifos]
        else:
            self.ifo_list = ifos
        if self.ifo_list is not None:
            self.ifo_string = ''.join(self.ifo_list)
        else:
            self.ifo_string = None
        self.cp = cp
        self.universe=universe
        self.container_cls = None
        self.container_type = None

        try:
            self.installed = cp.getboolean('pegasus_profile-%s' % name, 'pycbc|installed')
        except:
            self.installed = True

        self.name=name

        self.update_current_tags(tags)

        self.update_output_directory(out_dir=out_dir)

        # Determine the level at which output files should be kept
        self.update_current_retention_level(self.current_retention_level)

        # Determine if this executables should be run in a container
        try:
            self.container_type = cp.get('pegasus_profile-%s' % name,
                                         'container|type')
        except:
            pass

        if self.container_type is not None:
            self.container_img = cp.get('pegasus_profile-%s' % name,
                                        'container|image')
            try:
                self.container_site = cp.get('pegasus_profile-%s' % name,
                                             'container|image_site')
            except:
                self.container_site = 'local'

            try:
                self.container_mount = cp.get('pegasus_profile-%s' % name,
                                             'container|mount').split(',')
            except:
                self.container_mount = None


            self.container_cls = Pegasus.DAX3.Container("{}-container".format(
                                                    name),
                                                    self.container_type,
                                                    self.container_img,
                                                    imagesite=self.container_site,
                                                    mount=self.container_mount)

            super(Executable, self).__init__(self.tagged_name,
                                             installed=self.installed,
                                             container=self.container_cls)

        else:
            super(Executable, self).__init__(self.tagged_name,
                                             installed=self.installed)

        self._set_pegasus_profile_options()

        # Check that the executable actually exists locally or
        # looks like a URL, in which case trust Pegasus to be
        # able to fetch it.
        exe_path = cp.get('executables', name)
        self.needs_fetching = False

        exe_url = urllib.parse.urlparse(exe_path)

        # See if the user specified a list of sites for the executable
        try:
            exe_site_list = cp.get('pegasus_profile-%s' % name, 'pycbc|site')
        except:
            exe_site_list = 'local'

        for s in exe_site_list.split(','):
            exe_site = s.strip()

            if exe_url.scheme in ['', 'file']:
                if exe_site is 'local':
                    # Check that executables at file urls
                    #  on the local site exist
                    if os.path.isfile(exe_url.path) is False:
                        raise TypeError("Failed to find %s executable "
                                        "at %s on site %s" % (name, exe_path,
                                        exe_site))
            else:
                # Could be http, gsiftp, etc. so it needs fetching if run now
                self.needs_fetching = True

            self.add_pfn(exe_path, site=exe_site)
            logging.debug("Using %s executable "
                          "at %s on site %s" % (name, exe_url.path, exe_site))

        # Determine the condor universe if we aren't given one
        if self.universe is None:
            if is_condor_exec(exe_path):
                self.universe = 'standard'
            else:
                self.universe = 'vanilla'

        logging.debug("%s executable will run as %s universe"
                     % (name, self.universe))

        self.set_universe(self.universe)

        if hasattr(self, "group_jobs"):
            self.add_profile('pegasus', 'clusters.size', self.group_jobs)

    @property
    def ifo(self):
        """Return the ifo.

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
        """Add profile from configuration file.

        Parameters
        -----------
        cp : ConfigParser object
            The ConfigParser object holding the workflow configuration settings
        sec : string
            The section containing options for this job.
        """
        for opt in cp.options(sec):
            namespace = opt.split('|')[0]
            if namespace == 'pycbc' or namespace == 'container':
                continue

            value = cp.get(sec, opt).strip()
            key = opt.split('|')[1]
            self.add_profile(namespace, key, value, force=True)

            # Remove if Pegasus can apply this hint in the TC
            if namespace == 'hints' and key == 'execution.site':
                self.execution_site = value

    def add_ini_opts(self, cp, sec):
        """Add job-specific options from configuration file.

        Parameters
        -----------
        cp : ConfigParser object
            The ConfigParser object holding the workflow configuration settings
        sec : string
            The section containing options for this job.
        """
        for opt in cp.options(sec):
            value = cp.get(sec, opt).strip()
            opt = '--%s' %(opt,)
            if opt in self.file_input_options:
                # This now expects the option to be a file
                # Check is we have a list of files
                values = [path for path in value.split(' ') if path]

                self.common_raw_options.append(opt)
                self.common_raw_options.append(' ')

                # Get LFN and PFN
                for path in values:
                    # Here I decide if the path is URL or
                    # IFO:/path/to/file or IFO:url://path/to/file
                    # That's somewhat tricksy as we used : as delimiter
                    split_path = path.split(':', 1)
                    if len(split_path) == 1:
                        ifo = None
                        path = path
                    else:
                        # Have I split a URL or not?
                        if split_path[1].startswith('//'):
                            # URL
                            ifo = None
                            path = path
                        else:
                            #IFO:path or IFO:URL
                            ifo = split_path[0]
                            path = split_path[1]

                    curr_lfn = os.path.basename(path)

                    # If the file exists make sure to use the
                    # fill path as a file:// URL
                    if os.path.isfile(path):
                        curr_pfn = urljoin('file:',
                                           pathname2url(os.path.abspath(path)))
                    else:
                        curr_pfn = path

                    if curr_lfn in file_input_from_config_dict.keys():
                        file_pfn = file_input_from_config_dict[curr_lfn][2]
                        assert(file_pfn == curr_pfn)
                        curr_file = file_input_from_config_dict[curr_lfn][1]
                    else:
                        local_file_path = resolve_url(curr_pfn)
                        curr_file = File.from_path(local_file_path)
                        tuple_val = (local_file_path, curr_file, curr_pfn)
                        file_input_from_config_dict[curr_lfn] = tuple_val
                    self.common_input_files.append(curr_file)
                    if ifo:
                        self.common_raw_options.append(ifo + ':')
                        self.common_raw_options.append(curr_file.dax_repr)
                    else:
                        self.common_raw_options.append(curr_file.dax_repr)
                    self.common_raw_options.append(' ')
            else:
                self.common_options += [opt, value]

    def add_opt(self, opt, value=None):
        """Add option to job.

        Parameters
        -----------
        opt : string
            Name of option (e.g. --output-file-format)
        value : string, (default=None)
            The value for the option (no value if set to None).
        """
        if value is None:
            self.common_options += [opt]
        else:
            self.common_options += [opt, value]

    def get_opt(self, opt):
        """Get value of option from configuration file

        Parameters
        -----------
        opt : string
            Name of option (e.g. output-file-format)

        Returns
        --------
        value : string
            The value for the option. Returns None if option not present.
        """
        for sec in self.sections:
            try:
                key = self.cp.get(sec, opt)
                if key:
                    return key
            except ConfigParser.NoOptionError:
                pass

        return None

    def has_opt(self, opt):
        """Check if option is present in configuration file

        Parameters
        -----------
        opt : string
            Name of option (e.g. output-file-format)
        """
        for sec in self.sections:
            val = self.cp.has_option(sec, opt)
            if val:
                return val

        return False

    def create_node(self):
        """Default node constructor.

        This is usually overridden by subclasses of Executable.
        """
        return Node(self)

    def update_current_retention_level(self, value):
        """Set a new value for the current retention level.

        This updates the value of self.retain_files for an updated value of the
        retention level.

        Parameters
        -----------
        value : int
            The new value to use for the retention level.
        """
        # Determine the level at which output files should be kept
        self.current_retention_level = value
        try:
            global_retention_level = \
                self.cp.get_opt_tags("workflow", "file-retention-level",
                                   self.tags+[self.name])
        except ConfigParser.Error:
            msg="Cannot find file-retention-level in [workflow] section "
            msg+="of the configuration file. Setting a default value of "
            msg+="retain all files."
            logging.warn(msg)
            self.retain_files = True
            self.global_retention_threshold = 1
            self.cp.set("workflow", "file-retention-level", "all_files")
        else:
            # FIXME: Are these names suitably descriptive?
            retention_choices = {
                                 'all_files' : 1,
                                 'all_triggers' : 2,
                                 'merged_triggers' : 3,
                                 'results' : 4
                                }
            try:
                self.global_retention_threshold = \
                      retention_choices[global_retention_level]
            except KeyError:
                err_msg = "Cannot recognize the file-retention-level in the "
                err_msg += "[workflow] section of the ini file. "
                err_msg += "Got : {0}.".format(global_retention_level)
                err_msg += "Valid options are: 'all_files', 'all_triggers',"
                err_msg += "'merged_triggers' or 'results' "
                raise ValueError(err_msg)
            if self.current_retention_level == 5:
                self.retain_files = True
                if type(self).__name__ in Executable._warned_classes_list:
                    pass
                else:
                    warn_msg = "Attribute current_retention_level has not "
                    warn_msg += "been set in class {0}. ".format(type(self))
                    warn_msg += "This value should be set explicitly. "
                    warn_msg += "All output from this class will be stored."
                    logging.warn(warn_msg)
                    Executable._warned_classes_list.append(type(self).__name__)
            elif self.global_retention_threshold > self.current_retention_level:
                self.retain_files = False
            else:
                self.retain_files = True

    def update_current_tags(self, tags):
        """Set a new set of tags for this executable.

        Update the set of tags that this job will use. This updated default
        file naming and shared options. It will *not* update the pegasus
        profile, which belong to the executable and cannot be different for
        different nodes.

        Parameters
        -----------
        tags : list
            The new list of tags to consider.
        """
        if tags is None:
            tags = []
        tags = [tag.upper() for tag in tags]
        self.tags = tags

        if len(tags) > 6:
            warn_msg = "This job has way too many tags. "
            warn_msg += "Current tags are {}. ".format(' '.join(tags))
            warn_msg += "Current executable {}.".format(self.name)
            logging.info(warn_msg)

        if len(tags) != 0:
            self.tagged_name = "{0}-{1}".format(self.name, '_'.join(tags))
        else:
            self.tagged_name = self.name
        if self.ifo_string is not None:
            self.tagged_name = "{0}-{1}".format(self.tagged_name,
                                                self.ifo_string)


        # Determine the sections from the ini file that will configure
        # this executable
        sections = [self.name]
        if self.ifo_list is not None:
            if len(self.ifo_list) > 1:
                sec_tags = tags + self.ifo_list + [self.ifo_string]
            else:
                sec_tags = tags + self.ifo_list
        else:
            sec_tags = tags
        for sec_len in range(1, len(sec_tags)+1):
            for tag_permutation in permutations(sec_tags, sec_len):
                joined_name = '-'.join(tag_permutation)
                section = '{0}-{1}'.format(self.name, joined_name.lower())
                if self.cp.has_section(section):
                    sections.append(section)

        self.sections = sections

        # Do some basic sanity checking on the options
        for sec1, sec2 in combinations(sections, 2):
            self.cp.check_duplicate_options(sec1, sec2, raise_error=True)

        # collect the options and profile information
        # from the ini file section(s)
        self.common_options = []
        self.common_raw_options = []
        self.common_input_files = []
        for sec in sections:
            if self.cp.has_section(sec):
                self.add_ini_opts(self.cp, sec)
            else:
                warn_string = "warning: config file is missing section "
                warn_string += "[{0}]".format(sec)
                logging.warn(warn_string)

    def update_output_directory(self, out_dir=None):
        """Update the default output directory for output files.

        Parameters
        -----------
        out_dir : string (optional, default=None)
            If provided use this as the output directory. Else choose this
            automatically from the tags.
        """
        # Determine the output directory
        if out_dir is not None:
            self.out_dir = out_dir
        elif len(self.tags) == 0:
            self.out_dir = self.name
        else:
            self.out_dir = self.tagged_name
        if not os.path.isabs(self.out_dir):
            self.out_dir = os.path.join(os.getcwd(), self.out_dir)

    def _set_pegasus_profile_options(self):
        """Set the pegasus-profile settings for this Executable.

        These are a property of the Executable and not of nodes that it will
        spawn. Therefore it *cannot* be updated without also changing values
        for nodes that might already have been created. Therefore this is
        only called once in __init__. Second calls to this will fail.
        """
        # Add executable non-specific profile information
        if self.cp.has_section('pegasus_profile'):
            self.add_ini_profile(self.cp, 'pegasus_profile')

        # Executable- and tag-specific profile information
        for sec in self.sections:
            if self.cp.has_section('pegasus_profile-{0}'.format(sec)):
                self.add_ini_profile(self.cp,
                                     'pegasus_profile-{0}'.format(sec))

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
        if self.in_workflow is not False:
            name = self.name + '.map'
        else:
            name = 'output.map'
        path =  os.path.join(os.getcwd(), name)
        return path

    @property
    def transformation_catalog(self):
        if self.in_workflow is not False:
            name = self.name + '.tc.txt'
        else:
            name = 'tc.txt'
        path =  os.path.join(os.getcwd(), name)
        return path

    @property
    def staging_site(self):
        if self.in_workflow is not False:
            workflow_section = 'workflow-%s' % self.name
        else:
            workflow_section = 'workflow'
        try:
            staging_site = self.cp.get(workflow_section,'staging-site')
        except:
            staging_site = None
        return staging_site

    def execute_node(self, node, verbatim_exe = False):
        """ Execute this node immediately on the local machine
        """
        node.executed = True

        # Check that the PFN is for a file or path
        if node.executable.needs_fetching:
            try:
                # The pfn may have been marked local...
                pfn = node.executable.get_pfn()
            except:
                # or it may have been marked nonlocal.  That's
                # fine, we'll resolve the URL and make a local
                # entry.
                pfn = node.executable.get_pfn('nonlocal')

            resolved = resolve_url(pfn, permissions=stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)
            node.executable.clear_pfns()
            node.executable.add_pfn(urljoin('file:', pathname2url(resolved)),
                                    site='local')

        cmd_list = node.get_command_line()

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
            fil.PFN(urljoin('file:', pathname2url(fil.storage_path)),
                    site='local')

    @staticmethod
    def set_job_properties(job, output_map_file, transformation_catalog_file,
                           staging_site=None):

        job.addArguments('-Dpegasus.dir.storage.mapper.replica.file=%s' %
                         os.path.basename(output_map_file.name))
        job.uses(output_map_file, link=Pegasus.DAX3.Link.INPUT)
        job.addArguments('-Dpegasus.dir.storage.mapper.replica=File')

        # FIXME this is an ugly hack to connect the right transformation
        # catalog to the right DAX beacuse Pegasus 4.9 does not support
        # the full transformation catalog syntax in the DAX. This will go
        # away in Pegasus 5.x when this code is re-written.

        job.addArguments('-Dpegasus.catalog.transformation.file=%s' %
                         os.path.basename(transformation_catalog_file.name))
        job.uses(transformation_catalog_file, link=Pegasus.DAX3.Link.INPUT)

        job.addArguments('--output-site local')
        job.addArguments('--cleanup inplace')
        job.addArguments('--cluster label,horizontal')
        job.addArguments('-vvv')

        # FIXME _reuse_cache needs to be fixed to use PFNs properly. This will
        # work as pegasus-plan is currently invoked on the local site so has
        # access to a file in os.getcwd() but this code is fragile.
        job.addArguments('--cache %s' % os.path.join(os.getcwd(), '_reuse.cache'))

        if staging_site:
            job.addArguments('--staging-site %s' % staging_site)

    def save(self, filename=None, output_map_path=None,
             transformation_catalog_path=None, staging_site=None):

        if output_map_path is None:
            output_map_path = self.output_map
        output_map_file = Pegasus.DAX3.File(os.path.basename(output_map_path))
        output_map_file.addPFN(Pegasus.DAX3.PFN(output_map_path, 'local'))
        if self.in_workflow is not False:
            self.in_workflow._adag.addFile(output_map_file)

        if transformation_catalog_path is None:
            transformation_catalog_path = self.transformation_catalog
        transformation_catalog_file = Pegasus.DAX3.File(os.path.basename(
                                                        transformation_catalog_path))
        transformation_catalog_file.addPFN(Pegasus.DAX3.PFN(
            transformation_catalog_path, 'local'))
        if self.in_workflow is not False:
            self.in_workflow._adag.addFile(transformation_catalog_file)

        if staging_site is None:
            staging_site = self.staging_site

        Workflow.set_job_properties(self.as_job, output_map_file,
                                    transformation_catalog_file,
                                    staging_site)

        # add executable pfns for local site to dax
        for exe in self._executables:
            exe.insert_into_dax(self._adag)

        # add workflow input files pfns for local site to dax
        for fil in self._inputs:
            fil.insert_into_dax(self._adag)

        # save the configuration file
        ini_file = os.path.abspath(self.name + '.ini')

        # This shouldn't already exist, but just in case
        if os.path.isfile(ini_file):
            err_msg = "Refusing to overwrite configuration file that "
            err_msg += "shouldn't be there: "
            err_msg += os.path.join(os.getcwd(), ini_file)
            raise ValueError(err_msg)

        fp = open(ini_file, 'w')
        self.cp.write(fp)
        fp.close()

        # save the dax file
        super(Workflow, self).save(filename=filename,
                                   tc=transformation_catalog_path)

        # add workflow storage locations to the output mapper
        f = open(output_map_path, 'w')
        for out in self._outputs:
            try:
                f.write(out.output_map_str() + '\n')
            except ValueError:
                # There was no storage path
                pass

    def save_config(self, fname, output_dir, cp=None):
        """ Writes configuration file to disk and returns a pycbc.workflow.File
        instance for the configuration file.

        Parameters
        -----------
        fname : string
            The filename of the configuration file written to disk.
        output_dir : string
            The directory where the file is written to disk.
        cp : ConfigParser object
            The ConfigParser object to write. If None then uses self.cp.

        Returns
        -------
        FileList
            The FileList object with the configuration file.
        """
        cp = self.cp if cp is None else cp
        ini_file_path = os.path.join(output_dir, fname)
        with open(ini_file_path, "wb") as fp:
            cp.write(fp)
        ini_file = FileList([File(self.ifos, "",
                                  self.analysis_time,
                                  file_url="file://" + ini_file_path)])
        return ini_file

class Node(pegasus_workflow.Node):
    def __init__(self, executable):
        super(Node, self).__init__(executable)
        self.executed = False
        self.set_category(executable.name)

        if executable.universe == 'vanilla' and executable.installed:
            self.add_profile('condor', 'getenv', 'True')

        if hasattr(executable, 'execution_site'):
            self.add_profile('hints', 'execution.site', executable.execution_site)

        self._options += self.executable.common_options
        self._raw_options += self.executable.common_raw_options
        for inp in self.executable.common_input_files:
            self._add_input(inp)

    def get_command_line(self):
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
        exe_path = urllib.parse.urlsplit(self.executable.get_pfn()).path

        return [exe_path] + arglist

    def new_output_file_opt(self, valid_seg, extension, option_name, tags=None,
                            store_file=None, use_tmp_subdirs=False):
        """
        This function will create a workflow.File object corresponding to the given
        information and then add that file as output of this node.

        Parameters
        -----------
        valid_seg : ligo.segments.segment
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
        if tags is None:
            tags = []

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

    def add_multiifo_input_list_opt(self, opt, inputs):
        """ Add an option that determines a list of inputs from multiple
            detectors. Files will be supplied as --opt ifo1:input1 ifo2:input2
            .....
        """
        # NOTE: Here we have to use the raw arguments functionality as the
        #       file and ifo are not space separated.
        self.add_raw_arg(opt)
        self.add_raw_arg(' ')
        for infile in inputs:
            self.add_raw_arg(infile.ifo)
            self.add_raw_arg(':')
            self.add_raw_arg(infile.name)
            self.add_raw_arg(' ')
            self._add_input(infile)

    def add_multiifo_output_list_opt(self, opt, outputs):
        """ Add an option that determines a list of outputs from multiple
            detectors. Files will be supplied as --opt ifo1:input1 ifo2:input2
            .....
        """
        # NOTE: Here we have to use the raw arguments functionality as the
        #       file and ifo are not space separated.
        self.add_raw_arg(opt)
        self.add_raw_arg(' ')
        for outfile in outputs:
            self.add_raw_arg(outfile.ifo)
            self.add_raw_arg(':')
            self.add_raw_arg(outfile.name)
            self.add_raw_arg(' ')
            self._add_output(outfile)

    def new_multiifo_output_list_opt(self, opt, ifos, analysis_time, extension,
                                     tags=None, store_file=None,
                                     use_tmp_subdirs=False):
        """ Add an option that determines a list of outputs from multiple
            detectors. Files will be supplied as --opt ifo1:input1 ifo2:input2
            .....
            File names are created internally from the provided extension and
            analysis time.
        """
        if tags is None:
            tags = []
        all_tags = copy.deepcopy(self.executable.tags)
        for tag in tags:
            if tag not in all_tags:
                all_tags.append(tag)

        output_files = FileList([])
        store_file = store_file if store_file is not None \
                                              else self.executable.retain_files

        for ifo in ifos:
            curr_file = File(ifo, self.executable.name, analysis_time,
                             extension=extension, store_file=store_file,
                             directory=self.executable.out_dir, tags=all_tags,
                             use_tmp_subdirs=use_tmp_subdirs)
            output_files.append(curr_file)
        self.add_multiifo_output_list_opt(opt, output_files)


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

    >> c = File("H1", "INSPIRAL_S6LOWMASS", segments.segment(815901601, 815902001), file_url="file://localhost/home/spxiwh/H1-INSPIRAL_S6LOWMASS-815901601-400.xml.gz" )

    another where the file url is generated from the inputs:

    >> c = File("H1", "INSPIRAL_S6LOWMASS", segments.segment(815901601, 815902001), directory="/home/spxiwh", extension="xml.gz" )

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
        if isinstance(ifos, string_types):
            self.ifo_list = [ifos]
        else:
            self.ifo_list = ifos
        self.ifo_string = ''.join(self.ifo_list)
        self.description = exe_name

        if isinstance(segs, segments.segment):
            self.segment_list = segments.segmentlist([segs])
        elif isinstance(segs, (segments.segmentlist)):
            self.segment_list = segs
        else:
            err = "segs input must be either ligo.segments.segment or "
            err += "segments.segmentlist. Got %s." %(str(type(segs)),)
            raise ValueError(err)
        if tags is not None:
            self.tags = [t.upper() for t in tags]
        else:
            self.tags = []
        if len(self.tags):
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
            file_url = urllib.parse.urlunparse(['file', 'localhost', path,
                                                None, None, None])

        # Let's do a test here
        if use_tmp_subdirs and len(self.segment_list):
            pegasus_lfn = str(int(self.segment_list.extent()[0]))[:-4]
            pegasus_lfn = pegasus_lfn + '/' + os.path.basename(file_url)
        else:
            pegasus_lfn = os.path.basename(file_url)
        super(File, self).__init__(pegasus_lfn)

        if store_file:
            self.storage_path = urllib.parse.urlsplit(file_url).path
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

        file_url = urllib.parse.urlunparse(['file', 'localhost',
                                            self.storage_path, None,
                                            None, None])
        cache_entry = lal.utils.CacheEntry(self.ifo_string,
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
        '''Returns one File most appropriate at the given time/time range.

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
        pycbc_file : pycbc.workflow.File instance
           The File that corresponds to the time or time range
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
        currsegment_list = segments.segmentlist([segments.segment(start, end)])

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
                if idx < 0 or idx >= self._splitListsNum:
                    continue
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
        return FileList([i for i in self if tag not in i.tags])

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
        Return all files in this object as a glue.lal.Cache object
        """
        lal_cache = gluelal.Cache([])
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

        for currFile in self:
            segExtent = currFile.segment_list.extent()
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
        self._splitListsNum = numSubLists
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
    def load(cls, filename):
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

class SegFile(File):
    '''
    This class inherits from the File class, and is designed to store
    workflow output files containing a segment dict. This is identical in
    usage to File except for an additional kwarg for holding the
    segment dictionary, if it is known at workflow run time.
    '''
    def __init__(self, ifo_list, description, valid_segment,
                 segment_dict=None, seg_summ_dict=None, **kwargs):
        """
        See File.__init__ for a full set of documentation for how to
        call this class. The only thing unique and added to this class is
        the optional segment_dict. NOTE that while segment_dict is a
        ligo.segments.segmentlistdict rather than the usual dict[ifo]
        we key by dict[ifo:name].

        Parameters
        ------------
        ifo_list : string or list (required)
            See File.__init__
        description : string (required)
            See File.__init__
        segment : ligo.segments.segment or ligo.segments.segmentlist
            See File.__init__
        segment_dict : ligo.segments.segmentlistdict (optional, default=None)
            A ligo.segments.segmentlistdict covering the times covered by the
            segmentlistdict associated with this file.
            Can be added by setting self.segment_dict after initializing an
            instance of the class.

        """
        super(SegFile, self).__init__(ifo_list, description, valid_segment,
                                      **kwargs)
        # To avoid confusion with the segment_list property of the parent class
        # we refer to this as valid_segments here
        self.valid_segments = self.segment_list
        self.segment_dict = segment_dict
        self.seg_summ_dict = seg_summ_dict

    @classmethod
    def from_segment_list(cls, description, segmentlist, name, ifo,
                          seg_summ_list=None, **kwargs):
        """ Initialize a SegFile object from a segmentlist.

        Parameters
        ------------
        description : string (required)
            See File.__init__
        segmentlist : ligo.segments.segmentslist
            The segment list that will be stored in this file.
        name : str
            The name of the segment lists to be stored in the file.
        ifo : str
            The ifo of the segment lists to be stored in this file.
        seg_summ_list : ligo.segments.segmentslist (OPTIONAL)
            Specify the segment_summary segmentlist that goes along with the
            segmentlist. Default=None, in this case segment_summary is taken
            from the valid_segment of the SegFile class.
        """
        seglistdict = segments.segmentlistdict()
        seglistdict[ifo + ':' + name] = segmentlist
        if seg_summ_list is not None:
            seg_summ_dict = segments.segmentlistdict()
            seg_summ_dict[ifo + ':' + name] = seg_summ_list
        else:
            seg_summ_dict = None
        return cls.from_segment_list_dict(description, seglistdict,
                                          seg_summ_dict=None, **kwargs)

    @classmethod
    def from_multi_segment_list(cls, description, segmentlists, names, ifos,
                                seg_summ_lists=None, **kwargs):
        """ Initialize a SegFile object from a list of segmentlists.

        Parameters
        ------------
        description : string (required)
            See File.__init__
        segmentlists : List of ligo.segments.segmentslist
            List of segment lists that will be stored in this file.
        names : List of str
            List of names of the segment lists to be stored in the file.
        ifos : str
            List of ifos of the segment lists to be stored in this file.
        seg_summ_lists : ligo.segments.segmentslist (OPTIONAL)
            Specify the segment_summary segmentlists that go along with the
            segmentlists. Default=None, in this case segment_summary is taken
            from the valid_segment of the SegFile class.
        """
        seglistdict = segments.segmentlistdict()
        for name, ifo, segmentlist in zip(names, ifos, segmentlists):
            seglistdict[ifo + ':' + name] = segmentlist
        if seg_summ_lists is not None:
            seg_summ_dict = segments.segmentlistdict()
            for name, ifo, seg_summ_list in zip(names, ifos, seg_summ_lists):
                seg_summ_dict[ifo + ':' + name] = seg_summ_list
        else:
            seg_summ_dict = None

        return cls.from_segment_list_dict(description, seglistdict,
                                         seg_summ_dict=seg_summ_dict, **kwargs)

    @classmethod
    def from_segment_list_dict(cls, description, segmentlistdict,
                               ifo_list=None, valid_segment=None,
                               file_exists=False, seg_summ_dict=None,
                               **kwargs):
        """ Initialize a SegFile object from a segmentlistdict.

        Parameters
        ------------
        description : string (required)
            See File.__init__
        segmentlistdict : ligo.segments.segmentslistdict
            See SegFile.__init__
        ifo_list : string or list (optional)
            See File.__init__, if not given a list of all ifos in the
            segmentlistdict object will be used
        valid_segment : ligo.segments.segment or ligo.segments.segmentlist
            See File.__init__, if not given the extent of all segments in the
            segmentlistdict is used.
        file_exists : boolean (default = False)
            If provided and set to True it is assumed that this file already
            exists on disk and so there is no need to write again.
        seg_summ_dict : ligo.segments.segmentslistdict
            Optional. See SegFile.__init__.
        """
        if ifo_list is None:
            ifo_set = set([i.split(':')[0] for i in segmentlistdict.keys()])
            ifo_list = list(ifo_set)
            ifo_list.sort()
        if valid_segment is None:
            if seg_summ_dict and \
                    numpy.any([len(v) for _, v in seg_summ_dict.items()]):
                # Only come here if seg_summ_dict is supplied and it is
                # not empty.
                valid_segment = seg_summ_dict.extent_all()
            else:
                try:
                    valid_segment = segmentlistdict.extent_all()
                except:
                    # Numpty probably didn't supply a glue.segmentlistdict
                    segmentlistdict=segments.segmentlistdict(segmentlistdict)
                    try:
                        valid_segment = segmentlistdict.extent_all()
                    except ValueError:
                        # No segment_summary and segment list is empty
                        # Setting valid segment now is hard!
                        warn_msg = "No information with which to set valid "
                        warn_msg += "segment."
                        logging.warn(warn_msg)
                        valid_segment = segments.segment([0,1])
        instnc = cls(ifo_list, description, valid_segment,
                     segment_dict=segmentlistdict, seg_summ_dict=seg_summ_dict,
                     **kwargs)
        if not file_exists:
            instnc.to_segment_xml()
        else:
            instnc.PFN(urljoin('file:', pathname2url(instnc.storage_path)),
                       site='local')
        return instnc

    @classmethod
    def from_segment_xml(cls, xml_file, **kwargs):
        """
        Read a ligo.segments.segmentlist from the file object file containing an
        xml segment table.

        Parameters
        -----------
        xml_file : file object
            file object for segment xml file
        """
        # load xmldocument and SegmentDefTable and SegmentTables
        fp = open(xml_file, 'rb')
        xmldoc, _ = ligolw_utils.load_fileobj(fp,
                                              gz=xml_file.endswith(".gz"),
                                              contenthandler=ContentHandler)

        seg_def_table = table.get_table(xmldoc,
                                        lsctables.SegmentDefTable.tableName)
        seg_table = table.get_table(xmldoc, lsctables.SegmentTable.tableName)
        seg_sum_table = table.get_table(xmldoc,
                                        lsctables.SegmentSumTable.tableName)

        segs = segments.segmentlistdict()
        seg_summ = segments.segmentlistdict()

        seg_id = {}
        for seg_def in seg_def_table:
            # Here we want to encode ifo and segment name
            full_channel_name = ':'.join([str(seg_def.ifos),
                                          str(seg_def.name)])
            seg_id[int(seg_def.segment_def_id)] = full_channel_name
            segs[full_channel_name] = segments.segmentlist()
            seg_summ[full_channel_name] = segments.segmentlist()

        for seg in seg_table:
            seg_obj = segments.segment(
                    lal.LIGOTimeGPS(seg.start_time, seg.start_time_ns),
                    lal.LIGOTimeGPS(seg.end_time, seg.end_time_ns))
            segs[seg_id[int(seg.segment_def_id)]].append(seg_obj)

        for seg in seg_sum_table:
            seg_obj = segments.segment(
                    lal.LIGOTimeGPS(seg.start_time, seg.start_time_ns),
                    lal.LIGOTimeGPS(seg.end_time, seg.end_time_ns))
            seg_summ[seg_id[int(seg.segment_def_id)]].append(seg_obj)

        for seg_name in seg_id.values():
            segs[seg_name] = segs[seg_name].coalesce()

        xmldoc.unlink()
        fp.close()
        curr_url = urllib.parse.urlunparse(['file', 'localhost', xml_file,
                                            None, None, None])

        return cls.from_segment_list_dict('SEGMENTS', segs, file_url=curr_url,
                                          file_exists=True,
                                          seg_summ_dict=seg_summ, **kwargs)

    def remove_short_sci_segs(self, minSegLength):
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
        for key, seglist in self.segment_dict.items():
            newsegment_list = segments.segmentlist()
            for seg in seglist:
                if abs(seg) > minSegLength:
                    newsegment_list.append(seg)
            newsegment_list.coalesce()
            self.segment_dict[key] = newsegment_list
        self.to_segment_xml(override_file_if_exists=True)

    def return_union_seglist(self):
        return self.segment_dict.union(self.segment_dict.keys())

    def parse_segdict_key(self, key):
        """
        Return ifo and name from the segdict key.
        """
        splt = key.split(':')
        if len(splt) == 2:
            return splt[0], splt[1]
        else:
            err_msg = "Key should be of the format 'ifo:name', got %s." %(key,)
            raise ValueError(err_msg)

    def to_segment_xml(self, override_file_if_exists=False):
        """
        Write the segment list in self.segmentList to self.storage_path.
        """
        # create XML doc and add process table
        outdoc = ligolw.Document()
        outdoc.appendChild(ligolw.LIGO_LW())
        process = ligolw_process.register_to_xmldoc(outdoc, sys.argv[0], {})

        for key, seglist in self.segment_dict.items():
            ifo, name = self.parse_segdict_key(key)
            # Ensure we have LIGOTimeGPS
            fsegs = [(lal.LIGOTimeGPS(seg[0]),
                      lal.LIGOTimeGPS(seg[1])) for seg in seglist]

            if self.seg_summ_dict is None:
                vsegs = [(lal.LIGOTimeGPS(seg[0]),
                          lal.LIGOTimeGPS(seg[1])) \
                         for seg in self.valid_segments]
            else:
                vsegs = [(lal.LIGOTimeGPS(seg[0]),
                          lal.LIGOTimeGPS(seg[1])) \
                         for seg in self.seg_summ_dict[key]]

            # Add using glue library to set all segment tables
            with ligolw_segments.LigolwSegments(outdoc, process) as x:
                x.add(ligolw_segments.LigolwSegmentList(active=fsegs,
                                    instruments=set([ifo]), name=name,
                                    version=1, valid=vsegs))

        # write file
        url = urljoin('file:', pathname2url(self.storage_path))
        if not override_file_if_exists or not self.has_pfn(url, site='local'):
            self.PFN(url, site='local')
        ligolw_utils.write_filename(outdoc, self.storage_path)


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
    science_segs : ifo-keyed dictionary of ligo.segments.segmentlist instances
        The list of times that are being analysed in this workflow.

    Returns
    --------
    fullSegment : ligo.segments.segment
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
