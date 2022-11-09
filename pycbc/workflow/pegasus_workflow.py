# Copyright (C) 2014  Alex Nitz
#
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
""" This module provides thin wrappers around Pegasus.DAX3 functionality that
provides additional abstraction and argument handling.
"""
import os
import shutil
import tempfile
from urllib.request import pathname2url
from urllib.parse import urljoin, urlsplit
import Pegasus.api as dax

PEGASUS_FILE_DIRECTORY = os.path.join(os.path.dirname(__file__),
                                      'pegasus_files')

GRID_START_TEMPLATE = '''#!/bin/bash

if [ -f /tmp/x509up_u`id -u` ] ; then
  unset X509_USER_PROXY
else
  if [ ! -z ${X509_USER_PROXY} ] ; then
    if [ -f ${X509_USER_PROXY} ] ; then
      cp -a ${X509_USER_PROXY} /tmp/x509up_u`id -u`
    fi
  fi
  unset X509_USER_PROXY
fi

# Check that the proxy is valid
grid-proxy-info -exists
RESULT=${?}
if [ ${RESULT} -eq 0 ] ; then
  PROXY_TYPE=`grid-proxy-info -type | tr -d ' '`
  if [ x${PROXY_TYPE} == 'xRFC3820compliantimpersonationproxy' ] ; then
    grid-proxy-info
  else
    cp /tmp/x509up_u`id -u` /tmp/x509up_u`id -u`.orig
    grid-proxy-init -cert /tmp/x509up_u`id -u`.orig -key /tmp/x509up_u`id -u`.orig
    rm -f /tmp/x509up_u`id -u`.orig
    grid-proxy-info
  fi
else
  echo "Error: Could not find a valid grid proxy to submit workflow."
  exit 1
fi
'''

class ProfileShortcuts(object):
    """ Container of common methods for setting pegasus profile information
    on Executables and nodes. This class expects to be inherited from
    and for a add_profile method to be implemented.
    """
    def set_memory(self, size):
        """ Set the amount of memory that is required in megabytes
        """
        self.add_profile('condor', 'request_memory', '%sM' % size)

    def set_storage(self, size):
        """ Set the amount of storage required in megabytes
        """
        self.add_profile('condor', 'request_disk', '%sM' % size)

    def set_num_cpus(self, number):
        self.add_profile('condor', 'request_cpus', number)

    def set_universe(self, universe):
        if universe == 'standard':
            self.add_profile("pegasus", "gridstart", "none")

        self.add_profile("condor", "universe", universe)

    def set_category(self, category):
        self.add_profile('dagman', 'category', category)

    def set_priority(self, priority):
        self.add_profile('dagman', 'priority', priority)

    def set_num_retries(self, number):
        self.add_profile("dagman", "retry", number)

    def set_execution_site(self, site):
        self.add_profile("selector", "execution_site", site)


class Executable(ProfileShortcuts):
    """ The workflow representation of an Executable
    """
    id = 0
    def __init__(self, name, os='linux',
                 arch='x86_64', installed=False,
                 container=None):
        self.logical_name = name + "_ID%s" % str(Executable.id)
        self.pegasus_name = name
        Executable.id += 1
        self.os = dax.OS(os)
        self.arch = dax.Arch(arch)
        self.installed = installed
        self.container = container
        self.in_workflow = False
        self.profiles = {}
        self.transformations = {}

    def create_transformation(self, site, url):
        transform = Transformation(
            self.logical_name,
            site=site,
            pfn=url,
            is_stageable=self.installed,
            arch=self.arch,
            os_type=self.os,
            container=self.container
        )
        transform.pycbc_name = self.pegasus_name
        for (namespace, key), value in self.profiles.items():
            transform.add_profiles(
                dax.Namespace(namespace),
                key=key,
                value=value
            )
        self.transformations[site] = transform

    def add_profile(self, namespace, key, value):
        """ Add profile information to this executable
        """
        if self.transformations:
            err_msg = "Need code changes to be able to add profiles "
            err_msg += "after transformations are created."
            raise ValueError(err_msg)
        self.profiles[(namespace, key)] = value


class Transformation(dax.Transformation):

    def is_same_as(self, other):
        test_vals = ['namespace', 'version']
        test_site_vals = ['arch', 'os_type', 'os_release',
                          'os_version', 'bypass', 'container']
        # Check for logical name first
        if not self.pycbc_name == other.pycbc_name:
            return False

        # Check the properties of the executable
        for val in test_vals:
            sattr = getattr(self, val)
            oattr = getattr(other, val)
            if not sattr == oattr:
                return False
        # Some properties are stored in the TransformationSite
        self_site = list(self.sites.values())
        assert len(self_site) == 1
        self_site = self_site[0]
        other_site = list(other.sites.values())
        assert len(other_site) == 1
        other_site = other_site[0]
        for val in test_site_vals:
            sattr = getattr(self_site, val)
            oattr = getattr(other_site, val)
            if not sattr == oattr:
                return False

        # Also check the "profile". This is things like Universe, RAM/disk/CPU
        # requests, execution site, getenv=True, etc.
        for profile in self.profiles:
            if profile not in other.profiles:
                return False
        for profile in other.profiles:
            if profile not in self.profiles:
                return False

        return True


class Node(ProfileShortcuts):
    def __init__(self, transformation):
        self.in_workflow = False
        self.transformation=transformation
        self._inputs = []
        self._outputs = []
        self._dax_node = dax.Job(transformation)
        # NOTE: We are enforcing one site per transformation. Therefore the
        #       transformation used indicates the site to be used.
        self.set_execution_site(list(transformation.sites.keys())[0])
        self._args = []
        # Each value in _options is added separated with whitespace
        # so ['--option','value'] --> "--option value"
        self._options = []
        # For _raw_options *NO* whitespace is added.
        # so ['--option','value'] --> "--optionvalue"
        # and ['--option',' ','value'] --> "--option value"
        self._raw_options = []

    def add_arg(self, arg):
        """ Add an argument
        """
        if not isinstance(arg, File):
            arg = str(arg)

        self._args += [arg]

    def add_raw_arg(self, arg):
        """ Add an argument to the command line of this job, but do *NOT* add
            white space between arguments. This can be added manually by adding
            ' ' if needed
        """
        if not isinstance(arg, File):
            arg = str(arg)

        self._raw_options += [arg]

    def add_opt(self, opt, value=None):
        """ Add a option
        """
        if value is not None:
            if not isinstance(value, File):
                value = str(value)
            self._options += [opt, value]
        else:
            self._options += [opt]

    #private functions to add input and output data sources/sinks
    def _add_input(self, inp):
        """ Add as source of input data
        """
        self._inputs += [inp]
        self._dax_node.add_inputs(inp)

    def _add_output(self, out):
        """ Add as destination of output data
        """
        self._outputs += [out]
        out.node = self
        stage_out = out.storage_path is not None
        self._dax_node.add_outputs(out, stage_out=stage_out)

    # public functions to add options, arguments with or without data sources
    def add_input(self, inp):
        """Declares an input file without adding it as a command-line option.
        """
        self._add_input(inp)

    def add_output(self, inp):
        """Declares an output file without adding it as a command-line option.
        """
        self._add_output(inp)

    def add_input_opt(self, opt, inp):
        """ Add an option that determines an input
        """
        self.add_opt(opt, inp._dax_repr())
        self._add_input(inp)

    def add_output_opt(self, opt, out):
        """ Add an option that determines an output
        """
        self.add_opt(opt, out._dax_repr())
        self._add_output(out)

    def add_output_list_opt(self, opt, outputs):
        """ Add an option that determines a list of outputs
        """
        self.add_opt(opt)
        for out in outputs:
            self.add_opt(out)
            self._add_output(out)

    def add_input_list_opt(self, opt, inputs):
        """ Add an option that determines a list of inputs
        """
        self.add_opt(opt)
        for inp in inputs:
            self.add_opt(inp)
            self._add_input(inp)

    def add_list_opt(self, opt, values):
        """ Add an option with a list of non-file parameters.
        """
        self.add_opt(opt)
        for val in values:
            self.add_opt(val)

    def add_input_arg(self, inp):
        """ Add an input as an argument
        """
        self.add_arg(inp._dax_repr())
        self._add_input(inp)

    def add_output_arg(self, out):
        """ Add an output as an argument
        """
        self.add_arg(out._dax_repr())
        self._add_output(out)

    def new_output_file_opt(self, opt, name):
        """ Add an option and return a new file handle
        """
        fil = File(name)
        self.add_output_opt(opt, fil)
        return fil

    # functions to describe properties of this node
    def add_profile(self, namespace, key, value):
        """ Add profile information to this node at the DAX level
        """
        self._dax_node.add_profiles(
            dax.Namespace(namespace),
            key=key,
            value=value
        )

    def _finalize(self):
        if len(self._raw_options):
            raw_args = [''.join([str(a) for a in self._raw_options])]
        else:
            raw_args = []
        args = self._args + raw_args + self._options
        self._dax_node.add_args(*args)


class Workflow(object):
    """
    """
    def __init__(self, name='my_workflow', directory=None, cache_file=None,
                 dax_file_name=None):
        self.name = name
        self._rc = dax.ReplicaCatalog()
        self._tc = dax.TransformationCatalog()

        if directory is None:
            self.out_dir = os.getcwd()
        else:
            self.out_dir = os.path.abspath(directory)

        if cache_file is not None:
            cache_file = os.path.abspath(cache_file)
        self.cache_file = cache_file

        self._inputs = []
        self._outputs = []
        self._transformations = []
        self._containers = []
        self.in_workflow = False
        self.sub_workflows = []
        if dax_file_name is None:
            self.filename = self.name + '.dax'
        else:
            self.filename = dax_file_name
        self._adag = dax.Workflow(self.filename)

        # A pegasus job version of this workflow for use if it isncluded
        # within a larger workflow
        self._as_job = SubWorkflow(self.filename, is_planned=False,
                                   _id=self.name)
        self._swinputs = []

    def add_workflow(self, workflow):
        """ Add a sub-workflow to this workflow

        This function adds a sub-workflow of Workflow class to this workflow.
        Parent child relationships are determined by data dependencies

        Parameters
        ----------
        workflow : Workflow instance
            The sub-workflow to add to this one
        """
        workflow.in_workflow = self
        self.sub_workflows += [workflow]
        self._adag.add_jobs(workflow._as_job)
        return self

    def add_explicit_dependancy(self, parent, child):
        """
        Add an explicit dependancy between two Nodes in this workflow.

        Most dependencies (in PyCBC and Pegasus thinking) are added by
        declaring file linkages. However, there are some cases where you might
        want to override that and add an explicit dependancy.

        Parameters
        ----------
        parent : Node instance
            The parent Node.
        child : Node instance
            The child Node
        """
        self._adag.add_dependency(parent._dax_node, children=[child._dax_node])

    def add_subworkflow_dependancy(self, parent_workflow, child_workflow):
        """
        Add a dependency between two sub-workflows in this workflow

        This is done if those subworkflows are themselves declared as Workflows
        which are sub-workflows and not explicit SubWorkflows. (These Workflows
        contain SubWorkflows inside them .... Yes, the relationship between
        PyCBC and Pegasus becomes confusing here). If you are working with
        explicit SubWorkflows these can be added normally using File relations.

        Parameters
        ----------
        parent_workflow : Workflow instance
            The sub-workflow to use as the parent dependence.
            Must be a sub-workflow of this workflow.
        child_workflow : Workflow instance
            The sub-workflow to add as the child dependence.
            Must be a sub-workflow of this workflow.
        """
        self._adag.add_dependency(parent_workflow._as_job,
                                  children=[child_workflow._as_job])

    def add_transformation(self, tranformation):
        """ Add a transformation to this workflow

        Adds the input transformation to this workflow.

        Parameters
        ----------
        transformation : Pegasus.api.Transformation
            The transformation to be added.
        """
        self._tc.add_transformations(tranformation)

    def add_container(self, container):
        """ Add a container to this workflow

        Adds the input container to this workflow.

        Parameters
        ----------
        container : Pegasus.api.Container
            The container to be added.
        """
        self._tc.add_containers(container)

    def add_node(self, node):
        """ Add a node to this workflow

        This function adds nodes to the workflow. It also determines
        parent/child relations from the inputs to this job.

        Parameters
        ----------
        node : pycbc.workflow.pegasus_workflow.Node
            A node that should be executed as part of this workflow.
        """
        node._finalize()
        node.in_workflow = self

        # Record the executable that this node uses
        if node.transformation not in self._transformations:
            for tform in self._transformations:
                # Check if transform is already in workflow
                if node.transformation.is_same_as(tform):
                    node.transformation.in_workflow = True
                    node._dax_node.transformation = tform.name
                    node.transformation.name = tform.name
                    break
            else:
                self._transformations += [node.transformation]
                lgc = (hasattr(node, 'executable')
                       and node.executable.container is not None
                       and node.executable.container not in self._containers)
                if lgc:
                    self._containers.append(node.executable.container)

        # Add the node itself
        self._adag.add_jobs(node._dax_node)

        # Determine the parent child relationships based on the inputs that
        # this node requires.
        # In Pegasus5 this is mostly handled by pegasus, we just need to
        # connect files correctly if dealing with file management between
        # workflows/subworkflows
        for inp in node._inputs:
            if inp.node is not None and inp.node.in_workflow == self:
                # Standard case: File produced within the same workflow.
                # Don't need to do anything here.
                continue

            elif inp.node is not None and not inp.node.in_workflow:
                # This error should be rare, but can happen. If a Node hasn't
                # yet been added to a workflow, this logic breaks. Always add
                # nodes in order that files will be produced.
                raise ValueError('Parents of this node must be added to the '
                                 'workflow first.')

            elif inp.node is None:
                # File is external to the workflow (e.g. a pregenerated
                # template bank). (if inp.node is None)
                if inp not in self._inputs:
                    self._inputs += [inp]

            elif inp.node.in_workflow != self:
                # File is coming from a parent workflow, or other workflow
                # These needs a few extra hooks later, use _swinputs for this.
                if inp not in self._inputs:
                    self._inputs += [inp]
                    self._swinputs += [inp]
            else:
                err_msg = ("I don't understand how to deal with an input file "
                           "here. Ian doesn't think this message should be "
                           "possible, but if you get here something has gone "
                           "wrong and will need debugging!")
                raise ValueError(err_msg)

        # Record the outputs that this node generates
        self._outputs += node._outputs

        return self

    def __add__(self, other):
        if isinstance(other, Node):
            return self.add_node(other)
        elif isinstance(other, Workflow):
            return self.add_workflow(other)
        else:
            raise TypeError('Cannot add type %s to this workflow' % type(other))

    def traverse_workflow_io(self):
        """ If input is needed from another workflow within a larger
        hierarchical workflow, determine the path for the file to reach
        the destination and add the file to workflows input / output as
        needed.
        """
        def root_path(v):
            path = [v]
            while v.in_workflow:
                path += [v.in_workflow]
                v = v.in_workflow
            return path

        for inp in self._swinputs:
            workflow_root = root_path(self)
            input_root = root_path(inp.node.in_workflow)
            for step in workflow_root:
                if step in input_root:
                    common = step
                    break

            # Set our needed file as output so that it gets staged upwards
            # to a workflow that contains the job which needs it.
            for idx in range(input_root.index(common)):
                child_wflow = input_root[idx]
                parent_wflow = input_root[idx+1]
                if inp not in child_wflow._as_job.get_outputs():
                    child_wflow._as_job.add_outputs(inp, stage_out=True)
                    parent_wflow._outputs += [inp]

            # Set out needed file so it gets staged downwards towards the
            # job that needs it.
            for wf in workflow_root[:workflow_root.index(common)]:
                if inp not in wf._as_job.get_inputs():
                    wf._as_job.add_inputs(inp)

        for wf in self.sub_workflows:
            wf.traverse_workflow_io()

    def save(self, filename=None, submit_now=False, plan_now=False,
             output_map_path=None, root=True):
        """ Write this workflow to DAX file
        """
        if filename is None:
            filename = self.filename

        if output_map_path is None:
            output_map_path = 'output.map'

        # Handle setting up io for inter-workflow file use ahead of time
        # so that when daxes are saved the metadata is complete
        if root:
            self.traverse_workflow_io()

        for sub in self.sub_workflows:
            sub.save(root=False)
            # FIXME: If I'm now putting output_map here, all output_map stuff
            #        should move here.
            sub.output_map_file.insert_into_dax(self._rc, self.sites)
            sub_workflow_file = File(sub.filename)
            pfn = os.path.join(os.getcwd(), sub.filename)
            sub_workflow_file.add_pfn(pfn, site='local')
            sub_workflow_file.insert_into_dax(self._rc, self.sites)

        # add workflow input files pfns for local site to dax
        for fil in self._inputs:
            fil.insert_into_dax(self._rc, self.sites)

        self._adag.add_replica_catalog(self._rc)

        # Add TC into workflow
        self._adag.add_transformation_catalog(self._tc)

        with open(output_map_path, 'w') as f:
            for out in self._outputs:
                try:
                    f.write(out.output_map_str() + '\n')
                except ValueError:
                    # There was no storage path
                    pass

        # Pegasus requires that we write the DAX file into the local directory
        olddir = os.getcwd()
        os.chdir(self.out_dir)
        self._adag.write(filename)
        if not self.in_workflow:
            if submit_now or plan_now:
                self.plan_and_submit(submit_now=submit_now)
            else:
                with open('additional_planner_args.dat', 'w') as f:
                    stage_site_str = self.staging_site_str
                    exec_sites = self.exec_sites_str
                    # For now we don't include --config as this can be added to
                    # in submit_dax. We should add an option to add additional
                    # pegasus properties (through the config files?) here.
                    #prop_file = os.path.join(PEGASUS_FILE_DIRECTORY,
                    #                         'pegasus-properties.conf')
                    #f.write('--conf {} '.format(prop_file))
                    if self.cache_file is not None:
                        f.write('--cache {} '.format(self.cache_file))

                    f.write('--output-sites local ')
                    f.write('--sites {} '.format(exec_sites))
                    f.write('--staging-site {} '.format(stage_site_str))
                    f.write('--cluster label,horizontal ')
                    f.write('--cleanup inplace ')
                    f.write('--relative-dir work ')
                    # --dir is not being set here because it might be easier to
                    # set this in submit_dax still?
                    f.write('-q ')
                    f.write('--dax {}'.format(filename))
        os.chdir(olddir)

    def plan_and_submit(self, submit_now=True):
        """ Plan and submit the workflow now.
        """
        # New functionality, this might still need some work. Here's things
        # that this might want to do, that submit_dax does:
        # * Checks proxy (ignore this, user should already have this done)
        # * Pulls properties file in (DONE)
        # * Send necessary options to the planner (DONE)
        # * Some logging about hostnames (NOT DONE, needed?)
        # * Setup the helper scripts (start/debug/stop/status) .. (DONE)
        # * Copy some of the interesting files into workflow/ (DONE)
        # * Checks for dashboard URL (NOT DONE)
        # * Does something with condor_reschedule (NOT DONE, needed?)

        planner_args = {}
        planner_args['submit'] = submit_now

        # Get properties file - would be nice to add extra properties here.
        prop_file = os.path.join(PEGASUS_FILE_DIRECTORY,
                                 'pegasus-properties.conf')
        planner_args['conf'] = prop_file

        # Cache file, if there is one
        if self.cache_file is not None:
            planner_args['cache'] = [self.cache_file]

        # Not really sure what this does, but Karan said to use it. Seems to
        # matter for subworkflows
        planner_args['output_sites'] = ['local']

        # Site options
        planner_args['sites'] = self.sites
        planner_args['staging_sites'] = self.staging_site

        # Make tmpdir for submitfiles
        # default directory is the system default, but is overrideable
        # This should probably be moved to core.py?
        submit_opts = 'pegasus_profile', 'pycbc|submit-directory'
        submit_dir = None
        if self.cp.has_option(*submit_opts):
            submit_dir = self.cp.get(*submit_opts)
        submitdir = tempfile.mkdtemp(prefix='pycbc-tmp_', dir=submit_dir)
        os.chmod(submitdir, 0o755)
        try:
            os.remove('submitdir')
        except FileNotFoundError:
            pass
        os.symlink(submitdir, 'submitdir')
        planner_args['dir'] = submitdir

        # Other options
        planner_args['cluster'] = ['label,horizontal']
        planner_args['relative_dir'] = 'work'
        planner_args['cleanup'] = 'inplace'
        # This quietens the planner a bit. We cannot set the verbosity
        # directly, which would be better. So be careful, if changing the
        # pegasus.mode property, it will change the verbosity (a lot).
        planner_args['quiet'] = 1

        # FIXME: The location of output.map is hardcoded in the properties
        #        file. This is overridden for subworkflows, but is not for
        #        main workflows with submit_dax. If we ever remove submit_dax
        #        we should include the location explicitly here.
        self._adag.plan(**planner_args)

        # Set up convenience scripts
        with open('status', 'w') as fp:
            fp.write('pegasus-status --verbose ')
            fp.write('--long {}/work $@'.format(submitdir))

        with open('debug', 'w') as fp:
            fp.write('pegasus-analyzer -r ')
            fp.write('-v {}/work $@'.format(submitdir))

        with open('stop', 'w') as fp:
            fp.write('pegasus-remove {}/work $@'.format(submitdir))

        with open('start', 'w') as fp:
            if self.cp.has_option('pegasus_profile', 'pycbc|check_grid'):
                fp.write(GRID_START_TEMPLATE)
                fp.write('\n')
            fp.write('pegasus-run {}/work $@'.format(submitdir))

        os.chmod('status', 0o755)
        os.chmod('debug', 0o755)
        os.chmod('stop', 0o755)
        os.chmod('start', 0o755)

        os.makedirs('workflow/planning', exist_ok=True)

        shutil.copy2(prop_file, 'workflow/planning')
        shutil.copy2(os.path.join(submitdir, 'work', 'braindump.yml'),
                     'workflow/planning')

        if self.cache_file is not None:
            shutil.copy2(self.cache_file, 'workflow/planning')


class SubWorkflow(dax.SubWorkflow):
    """Workflow job representation of a SubWorkflow.

    This follows the Pegasus nomenclature where there are Workflows, Jobs and
    SubWorkflows. Be careful though! A SubWorkflow is actually a Job, not a
    Workflow. If creating a sub-workflow you would create a Workflow as normal
    and write out the necessary dax files. Then you would create a SubWorkflow
    object, which acts as the Job in the top-level workflow. Most of the
    special linkages that are needed for sub-workflows are then handled at that
    stage. We do add a little bit of functionality here.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.pycbc_planner_args = {}

    def add_into_workflow(self, container_wflow):
        """Add this Job into a container Workflow
        """
        self.add_planner_args(**self.pycbc_planner_args)

        # Set this to None so code will fail if more planner args are added
        self.pycbc_planner_args = None
        container_wflow._adag.add_jobs(self)

    def add_planner_arg(self, value, option):
        if self.pycbc_planner_args is None:
            err_msg = ("We cannot add arguments to the SubWorkflow planning "
                       "stage after this is added to the parent workflow.")
            raise ValueError(err_msg)

        self.pycbc_planner_args[value] = option

    def set_subworkflow_properties(self, output_map_file,
                                   staging_site,
                                   cache_file):

        self.add_planner_arg('pegasus.dir.storage.mapper.replica.file',
                             os.path.basename(output_map_file.name))
        self.add_inputs(output_map_file)

        # I think this is needed to deal with cases where the subworkflow file
        # does not exist at submission time.
        bname = os.path.splitext(os.path.basename(self.file))[0]
        self.add_planner_arg('basename',  bname)
        self.add_planner_arg('output_sites', ['local'])
        self.add_planner_arg('cleanup', 'inplace')
        self.add_planner_arg('cluster', ['label', 'horizontal'])
        self.add_planner_arg('verbose', 3)

        if cache_file:
            self.add_planner_arg('cache', [cache_file])

        if staging_site:
            self.add_planner_arg('staging_sites', staging_site)


class File(dax.File):
    """ The workflow representation of a physical file

    An object that represents a file from the perspective of setting up a
    workflow. The file may or may not exist at the time of workflow generation.
    If it does, this is represented by containing a physical file name (PFN).
    A storage path is also available to indicate the desired final
    destination of this file.
    """
    def __init__(self, name):
        self.name = name
        self.node = None
        dax.File.__init__(self, name)
        # Storage_path is where the file would be *output* to
        self.storage_path = None
        # Input_pfns is *input* locations of the file. This needs a site.
        self.input_pfns = []
        # Adding to a dax finalizes the File. Ensure that changes cannot be
        # made after doing this.
        self.added_to_dax = False

    def _dax_repr(self):
        return self

    @property
    def dax_repr(self):
        """Return the dax representation of a File."""
        return self._dax_repr()

    def output_map_str(self):
        if self.storage_path:
            return '%s %s pool="%s"' % (self.name, self.storage_path, 'local')
        else:
            raise ValueError('This file does not have a storage path')

    def add_pfn(self, url, site):
        """
        Associate a PFN with this file. Takes a URL and associated site.
        """
        self.input_pfns.append((url, site))

    def has_pfn(self, url, site='local'):
        """
        Check if the url, site is already associated to this File. If site is
        not provided, we will assume it is 'local'.
        """
        return (((url, site) in self.input_pfns)
                or ((url, 'all') in self.input_pfns))

    def insert_into_dax(self, rep_cat, sites):
        for (url, site) in self.input_pfns:
            if site == 'all':
                for curr_site in sites:
                    rep_cat.add_replica(curr_site, self, url)
            else:
                rep_cat.add_replica(site, self, url)

    @classmethod
    def from_path(cls, path):
        """Takes a path and returns a File object with the path as the PFN."""
        logging.warn("The from_path method in pegasus_workflow is deprecated. "
                     "Please use File.from_path (for output files) in core.py "
                     "or resolve_url_to_file in core.py (for input files) "
                     "instead.")
        urlparts = urlsplit(path)
        site = 'nonlocal'
        if (urlparts.scheme == '' or urlparts.scheme == 'file'):
            if os.path.isfile(urlparts.path):
                path = os.path.abspath(urlparts.path)
                path = urljoin('file:', pathname2url(path))
                site = 'local'

        fil = cls(os.path.basename(path))
        fil.add_pfn(path, site=site)
        return fil
