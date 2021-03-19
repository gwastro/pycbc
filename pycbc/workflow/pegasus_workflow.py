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
from six.moves.urllib.request import pathname2url
from six.moves.urllib.parse import urljoin, urlsplit
import Pegasus.api as dax
from .pegasus_sites import add_site

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
        if universe is 'standard':
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
                       arch='x86_64', installed=True,
                       container=None):
        self.logical_name = name + "_ID%s" % str(Executable.id)
        Executable.id += 1
        self.os = dax.OS(os)
        self.arch = dax.Arch(arch)
        self.installed = installed
        self.container = container
        self.in_workflow = False
        self.profiles = {}
        self.transformations = {}

    def get_transformation(self, site, url=None):
        if site in self.transformations:
            return self.transformations[site]
        else:
            return self.create_transformation(site, url=url)

    def create_transformation(self, site, url=None):
        if url is None:
            # Rely on URL being set in a child class
            # FIXME: Error handling niceness needed
            url = self.exe_pfns[site]

        # FIXME: Understand is_stageable. Do *not* want exes transferred in
        #        almost all cases. (Probably that's a site property)
        transform = dax.Transformation(
            self.logical_name,
            site=site,
            pfn=url,
            is_stageable=self.installed, # I think??
            arch=self.arch,
            os_type=self.os,
            container=self.container # Is it? There's some new container stuff
        )
        for (namespace,key), value in self.profiles.items():
            transform.add_profiles(
                dax.Namespace(namespace),
                key=key,
                value=value
            )
        self.transformations[site] = transform
        return transform

    def add_profile(self, namespace, key, value):
        """ Add profile information to this executable
        """
        if self.transformations:
            raise ValueError("Need code changes to be able to add profiles after transformations are created.")
        self.profiles[(namespace,key)] = value

    def is_same_as(self, other):
        test_vals = ['namespace', 'version', 'arch', 'os', 'osrelease',
                     'osversion', 'glibc', 'installed', 'container']
        # Check for logical name first
        if not self.pegasus_name == other.pegasus_name:
            return False

        # Check the properties of the executable
        for val in test_vals:
            sattr = getattr(self._dax_executable, val)
            oattr = getattr(other._dax_executable, val)
            if not sattr == oattr:
                return False
        # Also check the "profile". This is things like Universe, RAM/disk/CPU
        # requests, execution site, getenv=True, etc.
        for profile in self._dax_executable.profiles:
            if profile not in other._dax_executable.profiles:
                return False
        for profile in other._dax_executable.profiles:
            if profile not in self._dax_executable.profiles:
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

    def add_input(self, inp):
        """Declares an input file without adding it as a command-line option.
        """
        self._add_input(inp)

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
    def add_profile(self, namespace, key, value, force=False):
        """ Add profile information to this node at the DAX level
        """
        try:
            self._dax_node.add_profiles(
                dax.Namespace(namespace),
                key=key,
                value=value
            )
        except dax.errors.DuplicateError:
            if force:
                # FIXME: This definitely won't work. Not sure how to fix yet!
                raise
                # Replace with the new key
                self._dax_node.removeProfile(entry)
                self._dax_node.addProfile(entry)
            else:
                raise

    def _finalize(self):
        if len(self._raw_options):
            raw_args = [''.join([str(a) for a in self._raw_options])]
            print (raw_args)
        else:
            raw_args = []
        args = self._args + raw_args + self._options
        self._dax_node.add_args(*args)


class Workflow(object):
    """
    """
    def __init__(self, name='my_workflow', is_subworkflow=False):
        self.name = name
        self._rc = dax.ReplicaCatalog()
        self._tc = dax.TransformationCatalog()
        self._sc = dax.SiteCatalog()

        self._inputs = []
        self._outputs = []
        self._transformations = []
        self._containers = []
        self.in_workflow = False
        self.sub_workflows = []
        self._external_workflow_inputs = []
        self.filename = self.name + '.dax'
        self._adag = dax.Workflow(self.filename)
        if is_subworkflow:
            self._asdag = dax.SubWorkflow(self.filename, is_planned=False)
        else:
            self._asdag = None

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

        self._adag.add_jobs(workflow._asdag)

        for inp in workflow._external_workflow_inputs:
            self._adag.add_dependency(inp.node, children=[workflow._asdag])

        return self

    def add_subworkflow_dependancy(self, parent_workflow, child_workflow):
        """
        Add a dependency between two sub-workflows in this workflow

        Parameters
        ----------
        parent_workflow : Workflow instance
            The sub-workflow to use as the parent dependence.
            Must be a sub-workflow of this workflow.
        child_workflow : Workflow instance
            The sub-workflow to add as the child dependence.
            Must be a sub-workflow of this workflow.
        """
        self._adag.add_dependency(parent_workflow._asdag,
                                  children=[child_workflow._asdag])

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
        if not node.transformation in self._transformations:
            for exe in self._transformations:
                # Check if exe is already in self.executables properly
                # FIXME: Make this work again
                if 0: #node.executable.is_same_as(exe):
                    node.executable.in_workflow = True
                    node._dax_node.name = exe.logical_name
                    node.executable.logical_name = exe.logical_name
                    break
            else:
                #node.executable.in_workflow = True
                if node.transormation.site in self._tc.sites:
                    add_site(self._tc, node.transormation.site)
                self._transformations += [node.transformation]
                if hasattr(node, 'executable') and node.executable.container is not None and node.executable.container not in self._containers:
                    self._containers.append(node.executable.container)

        # Add the node itself
        self._adag.add_jobs(node._dax_node)

        # Determine the parent child relationships based on the inputs that
        # this node requires.
        added_nodes = []
        for inp in node._inputs:
            # FIXME: Make this do whatever it was supposed to do!
            # Breaking this loop for testing
            if inp.node is not None and inp.node.in_workflow == self:
                if inp.node not in added_nodes:
                    #parent = inp.node._dax_node
                    #child = node._dax_node
                    #dep = dax.Dependency(parent=parent, child=child)
                    #self._adag.addDependency(dep)
                    added_nodes.append(inp.node)

            elif inp.node is not None and not inp.node.in_workflow:
                raise ValueError('Parents of this node must be added to the '
                                 'workflow first.')

            elif inp.node is None and not inp.workflow_input:
                self._inputs += [inp]
                inp.workflow_input = True

            elif inp.node is not None and inp.node.in_workflow != self and inp not in self._inputs:
                # FIXME: Check if we still need this complication
                self._inputs += [inp]
                self._external_workflow_inputs += [inp]

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


    def save(self, filename=None, transformation_catalog_path=None):
        """ Write this workflow to DAX file
        """
        if filename is None: 
            filename = self.filename
        if transformation_catalog_path is None:
            transformation_catalog_path = self.transformation_catalog

        for sub in self.sub_workflows:
            sub.save()
            sub.transformation_catalog_file.insert_into_dax(self._rc)
            sub.output_map_file.insert_into_dax(self._rc)
            sub_workflow_file = File(sub.filename)
            pfn = os.path.join(os.getcwd(), sub.filename)
            sub_workflow_file.add_pfn(pfn, site='local')
            sub_workflow_file.insert_into_dax(self._rc)

        self._adag.add_replica_catalog(self._rc)
        # FIXME: Cannot add TC into workflow
        #self._adag.add_transformation_catalog(self._tc)

        self._adag.write(filename)
        self._tc.write(transformation_catalog_path)


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
        self.workflow_input = False
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
        self.input_pfns.append((url,site))

    def has_pfn(self, url, site='local'):
        """ 
        Check if the url, site is already associated to this File. If site is
        not provided, we will assume it is 'local'.
        """
        return (url,site) in self.input_pfns

    def insert_into_dax(self, rep_cat):
        for (url, site) in self.input_pfns: 
            rep_cat.add_replica(site, self, url)

    @classmethod
    def from_path(cls, path):
        """Takes a path and returns a File object with the path as the PFN."""
        urlparts = urlsplit(path)
        site = 'nonlocal'
        if (urlparts.scheme == '' or urlparts.scheme == 'file'):
            if os.path.isfile(urlparts.path):
                path = os.path.abspath(urlparts.path)
                path = urljoin('file:', pathname2url(path)) 
                site = 'local'

        fil = File(os.path.basename(path))
        fil.PFN(path, site)
        return fil
