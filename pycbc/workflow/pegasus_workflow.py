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
import Pegasus.DAX3 as dax
import os

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

class Executable(ProfileShortcuts):
    """ The workflow representation of an Executable 
    """
    id = 0
    def __init__(self, name, namespace=None, os='linux', 
                       arch='x86_64', installed=True, version=None):
        self.name = name
        self.logical_name = self.name + "_ID%s" % str(Executable.id)
        Executable.id += 1
        self.namespace = namespace
        self.version = version
        self._dax_executable = dax.Executable(self.logical_name, 
                   namespace=self.namespace, version=version, os=os, 
                   arch=arch, installed=installed) 
        self.in_workflow = False
        self.pfns = {}
        
    def add_pfn(self, url, site='local'):
        self._dax_executable.PFN(url, site)
        self.pfns[site] = url
        
    def get_pfn(self, site='local'):
        return self.pfns[site]
        
    def insert_into_dax(self, dax):
        dax.addExecutable(self._dax_executable)
        
    def add_profile(self, namespace, key, value):
        """ Add profile information to this executable
        """
        try:
            entry = dax.Profile(namespace, key, value)
            self._dax_executable.addProfile(entry)  
        except dax.DuplicateError:
            pass
 
class Node(ProfileShortcuts):    
    def __init__(self, executable):
        self.in_workflow = False   
        self.executable=executable            
        self._inputs = []
        self._outputs = []        
        self._dax_node = dax.Job(name=executable.logical_name,
                                 version = executable.version, 
                                 namespace=executable.namespace)
        self._args = []
        self._options = []
        
    def add_arg(self, arg):
        """ Add an argument
        """
        if not isinstance(arg, File):
            arg = str(arg)
        
        self._args += [arg]
        
    def add_opt(self, opt, value=None):
        """ Add a option
        """
        if value:
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
        inp._set_as_input_of(self)
      
    def _add_output(self, out):
        """ Add as destination of output data
        """
        self._outputs += [out]
        out.node = self
        out._set_as_output_of(self)
        
    # public functions to add options, arguments with or without data sources
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
        try:
            entry = dax.Profile(namespace, key, value)
            self._dax_node.addProfile(entry)
        except dax.DuplicateError:
            pass    
        
    def _finalize(self):
        args = self._args + self._options
        self._dax_node.addArguments(*args)
        
class Workflow(object):
    """ 
    """
    def __init__(self, name='my_workflow'):
        self.name = name
        self._adag = dax.ADAG(name)
        
        self._inputs = []
        self._outputs = []
        self._executables = []
        self.in_workflow = False
        self.sub_workflows = []
        self._external_workflow_inputs = []
        self.filename = self.name + '.dax'
        
        self.as_job = dax.DAX(self.filename)
        
    def _make_root_dependency(self, inp):  
        def root_path(v):
            path = []
            while v.in_workflow:
                path += [v.in_workflow]
                v = v.in_workflow
            return path
                
        workflow_root = root_path(self)
        input_root = root_path(inp)
        
        for step in workflow_root:
            if step in input_root:
                common = step
                break
        
        dep = dax.Dependency(child=workflow_root[workflow_root.index(common)-1],
                             parent=input_root[input_root.index(common)-1])
        common._adag.addDependency(dep)
      
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
        
        node = workflow.as_job
        self._adag.addJob(node)
        
        node.file.PFN(os.path.join(os.getcwd(), node.file.name), site='local')
        self._adag.addFile(node.file)
        
        for inp in self._external_workflow_inputs:
            workflow._make_root_dependency(inp.node)
            
        return self
        
        
    def add_node(self, node):
        """ Add a node to this workflow
        
        This function adds nodes to the workflow. It also determines
        parent/child relations from the DataStorage inputs to this job. 
        
        Parameters
        ----------
        node : Node
            A node that should be exectuded as part of this workflow.
        """
        node._finalize()
        node.in_workflow = self
        self._adag.addJob(node._dax_node)
 
        # Determine the parent child relationships based on the inputs that
        # this node requires.
        added_nodes = []
        for inp in node._inputs:
            if inp.node is not None and inp.node.in_workflow == self:
                if inp.node not in added_nodes:
                    parent = inp.node._dax_node
                    child = node._dax_node
                    dep = dax.Dependency(parent=parent, child=child)
                    self._adag.addDependency(dep)    
                    added_nodes.append(inp.node)                  
                            
            elif inp.node is not None and not inp.node.in_workflow:
                raise ValueError('Parents of this node must be added to the '
                                 'workflow first.')   
                      
            elif inp.node is None and not inp.workflow_input:
                self._inputs += [inp]
                inp.workflow_input = True    
            
            elif inp.node is not None and inp.node.in_workflow != self and inp not in self._inputs:
                self._inputs += [inp]
                self._external_workflow_inputs += [inp]    

            
        # Record the outputs that this node generates                                   
        self._outputs += node._outputs  
            
        # Record the executable that this node uses
        if not node.executable.in_workflow:
            node.executable.in_workflow = True
            self._executables += [node.executable]
        
        return self
        
    def __add__(self, other):
        if isinstance(other, Node):
            return self.add_node(other)
        elif isinstance(other, Workflow):
            return self.add_workflow(other)
        else:
            raise TypeError('Cannot add type %s to this workflow' % type(other))
            
      
    def save(self):
        """ Write this workflow to DAX file 
        """        
        for sub in self.sub_workflows:
            sub.save()
        
        f = open(self.filename, "w")
        self._adag.writeXML(f)

class DataStorage(object):
    """ A workflow representation of a place to store and read data from. 
    
    The abstract representation of a place to store and read data from. This 
    can include files, database, or remote connections. This object is
    used as a handle to pass between functions, and is used a way to logically
    represent the order operation on the physical data. 
    """
    def __init__(self, name):
        self.name = name      
        self.node = None
        self.workflow_input = False
        
    def _set_as_node_input(self):
        pass
        
    def _set_as_node_output(self):
        pass
        
    def _dax_repr(self):
        return self.name
        
class File(DataStorage, dax.File):
    """ The workflow representation of a physical file
    
    An object that represents a file from the perspective of setting up a 
    workflow. The file may or may not exist at the time of workflow generation.
    If it does, this is represented by containing a physical file name (PFN).
    A storage path is also available to indicate the desired final 
    destination of this file.
    """
    def __init__(self, name):
        DataStorage.__init__(self, name)
        dax.File.__init__(self, name)
        self.storage_path = None

    def _dax_repr(self):
        return self
        
    def _set_as_input_of(self, node):
        node._dax_node.uses(self, link=dax.Link.INPUT, register=False, 
                                                       transfer=True)          
    def _set_as_output_of(self, node):
        if self.storage_path:
            transfer_file = True
        else:
            transfer_file = False
        node._dax_node.uses(self, link=dax.Link.OUTPUT, register=False, 
                                                        transfer=transfer_file)                                                       
    def output_map_str(self):
        if self.storage_path:
            return '%s %s pool="%s"' % (self.name, self.storage_path, 'local')
        else:
            raise ValueError('This file does not have a storage path')
        
    def insert_into_dax(self, dax):
        dax.addFile(self)
    
class Database(DataStorage):
    pass
