import copy, urlparse
import Pegasus.DAX3 as dax

class Executable(object):
    id = 0

    def __init__(self, name):
        self.name = name
        self.logical_name = self.name + str(Executable.id)
        Executable.id += 1
        self.namespace = self.logical_name
        self._dax_executable = dax.Executable(self.logical_name, 
                                              namespace=self.namespace,
                                              version="1.0",
                                              os="linux", arch="x86_64",
                                              installed=True) 
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
        entry = dax.Profile(namespace, key, value)
        self._dax_executable.addProfile(entry)
        
    def set_memory(self, size):
        self.add_profile('condor', 'request_memory', '%sM' % size)
         
    def set_storage(self, size):
        self.add_profile('condor', 'request_disk', '%sM' % size)
        
    def set_num_cpus(self, number):
        self.add_profile('condor', 'request_cpus', number)
        
    def set_universe(self, universe):
        self.add_profile("condor", "universe", universe)
        
    def set_category(self, category):
        self.add_profile('dagman', 'category', category)
        
    def set_priority(self, priority):
        self.add_profile('dagman', 'priority', priority)
        
    def set_num_retries(self, number):
        self.add_profile("dagman", "retry", number)    
 
#class DaxNode(Pegasus.DAX3.DAX):
#    def __init__(self, filename):
#        pass
        
class Node(object):    
    def __init__(self, executable):
        self.in_workflow = False   
        self.executable=executable            
        self._inputs = []
        self._outputs = []        
        self._dax_node = dax.Job(name=executable.logical_name,
                                 version = "1.0",
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
        """ Add an option that contains an input 
        """
        self.add_opt(opt, inp._dax_repr())
        self._add_input(inp)
        
    def add_output_opt(self, opt, out):
        """ Add an option that contains an output
        """
        self.add_opt(opt, out._dax_repr())
        self._add_output(out)
        
    def add_output_list_opt(self, opt, outputs):
        self.add_opt(opt)
        for out in outputs:
            self.add_opt(out)
            self._add_output(out)
            
    def add_input_list_opt(self, opt, inputs):
        self.add_opt(opt)
        for inp in inputs:
            self.add_opt(inp)
            self._add_input(inp)
        
    def add_input_arg(self, inp):
        self.add_arg(inp._dax_repr())
        self._add_input(inp)
        
    def add_output_arg(self, out):
        self.add_arg(out._dax_repr())
        self._add_output(out)
        
    def new_output_file_opt(self, opt, name):
        fil = File(name)
        self.add_output_opt(opt, fil)
        return fil

    # functions to describe properties of this node
    def add_profile(self, namespace, key, value):
        """ Add profile information to this node at the DAX level
        """
        entry = dax.Profile(namespace, key, value)
        self._dax_node.addProfile(entry)
    
    def set_memory(self, size):
        self.add_profile('condor', 'request_memory', '%sM' % size)
         
    def set_storage(self, size):
        self.add_profile('condor', 'request_disk', '%sM' % size)
        
    def set_num_cpus(self, number):
        self.add_profile('condor', 'request_cpus', number)
        
    def set_universe(self, universe):
        self.add_profile("condor", "universe", universe)
        
    def set_category(self, category):
        self.add_profile('dagman', 'category', category)
        
    def set_priority(self, priority):
        self.add_profile('dagman', 'priority', priority)
        
    def set_num_retries(self, number):
        self.add_profile("dagman", "retry", number)
        
    def _finalize(self):
        args = self._args + self._options
        self._dax_node.addArguments(*args)
        
class Workflow(object):
    def __init__(self, name='my_workflow'):
        self.name = name
        self._adag = dax.ADAG(name)
        
        self._inputs = []
        self._outputs = []
        self._executables = []
        
    def add_node(self, node):
        node._finalize()
        node.in_workflow = True
        self._adag.addJob(node._dax_node)
 
        for inp in node._inputs:
            if inp.node is not None and inp.node.in_workflow:
                parent = inp.node._dax_node
                child = node._dax_node
                dep = dax.Dependency(parent=parent, child=child)
                self._adag.addDependency(dep)    
                            
            elif inp.node is not None and not inp.node.in_workflow:
                raise ValueError('Parents of this node must be added to the '
                                 'workflow first.')   
                                  
            elif inp.node is None and inp.workflow_input is False:
                self._inputs += [inp]
                inp.workflow_input = True
                                               
        for out in node._outputs:
            self._outputs += node._outputs  
            
        if not node.executable.in_workflow:
            node.executable.in_workflow = True
            self._executables += [node.executable]
        
        return self
        
    __iadd__ = add_node
      
    def save(self, filename):
        """ Write this workflow to DAX file 
        """
        f = open(filename, "w")
        self._adag.writeXML(f)

class DataStorage(object):
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
    def __init__(self, name):
        DataStorage.__init__(self, name)
        dax.File.__init__(self, name)
        
        self._dax_file = dax.File(self.name)
        self.storage_path = None

    def _dax_repr(self):
        return self
        
    def _set_as_input_of(self, node):
        node._dax_node.uses(self, link=dax.Link.INPUT, register=False, 
                                                       transfer=True)          
    def _set_as_output_of(self, node):
        node._dax_node.uses(self, link=dax.Link.OUTPUT, register=False, 
                                                        transfer=True)                                                       
    def output_map_str(self):
        return '%s %s pool="%s"' % (self.name, self.storage_path, 'local') 
        
    def insert_into_dax(self, dax):
        dax.addFile(self)
    
class Database(DataStorage):
    pass
