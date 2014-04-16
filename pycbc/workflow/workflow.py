import Pegasus.DAX3
import copy

class Executable(object):
    def __init__(self, name):
        self.name = name

class Node(object):
    _id = 0
    
    def __init__(self, executable):
        self.in_workflow = False
    
        self.executable=executable 
            
        self._inputs = []
        self._outputs = []
        
        self._dax_node = Pegasus.DAX3.Job(name=executable.name, 
                                          id=self._get_node_id())
    
    def _get_node_id(self):
        Node._id += 1
        return "ID%06d" % self._id
        
    def add_arg(self, arg):
        """ Add an argument
        """
        self._dax_node.addArguments(arg)
        
        
    def add_opt(self, opt, value=None):
        """ Add a option
        """
        if value:
            self._dax_node.addArguments(opt, value)
        else:
            
        
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

    def new_output_file_arg(self, name):
        fil = File(name)
        self.add_output_arg(fil)
        return fil

class Workflow(object):
    def __init__(self, name='my_workflow'):
        self.name = name
        self._adag = Pegasus.DAX3.ADAG(name)
        
        self._inputs = []
        self._outputs = []
        
    def add_node(self, node):
        node.in_workflow = True
        self._adag.addJob(node._dax_node)
    
        for inp in node._inputs:
            if inp.node is not None and inp.node.in_workflow:
                parent = inp.node._dax_node
                child = node._dax_node
                dep = Pegasus.DAX3.Dependency(parent=parent, child=child)
                self._adag.addDependency(dep)    
                            
            elif inp.node is not None and not inp.node.in_workflow:
                raise ValueError('Parents of this node must be added to the '
                                 'workflow first.')   
                                  
            elif inp.node is None and inp.workflow_input is False:
                self._inputs += [inp]
                inp.workflow_input = True
                                               
        for out in node._outputs:
            self._outputs += node._outputs  
        
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
        self.workfow_input = False
        
    def _set_as_node_input(self):
        pass
        
    def _set_as_node_output(self):
        pass
        
    def _dax_repr(self):
        return self.name
        
class File(DataStorage):
    def __init__(self, name):
        DataStorage.__init__(self, name)


    def _dax_repr(self):
        return Pegasus.DAX3.File(self.name)

    def _set_as_input_of(self, node):
        fil = Pegasus.DAX3.File(self.name)
        node._dax_node.uses(fil, link=Pegasus.DAX3.Link.INPUT,
                                               register=False,
                                               transfer=True)
    
    def _set_as_output_of(self, node):
        fil = Pegasus.DAX3.File(self.name)
        node._dax_node.uses(fil, link=Pegasus.DAX3.Link.OUTPUT,
                                               register=False,
                                               transfer=True)
    
class Database(DataStorage):
    pass
