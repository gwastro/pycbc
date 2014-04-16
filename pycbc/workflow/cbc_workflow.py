from pycbc.workflow.workflow import Workflow, Node, File
from pycbc.ahope.configparserutils import AhopeConfigParser
import re
import Pegasus.DAX3 as dax

class AhopeWorkflow(Workflow):
    """
    This class manages an aHOPE style workflow. It provides convenience 
    functions for finding input files using time and keywords. It can also
    generate cache files from the inputs. It makes heavy use of the
    pipeline.CondorDAG class, which is instantiated under self.dag.
    """
    def __init__(self, args, name):
        """
        Create an aHOPE workflow
        
        Parameters
        ----------
        args : argparse.ArgumentParser
            The command line options to initialize an ahope workflow.
        """
        Workflow.__init__(self, name)
        
        # Parse ini file
        self.cp = AhopeConfigParser.from_args(args)
        
        # Dump the parsed config file
        ini_file = os.path.abspath(self.name + '_parsed.ini')
        if not os.path.isfile(iniFile):
            fp = open(iniFile, 'w')
            self.cp.write(fp)
            fp.close()
        else:
            logging.warn("Cowardly refusing to overwrite %s." %(ini_file))

        # Set global values
        start_time = int(self.cp.get("ahope", "start-time"))
        end_time = int(self.cp.get("ahope", "end-time"))
        self.analysis_time = segments.segment([start_time, end_time])

        # Set the ifos to analyse
        ifos = []
        for ifo in self.cp.options('ahope-ifos'):
            ifos.append(ifo.upper())
        self.ifos = ifos
        self.ifos.sort(key=str.lower)
        self.ifo_string = ''.join(self.ifos)
        
        # Set up input and output file lists for workflow
        self._inputs = AhopeFileList([])
        self._outputs = AhopeFileList([])
         
    def execute_node(self, node):
        """ Execute this node immediately on the local site and place its
        inputs and outputs into the workflow data lists. 
        """
        self.node.executed = True
        cmd_list = node.get_command_line()
        
        # Must execute in output directory.
        currDir = os.getcwd()
        os.chdir(job_dir)
        
        # Make call
        make_external_call(cmd_list, out_dir=os.path.join(job_dir, 'logs'),
                                     out_baseame=base_name) 
        # Change back
        os.chdir(currDir)
        
        self._outputs += node._outputs
            
    def save(self, basename):
        # add executable pfns for local site to dax
        
        # add workflow input files pfns for local site to dax
    
        # save the dax file
        Workflow.save(self, basename + '.dax')
    
class AhopeNode(Node):
    def __init__(self, name):
        Node.__init__(self, name)
        self.executed = False
    
    def get_command_line(self):
        arglist = self._dax_node.arguments
        arglist = [a.path if isinstance(a, dax.File) else a for a in arglist]
        print arglist
                        
        exe_path = urlparse.urlsplit(self.executable.get_pfn()).path
        return [exe_path] + arglist
        
    def new_output_file_opt(self, opt):
        
        fil = Node.new_output_file_opt(self, opt, name)
    
class AhopeFile(File):
    '''
    This class holds the details of an individual output file 
    This file(s) may be pre-supplied, generated from within the ahope
    command line script, or generated within the workflow. The important stuff
    is:

    * The ifo that the AhopeFile is valid for
    * The time span that the AhopeOutFile is valid for
    * A short description of what the file is
    * The extension that the file should have
    * The url where the file should be located

    An example of initiating this class:
    
    c = AhopeFile("H1", "INSPIRAL_S6LOWMASS", segments.segment(815901601, 815902001), file_url="file://localhost/home/spxiwh/H1-INSPIRAL_S6LOWMASS-815901601-400.xml.gz" )

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
            self.segment_list = segments.segmentlist([segs])
        elif isinstance(segs, (segments.segmentlist)):
            self.segment_list = segs
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
       
        cache_entry = lal.CacheEntry(self.ifoString,
                   self.tagged_description, self.segList.extent(), file_url)
        # Make a cyclical reference
        cache_entry.ahope_file = self
        self.cache_entry = cache_entry
        # This gets set if this becomes an input file in the workflow
        self.is_workflow_input = False
    
    @property
    def url(self):
        return self.cache_entry.url
       
    @property
    def path(self):
        """
        If only one file is contained in this instance this will be that path.
        Otherwise a TypeError is raised.
        """
        return self.cache_entry.path

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
        return basename(self.cache_entry.path)
        
    def _filename(self, ifo, description, extension, segment):
        """
        Construct the standard output filename. Should only be used internally
        of the AhopeFile class.
        """        
        if extension.startswith('.'):
            extension = extension[1:]
            
        # Follow the frame convention of using integer filenames,
        # but stretching to cover partially covered seconds.
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
