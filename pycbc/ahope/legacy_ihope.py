import os
import urlparse
from glue import pipeline
from glue import segments
from ahope_utils import Job, Node, Executable, AhopeFile

def legacy_get_valid_times(self, cp, ifo):
    """
    Return the length of data that the tmpltbank job will need to read and
    the part of that data that the template bank is valid for. In the case
    of lalapps_tmpltbank the following options are needed to set this up
    and will be used by the executable to figure this out:

    * --pad-data (seconds, amount of data used to pad the analysis region.
      This is needed as some data will be corrupted from the data
      conditioning process)

    * --segment-length (sample points, length of each analysis segment)

    * --sample-rate (Hz, number of sample points per second. The data will
      be resampled to this value if necessary

    * --number-of-segments (Number of analysis segments, note that
      overlapping segments are used for PSD estimation, so every data
      point will appear in two segments, except the first
      segment-length/4 and last segment-length/4 points.)

    Parameters
    ----------
    cp : ConfigParser object
        The ConfigParser object holding the ahope configuration settings
    ifo : string
        The interferometer being setup. It is possible to use different
        configuration settings for each ifo.

    Returns
    -------
    dataLength : float (seconds)
        The length of data that the job will need
    validChunk : glue.glue.segments.segment
        The start and end of the dataLength that is valid for the template
        bank.
    """
    # FIXME: This is only valid for templateBank not inspiral!
    # FIXME: Suggest making a separate inspiral function.

    # Read in needed options. This will fail if options not present
    # It will search relevant sub-sections for the option, so this can be
    # set differently for each ifo.
    padData = int(cp.get_opt_ifo(self.exe_name, 'pad-data', ifo))
    segmentLength = float(cp.get_opt_ifo(self.exe_name,
                                         'segment-length', ifo))
    sampleRate = float(cp.get_opt_ifo(self.exe_name,'sample-rate', ifo))
    numSegments = int(cp.get_opt_ifo(self.exe_name,
                                     'number-of-segments', ifo))
    # Calculate total valid duration
    analysisDur = int(segmentLength/sampleRate) * (numSegments + 1)/2
    if (segmentLength % sampleRate):
        errString = "In tmpltbank, when running lalapps_tmpltbank "
        errString += "segment-length must be a multiple of sample-rate."
        raise ValueError(errString)
    # Set the segments
    dataLength = analysisDur + 2*padData
    validStart = padData
    validEnd = analysisDur + padData
    # If this is inspiral we lose segment-length/4 on start and end
    if self.exe_name == 'inspiral':
        # Don't think inspiral will do well if segmentLength/4 is not
        # an integer
        validStart = validStart + int(segmentLength/(sampleRate * 4))
        validEnd = validEnd - int(segmentLength / (sampleRate * 4))
    validChunk = segments.segment([validStart,validEnd])

    return dataLength, validChunk

class LegacyAnalysisNode(Node, pipeline.AnalysisNode):
    set_jobnum_tag = pipeline.AnalysisNode.set_user_tag
    
        
class LegacyAnalysisJob(Job):
    def create_node(self, data_seg, valid_seg, parent=None, dfParents=None):
        node = LegacyAnalysisNode(self)
        
        if not dfParents or len(dfParents) != 1: 
            raise ValueError("%s must be supplied with a single cache file" 
                              %(self.exe_name))  
        
        pad_data = int(self.get_opt('pad-data'))
        if pad_data is None:
            raise ValueError("The option pad-data is a required option of "
                             "%s. Please check the ini file." % self.exe_name)                           
              
        node.set_start(data_seg[0] + pad_data)
        node.set_end(data_seg[1] - pad_data)
         
        cache_file = dfParents[0]       
        
        #check the extension       
        extension = '.xml'
        gzipped = self.get_opt('write-compress')
        if gzipped is not None:
            extension += '.gz'
        
        #create the ouptut file for this job
        name_segment = segments.segment([node.get_start(), node.get_end()])
        out_file = AhopeFile(self.ifo, self.exe_name, 
                             extension=extension,
                             segment=name_segment,
                             directory=self.out_dir)
        out_file.segment = valid_seg
        node.add_output(out_file)
        node.add_input(cache_file, opts='frame-cache')         
        return node
        
class LegacyInspiralJob(LegacyAnalysisJob):
    def create_node(self, data_seg, valid_seg, parent=None, dfParents=None):
        node = LegacyAnalysisJob.create_node(self, data_seg, valid_seg, 
                                                   parent, dfParents)
        node.set_trig_start(valid_seg[0])
        node.set_trig_end(valid_seg[1])  
        node.add_input(parent, opts='bank-file')    
        return node

class LegacyTmpltbankExec(Executable):
    def __init__(self, exe_name):
        if exe_name != 'tmpltbank':
            raise ValueError('lalapps_tmpltbank does not support setting '
                             'the exe_name to anything but "tmpltbank"')
                           
        Executable.__init__(self, 'tmpltbank', 'standard')

    def create_job(self, cp, ifo, out_dir=None):
        return LegacyAnalysisJob(cp, self.exe_name, self.condor_universe,
                                 ifo=ifo, out_dir=out_dir)
    
    get_valid_times = legacy_get_valid_times
        
class LegacyInspiralExec(Executable):
    def __init__(self, exe_name):
        if exe_name != 'inspiral':
            raise ValueError('lalapps_tmpltbank does not support setting '
                             'the exe_name to anything but "inspiral"')
        Executable.__init__(self, 'inspiral', 'standard')

    def create_job(self, cp, ifo, out_dir=None):
        return LegacyInspiralJob(cp, self.exe_name, self.condor_universe, 
                                 ifo=ifo, 
                                 out_dir=out_dir)

    get_valid_times = legacy_get_valid_times

class LegacySplitBankExec(Executable):
    """This class holds the function for lalapps_splitbank 
    usage following the old ihope specifications.
    """
    def __init__(self, exe_name):
        Executable.__init__(self, exe_name, 'standard')
    
    def create_job(self, cp, ifo, out_dir=None):
        return LegacySplitBankJob(cp, self.exe_name, self.condor_universe, 
                                  ifo=ifo,
                                  out_dir=out_dir)

class PyCBCInspiralJob(Job):
    def create_node(self, data_seg, valid_seg, parent=None, dfParents=None):
        node = LegacyAnalysisNode(self)
        
        pad_data = int(self.get_opt('pad-data'))
        if pad_data is None:
            raise ValueError("The option pad-data is a required option of "
                             "%s. Please check the ini file." % self.exe_name)
                             
        if not dfParents or len(dfParents) != 1: 
            raise ValueError("%s must be supplied with a single cache file" 
                              %(self.exe_name))     
        
        # set remaining options flags   
        node.set_start(data_seg[0] + pad_data)
        node.set_end(data_seg[1] - pad_data)       
        node.set_trig_start(valid_seg[0])
        node.set_trig_end(valid_seg[1])   
        
        cache_file = dfParents[0]                 
                            
        # FIXME add control for output type         
        insp = AhopeFile(self.ifo, self.exe_name, extension='.xml.gz',
                         segment=valid_seg,
                         directory=self.out_dir)
        
        # set the input and output files        
        node.add_output(insp, opts='output')
        node.add_input(cache_file, opts='frame-cache')    
        node.add_input(parent, opts='bank-file')                
        return node

class PyCBCInspiralExec(Executable):
    def __init__(self, exe_name):
        Executable.__init__(self, exe_name, 'vanilla')
    
    def create_job(self, cp, ifo, out_dir=None):
        return PyCBCInspiralJob(cp, self.exe_name, self.condor_universe, 
                                ifo=ifo, out_dir=out_dir)  
                                
    def get_valid_times(self, cp, ifo):
        pad_data = int(cp.get_opt_ifo(self.exe_name, 'pad-data', ifo))
        start_pad = int(cp.get_opt_ifo(self.exe_name, 'segment-start-pad', ifo))
        end_pad = int(cp.get_opt_ifo(self.exe_name, 'segment-end-pad', ifo))
        
        #FIXME this should not be hard coded (can be any integer > 
        #segment_length with pycbc_inspiral)
        data_length = 2048 + pad_data * 2
        start = pad_data + start_pad
        end = data_length - pad_data - end_pad
        return data_length, segments.segment(start, end)
        
class PyCBCTmpltbankJob(Job):
    def create_node(self, data_seg, valid_seg, parent=None, dfParents=None):
        node = LegacyAnalysisNode(self)
               
        if not dfParents or len(dfParents) != 1: 
            raise ValueError("%s must be supplied with a single cache file" 
                              %(self.exe_name)) 
                              
        pad_data = int(self.get_opt('pad-data'))                      
        if pad_data is None:
            raise ValueError("The option pad-data is a required option of "
                             "%s. Please check the ini file." % self.exe_name)
                             
        # set the remaining option flags
        node.set_start(data_seg[0] + pad_data)
        node.set_end(data_seg[1] - pad_data)       
        
        cache_file = dfParents[0]    
                            
        # FIXME add control for output type                       
        insp = AhopeFile(self.ifo, self.exe_name, extension='.xml.gz',
                         segment=valid_seg,
                         directory=self.out_dir)
        
        # set the input and output files      
        node.add_output(insp, opts='output-file')
        node.add_input(cache_file, opts='frame-cache')    
        return node

class PyCBCTmpltbankExec(Executable):
    def __init__(self, exe_name):
        Executable.__init__(self, exe_name, 'vanilla')
    
    def create_job(self, cp, ifo, out_dir=None):
        return PyCBCTmpltbankJob(cp, self.exe_name, self.condor_universe, 
                                 ifo=ifo, out_dir=out_dir)  
                                
    def get_valid_times(self, cp, ifo):
        pad_data = int(cp.get_opt_ifo(self.exe_name, 'pad-data', ifo))
        
        #FIXME this should not be hard coded 
        data_length = 2048 + pad_data * 2
        start = pad_data 
        end = data_length - pad_data
        return data_length, segments.segment(start, end)
    
class LegacySplitBankJob(Job):    
    def create_node(self, bank):
        """
        Set up a CondorDagmanNode class to run lalapps_splitbank code

        Parameters
        ----------
        bank : AhopeOutFile 
            The AhopeOutFile containing the template bank to be split

        Returns
        --------
        ode : Node
            The node to run the job
        """
        node = LegacyAnalysisNode(self)
        node.add_input(bank, opts='bank-file')
        
        # Get the output (taken from inspiral.py)
        url_list = []
        x = bank.path.split('-')
        num_banks = int(self.get_opt('number-of-banks'))
        for i in range( 0, num_banks):
            out_file = "%s-%s_%2.2d-%s-%s" %(x[0], x[1], i, x[2], x[3])
            out_url = urlparse.urlunparse(['file', 'localhost',
                                          os.path.join(self.out_dir, out_file),
                                          None, None, None])
            url_list.append(out_url)
        
        job_tag = bank.description + "_" + self.exe_name.upper()
        out_file_group = AhopeFile(bank.ifo, job_tag, bank.segment, 
                                   file_url=url_list)
        node.add_output(out_file_group)
        return node

