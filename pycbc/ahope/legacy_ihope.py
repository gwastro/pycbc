import os
import urlparse
from glue import pipeline
from glue import segments
from ahope_utils import Job, Node, Executable, AhopeFile

def legacy_get_valid_times(self):
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
    padData = int(self.get_opt('pad-data'))
    segmentLength = float(self.get_opt('segment-length'))
    sampleRate = float(self.get_opt('sample-rate'))
    numSegments = int(self.get_opt('number-of-segments'))
    
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
    # FIXME: This should probably be pulled into the Node class. It is not
    # specific to the Legacy analysis codes
    set_jobnum_tag = pipeline.AnalysisNode.set_user_tag
    
        
class LegacyAnalysisJob(Job):
    def __init__(self, cp, exe_name, universe, ifo=None, out_dir=None):
        Job.__init__(self, cp, exe_name, universe, ifo, out_dir)

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
        node.add_input(cache_file, opt='frame-cache')         
        return node

    get_valid_times = legacy_get_valid_times
        
class LegacyInspiralJob(LegacyAnalysisJob):
    def __init__(self, cp, exe_name, universe, ifo=None, injection_file=None, 
                       out_dir=None):
        LegacyAnalysisJob.__init__(cp, exe_name, universe, ifo, out_dir)
        self.injection_file = injection_file 

    def create_node(self, data_seg, valid_seg, parent=None, dfParents=None):
        node = LegacyAnalysisJob.create_node(self, data_seg, valid_seg, 
                                                   parent, dfParents)
        node.set_trig_start(valid_seg[0])
        node.set_trig_end(valid_seg[1])  
        node.add_input(parent, opt='bank-file')    
        
        if self.injection_file is not None:
            node.add_input(self.injection_file, 'injection-file')
        return node

class LegacyTmpltbankExec(Executable):
    def __init__(self, exe_name):
        if exe_name != 'tmpltbank':
            raise ValueError('lalapps_tmpltbank does not support setting '
                             'the exe_name to anything but "tmpltbank"')                           
        Executable.__init__(self, 'tmpltbank')

    def create_job(self, cp, ifo, out_dir=None):
        return LegacyAnalysisJob(cp, self.exe_name, self.condor_universe,
                                 ifo=ifo, out_dir=out_dir)   
        
class LegacyInspiralExec(Executable):
    def __init__(self, exe_name):
        if exe_name != 'inspiral':
            raise ValueError('lalapps_tmpltbank does not support setting '
                             'the exe_name to anything but "inspiral"')
        Executable.__init__(self, 'inspiral')

    def create_job(self, cp, ifo, out_dir=None):
        return LegacyInspiralJob(cp, self.exe_name, self.condor_universe, 
                                 ifo=ifo, 
                                 out_dir=out_dir)

class LegacySplitBankExec(Executable):
    """This class holds the function for lalapps_splitbank 
    usage following the old ihope specifications.
    """
    def create_job(self, cp, ifo, out_dir=None):
        return LegacySplitBankJob(cp, self.exe_name, self.condor_universe, 
                                  ifo=ifo,
                                  out_dir=out_dir)

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
        node : Node
            The node to run the job
        """
        node = LegacyAnalysisNode(self)
        node.add_input(bank, opt='bank-file')
        
        # Get the output (taken from inspiral.py)
        url_list = []
        x = bank.filename.split('-')
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

