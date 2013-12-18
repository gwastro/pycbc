import os
import urlparse
from glue import pipeline
from glue import segments
from ahope_utils import Job, Node, Executable, AhopeFile

class LegacyAnalysisNode(Node, pipeline.AnalysisNode):
    pass
        
class LegacyAnalysisJob(Job):
    def create_node(self, data_seg, valid_seg, parents=None, dfparents=None):
        node = LegacyAnalysisNode(self)
        
        pad_data = int(self.get_opt('pad-data'))
        if pad_data is None:
            raise ValueError("The option pad-data is a required option of "
                             "%s. Please check the ini file." % self.exe_name)
                
        currNode.set_start(bankDataSeg[0] + pad_data)
        currNode.set_end(bankDataSeg[1] - pad_data)
        
        if not dfparents or len(dfparents) != 1: 
            raise ValueError("%s must be supplied with a single cache file" 
                              %(self.exe_name))   
                              
        extension = '.xml'
        gzipped = self.get_opt('write-compress')
        if gzipped:
            extension += '.gz'
              
        bank = AhopeFile(self.ifo, self.exe_name, 
                         extension=extension,
                         directory=self.out_dir)
        node.add_output(bank)
        node.add_input(cache_file, opt='cache-file')         
        return node
        
class LegacyInspiralJob(LegacyAnalysisJob):
    def create_node(self, data_seg, valid_seg, parents=None, dfparents=None):
        node = LegacyAnalysisJob.create_node(self, data_seg, valid_seg, parents, dfparents)
        self.set_trig_start(jobValidSeg[0])
        self.set_trig_end(jobValidSeg[1])        
        return node


class LegacyTmpltbankExec(Executable, LegacyValidTimes):
    def __init__(self, exe_name):
        if exe_name != 'tmpltbank':
            raise ValueError('lalapps_tmpltbank does not support setting '
                             'the exe_name to anything but "tmpltbank"')
                           
        Executable.__init__('tmpltbank', 'standard')

    def create_job(self, cp, ifo, out_dir=None):
        return LegacyAnalysisJob(cp, self.exe_name, self.condor_universe,
                                 ifo=ifo, out_dir=out_dir)
        
class LegacyInspiralExec(Executable, LegacyValidTimes):
    def __init__(self, exe_name):
        if exe_name != 'inspiral':
            raise ValueError('lalapps_tmpltbank does not support setting '
                             'the exe_name to anything but "inspiral"')
        Executable.__init__('inspiral', 'standard')

    def create_job(self, cp, ifo, out_dir=None):
        return LegacyInspiralJob(cp, self.exe_name, self.condor_universe, ifo=ifo, 
                                 out_dir=out_dir)

class LegacyValidTimes(object):
    def legacy_valid_times(cp, ifo):
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
        padData = int(cp.get_opt_ifo(self.exename, 'pad-data', ifo))
        self.padData = 8
        segmentLength = float(cp.get_opt_ifo(self.exename,
                                             'segment-length', ifo))
        sampleRate = float(cp.get_opt_ifo(self.exename,'sample-rate', ifo))
        numSegments = int(cp.get_opt_ifo(self.exename,
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
        if self.exename == 'inspiral':
            # Don't think inspiral will do well if segmentLength/4 is not
            # an integer
            validStart = validStart + int(segmentLength/(sampleRate * 4))
            validEnd = validEnd - int(segmentLength / (sampleRate * 4))
        validChunk = segments.segment([validStart,validEnd])

        return dataLength, validChunk

class legacy_splitbank_job_utils(Executable):
    """This class holds the function for lalapps_splitbank 
    usage following the old ihope specifications.
    """
    def __init__(self,exeName):
        self.exeName = exeName

    def create_condornode(self, ahopeDax, currJob, numBanks, parent):
        """
        Set up a CondorDagmanNode class to run lalapps_splitbank code

        Parameters
        ----------
        ahopeDax : pipeline.CondorDAG instance
            The workflow to hold of the ahope jobs.
        currJob : pipeline.CondorDagmanJob
            The CondorDagmanJob to use when setting up the individual nodes.
        numBanks : int
            Number of parts to split template bank into.
        parent : AhopeOutFile (optional, kwarg, default=None)
            The AhopeOutFile containing the job that is parent to the one being
            set up.

        Returns
        --------
        tmpltBankNode : pipeline.CondorDagmanNode
            The node to run the job
        list of strings
            The output files
        """
        currNode = LegacyInspiralAnalysisNode(currJob)
        currNode.set_category(self.exeName)
        # Does this need setting?: currNode.set_priority(?)
        currNode.set_bank(parent.path)
        # Set the number of banks
        currNode.add_var_opt('number-of-banks',numBanks)
        # Get the output (taken from inspiral.py)
        outUrlList = []
        x = parent.path.split('-')
        for i in range( 0, numBanks ):
            outFile = "%s-%s_%2.2d-%s-%s" %(x[0], x[1], i, x[2], x[3])
            outUrl = urlparse.urlunparse(['file', 'localhost',\
                                          os.path.join(self.outDir, outFile),\
                                          None, None, None])
            outUrlList.append(outUrl)
        parentJob = parent.job
        if parentJob:
            currNode.add_parent(parentJob)
        currNode.finalize()
        ahopeDax.add_node(currNode)

        return currNode, outUrlList

