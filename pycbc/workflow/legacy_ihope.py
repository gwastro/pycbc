# Copyright (C) 2013  Ian Harry
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
This library code contains functions and classes that are used to set up
and add jobs/nodes from legacy lalapps C-code to a pycbc workflow.
For details about pycbc.workflow see here:
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/ahope.html
"""

import os
import urlparse
from glue import segments
from pycbc.workflow.core import Executable, File, FileList, Node

def legacy_get_valid_times(self):
    """
    Return the length of data that the tmpltbank job will need to read and
    the part of that data that the template bank is valid for. In the case
    of lalapps_tmpltbank the following options are needed to set this up
    and will be used by the Executable to figure this out:

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

    Returns
    -------
    dataLength : float (seconds)
        The length of data that the job will need
    validChunk : glue.glue.segments.segment
        The start and end of the dataLength that is valid for the template
        bank.
    """
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
    if self.name == 'inspiral':
        # Don't think inspiral will do well if segmentLength/4 is not
        # an integer
        validStart = validStart + int(segmentLength/(sampleRate * 4))
        validEnd = validEnd - int(segmentLength / (sampleRate * 4))
    validChunk = segments.segment([validStart,validEnd])

    return [dataLength], [validChunk]

        
class LegacyAnalysisNode(Node):
    # This is *ONLY* used by legacy codes where pycbc.workflow cannot directly
    # set the output name. Do not use elsewhere!
    def set_jobnum_tag(self, num):
        self.add_opt('--user-tag', num)
        
class LegacyAnalysisExecutable(Executable):
    """
    The class responsible for setting up jobs for legacy lalapps C-code
    Executables.
    """
    current_retention_level = Executable.CRITICAL
    def __init__(self, cp, name, universe=None, ifo=None, tags=[], out_dir=None):
        super(LegacyAnalysisExecutable, self).__init__(cp, name, universe, ifo, out_dir, tags=tags)

    def create_node(self, data_seg, valid_seg, parent=None, dfParents=None, tags=[]):
        node = LegacyAnalysisNode(self)
        
        if not dfParents: 
            raise ValueError("%s must be supplied with frame files" 
                              %(self.name))  
        
        pad_data = int(self.get_opt('pad-data'))
        if pad_data is None:
            raise ValueError("The option pad-data is a required option of "
                             "%s. Please check the ini file." % self.name)                                     
          
        # hide import here to avoid circular import
        from pycbc.workflow import int_gps_time_to_str
        node.add_opt('--gps-start-time', int_gps_time_to_str(data_seg[0] + pad_data))
        node.add_opt('--gps-end-time', int_gps_time_to_str(data_seg[1] - pad_data))   
         
        cache_file = dfParents[0]       
        
        #check the extension       
        extension = '.xml'
        gzipped = self.has_opt('write-compress')
        if gzipped is not None:
            extension += '.gz'
        
        #create the output file for this job 
        out_file = File(self.ifo, self.name, valid_seg,
                             extension=extension,
                             directory=self.out_dir,
                             tags=self.tags + tags,
                             store_file=self.retain_files)
 
        node.add_output_opt('--output-file', out_file)
        node.add_input_list_opt('--frame-files', dfParents)
        return node

    get_valid_times = legacy_get_valid_times
    
class LegacyTmpltbankExecutable(LegacyAnalysisExecutable):
    current_retention_level = Executable.CRITICAL
        
class LegacyInspiralExecutable(LegacyAnalysisExecutable):
    """
    The class responsible for setting up jobs for legacy lalapps_inspiral
    Executable.
    """
    current_retention_level = Executable.CRITICAL
    def __init__(self, cp, name, universe=None, ifo=None, injection_file=None, 
                       gate_files=None, out_dir=None, tags=[]):
        super(LegacyInspiralExecutable, self).__init__(cp, name, universe, ifo, 
                                    out_dir=out_dir, tags=tags)
        self.injection_file = injection_file 

    def create_node(self, data_seg, valid_seg, parent=None, dfParents=None, tags=[]):
        node = LegacyAnalysisExecutable.create_node(self, data_seg, valid_seg, 
                                                   parent, dfParents, tags=tags)
        node.add_opt('--trig-start-time', valid_seg[0])
        node.add_opt('--trig-end-time', valid_seg[1])  
        node.add_input_opt('--bank-file', parent)    
        
        if self.injection_file is not None:
            node.add_input_opt('--injection-file', self.injection_file)
        return node
    
class LegacySplitBankExecutable(Executable):    
    """
    The class responsible for creating jobs for lalapps_splitbank.
    """
    current_retention_level = Executable.NON_CRITICAL
    def __init__(self, cp,name, universe, numBanks,
                 ifo=None, out_dir=None, tags=[]):
        Job.__init__(self, cp, name, universe, ifo, out_dir, tags=tags)
        self.num_banks = int(numBanks)
        tmpNumBanks = self.get_opt("number-of-banks")
        if tmpNumBanks:
            # Option is in ini file, check it agrees
            if not int(tmpNumBanks) == self.num_banks:
                errMsg = "Number of banks provided in [workflow-splitbank] "
                errMsg += "section does not agree with value given in the "
                errMsg += "[%s] section of the ini file. " %(name)
                raise ValueError(errMsg)
        else:
            # Option not given, so add it
           self.add_opt("number-of-banks", self.num_banks)
            

    def create_node(self, bank):
        """
        Set up a CondorDagmanNode class to run lalapps_splitbank code

        Parameters
        ----------
        bank : pycbc.workflow.core.File 
            The OutFile containing the template bank to be split

        Returns
        --------
        node : pycbc.workflow.core.Node
            The node to run the job
        """
        node = Node(self)
        # FIXME: This is a hack because SplitBank fails if given an input file
        # whose path contains the character '-' or if the input file is not in
        # the same directory as the output. Therefore we just set the path to
        # be the local path
        fullPath = bank.cache_entry.path
        bank.cache_entry.path = os.path.basename(fullPath)
        node.add_input_opt('--bank-file', bank)
        # FIXME: Set the path back to what it was. This is part of the hack
        #        above and should be removed if possible.
        bank.cache_entry.path = fullPath
        
        # Get the output (taken from inspiral.py)
        url_list = []
        x = bank.filename.split('-')
        if len(x) != 4:
            errMsg = "Input file name is not compatible with splitbank. Name "
            errMsg += "must follow the lal cache standard, for example "
            errMsg += "H1-TMPLTBANK-900000000-1000.xml. "
            errMsg += "Got %s." %(bank.filename,)
            raise ValueError(errMsg)
        for i in range( 0, self.num_banks):
            out_file = "%s-%s_%2.2d-%s-%s" %(x[0], x[1], i, x[2], x[3])
            out_url = urlparse.urlunparse(['file', 'localhost',
                                          os.path.join(self.out_dir, out_file),
                                          None, None, None])
            url_list.append(out_url)
                
            job_tag = bank.description + "_" + self.name.upper()
            out_file = File(bank.ifo, job_tag, bank.segment, 
                            file_url=out_url, tags=bank.tags,
                            store_file=self.retain_files)
            node._add_output(out_file)
        return node


class LegacyCohPTFInspiralExecutable(LegacyAnalysisExecutable):
    """
    The class responsible for setting up jobs for legacy
    lalapps_coh_PTF_inspiral executable.
    """
    current_retention_level = Executable.CRITICAL
    def __init__(self, cp, name, universe=None, ifo=None, injection_file=None,
                 gate_files=None, out_dir=None, tags=[]):
        super(LegacyCohPTFInspiralExecutable, self).__init__(cp, name, universe,
                ifo, out_dir=out_dir, tags=tags)
        self.cp = cp
        self.injection_file = injection_file
        self.data_seg = segments.segment(int(cp.get('workflow', 'start-time')),
                                         int(cp.get('workflow', 'end-time')))
        self.num_threads = 1
 
    def create_node(self, data_seg, valid_seg, parent=None, inj_file=None,
                    dfParents=None, bankVetoBank=None, tags=[]):
        node = Node(self)

        if not dfParents:
            raise ValueError("%s must be supplied with frame files"
                              % self.name)

        pad_data = self.get_opt('pad-data')
        if pad_data is None:
            raise ValueError("The option pad-data is a required option of "
                             "%s. Please check the ini file." % self.name)
        
        # Feed in bank_veto_bank.xml
        if self.cp.has_option('inspiral', 'do-bank-veto'):
            if not bankVetoBank:
                raise ValueError("%s must be given a bank veto file if the"
                                 "argument 'do-bank-veto' is given"
                                 % self.name)
            node.add_input_opt('--bank-veto-templates', bankVetoBank)
        
        node.add_opt('--gps-start-time', data_seg[0] + int(pad_data))
        node.add_opt('--gps-end-time', data_seg[1] - int(pad_data))
        node.add_opt('--trig-start-time', valid_seg[0])
        node.add_opt('--trig-end-time', valid_seg[1])

        node.add_profile('condor', 'request_cpus', self.num_threads)

        # Set the input and output files
        node.new_output_file_opt(data_seg, '.xml.gz', '--output-file',
                                 tags=tags, store_file=self.retain_files)
        node.add_input_opt('--non-spin-bank', parent, )

        for frameCache in dfParents:
            node.add_input_opt('--%s-frame-cache' % frameCache.ifo.lower(),
                               frameCache)
            node.add_arg('--%s-data' % frameCache.ifo.lower())
            node.add_opt('--%s-channel-name' % frameCache.ifo.lower(),
                         self.cp.get('workflow',
                                     '%s-channel-name' %frameCache.ifo.lower()))

        if inj_file is not None:
            node.add_input_opt('--injection-file', inj_file)

        return node

    def get_valid_times(self):
        overlap = int(self.get_opt('segment-duration')) / 4
        pad_data = int(self.get_opt('pad-data'))
        
        valid_start = self.data_seg[0] + pad_data + overlap
        valid_end = self.data_seg[1] - pad_data - overlap
        return self.data_seg, segments.segment(valid_start, valid_end)

class LegacyCohPTFTrigCombiner(LegacyAnalysisExecutable):
    """
    The class responsible for setting up jobs for legacy coh_PTF_trig_combiner
    executable.
    """
    current_retention_level = Executable.INTERMEDIATE_PRODUCT
    def __init__(self, cp, name, universe=None, ifo=None, injection_file=None,
                 out_dir=None, tags=[]):
        super(LegacyCohPTFTrigCombiner, self).__init__(cp, name, universe,
              ifo=ifo, out_dir=out_dir, tags=tags)
        self.cp = cp
        self.ifos = ifo
        self.num_threads = 1

    def create_node(self, trig_files=None, segment_dir=None, out_tags=[],
                    tags=[]):
        node = Node(self)

        if not trig_files:
            raise ValueError("%s must be supplied with trigger files"
                              % self.name)

        # Data options
        pad_data = self.cp.get('inspiral', 'pad-data')
        if pad_data is None:
            raise ValueError("The option pad-data is a required option of "
                             "%s. Please check the ini file." % self.name)

        num_trials = int(self.cp.get("trig_combiner", "num-trials"))
        trig_name = self.cp.get('workflow', 'trigger-name')
        node.add_opt('--grb-name', trig_name)
        
        node.add_opt('--pad-data', pad_data)
        node.add_opt('--segment-length', self.cp.get('inspiral',
                                                     'segment-duration'))
        node.add_opt('--ifo-tag', self.ifos)
        node.add_opt('--user-tag', 'INSPIRAL')

        # Set input / output options
        node.add_input_list_opt('--input-files', trig_files)

        node.add_opt('--segment-dir', segment_dir)
        node.add_opt('--output-dir', self.out_dir)

        out_files = FileList([])
        for out_tag in out_tags:
            out_file = File(self.ifos, 'INSPIRAL', trig_files[0].segment,
                            directory=self.out_dir, extension='xml.gz',
                            tags=["GRB%s" % trig_name, out_tag],
                            store_file=self.retain_files)
            #out_file.PFN(out_file.cache_entry.path, site="local")
            out_files.append(out_file)

        for trial in range(1, num_trials + 1):
            out_file = File(self.ifos, 'INSPIRAL', trig_files[0].segment,
                            directory=self.out_dir, extension='xml.gz',
                            tags=["GRB%s" % trig_name, "OFFTRIAL_%d" % trial],
                            store_file=self.retain_files)
            #out_file.PFN(out_file.cache_entry.path, site="local")
            out_files.append(out_file)

        node.add_profile('condor', 'request_cpus', self.num_threads)

        return node, out_files


class LegacyCohPTFTrigCluster(LegacyAnalysisExecutable):
    """
    The class responsible for setting up jobs for legacy coh_PTF_trig_cluster
    executable.
    """
    current_retention_level = Executable.CRITICAL
    def __init__(self, cp, name, universe=None, ifo=None, injection_file=None,
                 out_dir=None, tags=[]):
        super(LegacyCohPTFTrigCluster, self).__init__(cp, name, universe,
              ifo=ifo, out_dir=out_dir, tags=tags)
        self.cp = cp
        self.ifos = ifo
        self.num_threads = 1
 
    def create_node(self, parent, tags=[]):
        node = Node(self)

        # Set input / output options
        node.add_opt('--trig-file', '%s' % parent.storage_path)
        node.add_opt('--output-dir', self.out_dir)

        node.add_profile('condor', 'request_cpus', self.num_threads)

        # Adding output files as pycbc.workflow.core.File objects
        out_file = File(self.ifos, 'INSPIRAL', parent.segment,
                        directory=self.out_dir, extension='xml.gz',
                        tags=[parent.tag_str, 'CLUSTERED'],
                        store_file=self.retain_files)
        #out_file.PFN(out_file.cache_entry.path, site="local")

        return node, FileList([out_file])


class LegacyCohPTFInjfinder(LegacyAnalysisExecutable):
    """
    The class responsible for setting up jobs for legacy coh_PTF_injfinder
    executable.
    """
    current_retention_level = Executable.CRITICAL
    def __init__(self, cp, name, universe=None, ifo=None, injection_file=None,
                 out_dir=None, tags=[]):
        super(LegacyCohPTFInjfinder, self).__init__(cp, name, universe,
              ifo=ifo, out_dir=out_dir, tags=tags)
        self.cp = cp
        self.ifos = ifo
        self.num_threads = 1

    def create_node(self, trig_files, inj_files, seg_dir, tags=[]):
        node = Node(self)

        # Set input / output options
        node.add_input_list_opt('--input-files', trig_files)
        node.add_input_list_opt('--inj-files', inj_files)

        node.add_opt('--ifo-tag', self.ifos)
        node.add_opt('--exclude-segments', '%s/bufferSeg.txt' % seg_dir)
        node.add_opt('--output-dir', self.out_dir)

        # Create output files as File objects
        name_string = inj_files[0].description
        seg = trig_files[0].segment

        f_file = File(self.ifos, name_string, seg, extension="xml",
                      directory=self.out_dir, store_file=self.retain_files,
                      tags=[inj_files[0].tag_str.replace("split0", "FOUND")])

        m_file = File(self.ifos, name_string, seg, extension="xml",
                      directory=self.out_dir, store_file=self.retain_files,
                      tags=[inj_files[0].tag_str.replace("split0", "MISSED")])

        return node, FileList([f_file, m_file])


class LegacyCohPTFInjcombiner(LegacyAnalysisExecutable):
    """
    The class responsible for setting up jobs for legacy coh_PTF_injcombiner
    executable.
    """
    current_retention_level = Executable.CRITICAL
    def __init__(self, cp, name, universe=None, ifo=None, injection_file=None,
                 out_dir=None, tags=[]):
        super(LegacyCohPTFInjcombiner, self).__init__(cp, name, universe,
              ifo=ifo, out_dir=out_dir, tags=tags)
        self.cp = cp
        self.ifos = ifo
        self.num_threads = 1

    def create_node(self, parent, inj_trigs, inj_string, max_inc, segment):
        node = Node(self)

        trig_name = self.cp.get('workflow', 'trigger-name')
        node.add_opt('--inj-string', inj_string)
        node.add_opt('--max-inclination', max_inc)
        node.add_opt('--inj-cache', '%s' % parent.storage_path)

        out_files = FileList([])
        for inj_trig in inj_trigs:
            out_file_tag = [inj_string, "FILTERED", max_inc,
                            inj_trig.tag_str.rsplit('_', 1)[-1]]
            out_file = File(self.ifos, inj_trig.description,
                            inj_trig.segment, extension="xml",
                            directory=self.out_dir, tags=out_file_tag)
            out_file.PFN(out_file.cache_entry.path, site="local")
            out_files.append(out_file)

        node.add_opt('--output-dir', self.out_dir)

        return node, out_files


class LegacyCohPTFSbvPlotter(LegacyAnalysisExecutable):
    """
    The class responsible for setting up jobs for legacy coh_PTF_sbv_plotter
    executable.
    """
    current_retention_level = Executable.CRITICAL
    def __init__(self, cp, name, universe=None, ifo=None, injection_file=None,
                 out_dir=None, tags=[]):
        super(LegacyCohPTFSbvPlotter, self).__init__(cp, name, universe,
              ifo=ifo, out_dir=out_dir, tags=tags)
        self.cp = cp
        self.ifos = ifo
        self.num_threads = 1

    def create_node(self, parent=None, seg_dir=None, inj_file=None, tags=[]):
        node = Node(self)

        if not parent:
            raise ValueError("%s must be supplied with trigger files"
                             % self.name)

        if isinstance(parent, str) and tags[1]=="_unclustered":
            node.add_opt('--trig-file', '%s' % parent.storage_path)
            tags[1] = ""
        else:
            node.add_opt('--trig-file', '%s' % parent.storage_path)

        node.add_opt('--grb-name', self.cp.get('workflow', 'trigger-name'))

        # Set input / output options
        node.add_opt('--veto-directory', seg_dir)
        node.add_opt('--segment-dir', seg_dir)
        out_dir = "%s/output/%s/plots%s" % (self.out_dir, tags[0], tags[1])
        node.add_opt('--output-path', out_dir)

        if inj_file is not None:
            node.add_opt('--inj-file', inj_file.storage_path)

        node.add_profile('condor', 'request_cpus', self.num_threads)

        return node


class LegacyCohPTFEfficiency(LegacyAnalysisExecutable):
    """
    The class responsible for setting up jobs for legacy coh_PTF_efficiency
    executable.
    """
    current_retention_level = Executable.CRITICAL
    def __init__(self, cp, name, universe=None, ifo=None, injection_file=None,
                 out_dir=None, tags=[]):
        super(LegacyCohPTFEfficiency, self).__init__(cp, name, universe,
              ifo=ifo, out_dir=out_dir, tags=tags)
        self.cp = cp
        self.ifos = ifo
        self.num_threads = 1

    def create_node(self, parent=None, offsource_file=None, seg_dir=None,
                    found_file=None, missed_file=None, tags=[]):
        node = Node(self)

        if not parent:
            raise ValueError("%s must be supplied with trigger files"
                             % self.name)

        # Set input / output options
        node.add_opt('--onsource-file', '%s' % parent.storage_path)

        node.add_opt('--offsource-file', '%s' % offsource_file.storage_path)
        node.add_opt('--veto-directory', seg_dir)
        node.add_opt('--segment-dir', seg_dir)
        
        if found_file and missed_file:
            node.add_opt('--found-file', '%s' % found_file.storage_path)
            node.add_opt('--missed-file', '%s' % missed_file.storage_path)
            out_dir = "%s/output/%s/efficiency_%s" % (self.out_dir, tags[1],
                                                      tags[0])
            if self.cp.has_option_tag('injections', 'min-distance', tags[-1]):
                lower_dist = float(self.cp.get_opt_tag('injections',
                        'min-distance', tags[-1]))
                #Convert distance from kpc to Mpc then add as option
                lower_dist /= 1e3
                node.add_opt('--lower-inj-dist', lower_dist)
            if self.cp.has_option_tag('injections', 'max-distance', tags[-1]):
                upper_dist = float(self.cp.get_opt_tag('injections',
                        'max-distance', tags[-1]))
                #Convert distance from kpc to Mpc then add as option
                upper_dist /= 1e3
                node.add_opt('--upper-inj-dist', upper_dist)

        elif found_file or missed_file:
            if found_file:
                present = found_file
            else:
                present = missed_file
            raise ValueError("Must either be supplied with no injection files "
                             "or both missed and found injection files. "
                             "Received only %s" % present.name)
        else:
            out_dir = "%s/output/%s/efficiency" % (self.out_dir, tags[0])

        node.add_opt('--output-path', out_dir)
        node.add_profile('condor', 'request_cpus', self.num_threads)

        return node

class PyGRBMakeSummaryPage(LegacyAnalysisExecutable):
    """
    The class responsible for setting up the summary page generation job for
    the PyGRB workflow.
    """
    current_retention_level = Executable.CRITICAL
    def __init__(self, cp, name, universe=None, ifo=None, injection_file=None,
                 out_dir=None, tags=[]):
        super(PyGRBMakeSummaryPage, self).__init__(cp, name, universe, ifo=ifo,
              out_dir=out_dir, tags=tags)
        self.cp = cp
        self.ifos = ifo
        self.output_dir = out_dir
        self.num_threads = 1

    def create_node(self, parent=None, c_file=None, open_box=False,
                    tuning_tags=None, exclusion_tags=None, html_dir=None,
                    tags=[]):
        node = Node(self)

        node.add_opt('--grb-name', self.cp.get('workflow', 'trigger-name'))
        node.add_opt('--start-time', self.cp.get('workflow', 'trigger-time'))
        node.add_opt('--ra', self.cp.get('workflow', 'ra'))
        node.add_opt('--dec', self.cp.get('workflow', 'dec'))
        node.add_opt('--ifo-tag', self.ifos)

        if tuning_tags is not None:
            node.add_opt('--tuning-injections', ','.join(tuning_tags))

        if exclusion_tags is not None:
            node.add_opt('--exclusion-injections', ','.join(exclusion_tags))

        if open_box:
            node.add_opt('--open-box')

        if html_dir is not None:
            node.add_opt('--html-path', html_dir)

        # Set input / output options
        node.add_opt('--config-file', '%s' % c_file.storage_path) 
        node.add_opt('--output-path', "%s/output" % self.output_dir)

        node.add_profile('condor', 'request_cpus', self.num_threads)

        return node

