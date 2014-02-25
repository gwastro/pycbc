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
and add jobs/nodes to an ahope workflow. For details about ahope see here:
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/ahope.html
"""

import math
from glue import segments, pipeline
from pycbc.ahope.ahope_utils import *
from pycbc.ahope.legacy_ihope import *

def select_tmpltbankjob_instance(curr_exe, curr_section):
    """
    This function returns an instance of the class that is appropriate for
    creating a template bank within ihope.
    
    Parameters
    ----------
    curr_exe : string
        The name of the executable that is being used.
    curr_section : string
        The name of the section storing options for this executble

    Returns
    --------
    Instanced class : exe_class
        An instance of the class that holds the utility functions appropriate
        for the given executable. This class **must** contain
        * exe_class.create_job()
        and the job returned by this **must** contain
        * job.get_valid_times(ifo, )
        * job.create_node()
    """
    # This is basically a list of if statements

    if curr_exe == 'lalapps_tmpltbank':
        exe_class = LegacyTmpltbankExec(curr_section)
    elif curr_exe == 'pycbc_geom_nonspinbank':
        exe_class = PyCBCTmpltbankExec(curr_section)
    else:
        # Should we try some sort of default class??
        err_string = "No class exists for executable %s" %(curr_exe,)
        raise NotImplementedError(err_string)
    return exe_class

def select_matchedfilterjob_instance(curr_exe, curr_section):
    """
    This function returns an instance of the class that is appropriate for
    matched-filtering within ahope.
    
    Parameters
    ----------
    curr_exe : string
        The name of the executable that is being used.
    curr_section : string
        The name of the section storing options for this executble

    Returns
    --------
    Instanced class : exe_class
        An instance of the class that holds the utility functions appropriate
        for the given executable. This class **must** contain
        * exe_class.create_job()
        and the job returned by this **must** contain
        * job.get_valid_times(ifo, )
        * job.create_node()
    """
    # This is basically a list of if statements
    if curr_exe == 'lalapps_inspiral':
        exe_class = LegacyInspiralExec(curr_section)
    elif curr_exe == 'pycbc_inspiral':
        exe_class = PyCBCInspiralExec(curr_section)
    else:
        # Should we try some sort of default class??
        err_string = "No class exists for executable %s" %(curr_exe,)
        raise NotImplementedError(err_string)
        
    return exe_class

def select_splitfilejob_instance(curr_exe, curr_section):
    """
    This function returns an instance of the class that is appropriate for
    splitting an output file up within ahope (for e.g. splitbank).
    
    Parameters
    ----------
    curr_exe : string
        The name of the executable that is being used.
    curr_section : string
        The name of the section storing options for this executble

    Returns
    --------
    Instanced class : exe_class
        An instance of the class that holds the utility functions appropriate
        for the given executable. This class **must** contain
        * exe_class.create_job()
        and the job returned by this **must** contain
        * job.create_node()
    """
    # This is basically a list of if statements
    if curr_exe == 'lalapps_splitbank':
        exe_class = LegacySplitBankExec(curr_section)
    # Some elif statements
    else:
        # Should we try some sort of default class??
        err_string = "No class exists for executable %s" %(curr_exe,)
        raise NotImplementedError(errString)

    return exe_class

def sngl_ifo_job_setup(workflow, ifo, out_files, curr_exe_job, science_segs, 
                       datafind_outs, output_dir, parents=None, 
                       link_job_instance=None, allow_overlap=True):
    """
    This function sets up a set of single ifo jobs.

    Parameters
    -----------
    workflow: ahope.Workflow
        An instanced class that manages the constructed workflow.
    ifo : string
        The name of the ifo to set up the jobs for
    out_files : ahope.AhopeFileList
        The AhopeFileList containing the list of jobs. Jobs will be appended
        to this list, and it does not need to be empty when supplied.
    curr_exe_job : ahope.Job
        An instanced of the PyCBC Job class that has a get_valid times method.
    science_segs : glue.segments.segmentlist
        The list of times that the jobs should cover
    datafind_outs : ahope.AhopeFileList
        The file list containing the datafind files.
    output_dir : path string
        The directory where data products will be placed.
    parents : ahope.AhopeFileList (optional, kwarg, default=None)
        The AhopeFileList containing the list of jobs that are parents to
        the one being set up.
    link_job_instance : Job instance (optional),
        Coordinate the valid times with another executable.
    allow_overlap : boolean (optional, kwarg, default = True)
        If this is set the times that jobs are valid for will be allowed to
        overlap. This may be desired for template banks which may have some
        overlap in the times they cover. This may not be desired for inspiral
        jobs, where you probably want triggers recorded by jobs to not overlap
        at all.

    Returns
    --------
    out_files : ahope.AhopeFileList
        A list of the files that will be generated by this step in the
        ahope workflow.
    """
    cp = workflow.cp
    
    # Set up the condorJob class for the current executable
    data_length, valid_chunk = curr_exe_job.get_valid_times()
    
    # Begin by getting analysis start and end, and start and end of time
    # that the output file is valid for
    valid_length = abs(valid_chunk)

    data_chunk = segments.segment([0, data_length])
    job_tag = curr_exe_job.exe_name.upper()
    
    if link_job_instance:
        # FIXME: Should we remove this, after testing is complete??
        # EURGHH! What we are trying to do here is, if this option is given,
        # line up the template bank and inspiral jobs so that there is one
        # template bank for each inspiral job. This is proving a little messy
        # and probably still isn't perfect.

        # What data does the linked exe use?
        link_data_length,link_valid_chunk = link_job_instance.get_valid_times()
        # What data is lost at start and end from either job?
        start_data_loss = max(valid_chunk[0], link_valid_chunk[0])
        end_data_loss = max(data_length - valid_chunk[1],\
                            link_data_length - link_valid_chunk[1])
        # Correct the data for both jobs
        valid_chunk = segments.segment(start_data_loss, \
                                       data_length - end_data_loss)
        link_valid_chunk = segments.segment(start_data_loss, \
                                       link_data_length - end_data_loss)
        valid_length = abs(valid_chunk)
        link_valid_length = abs(link_valid_chunk)

        # Which one is now longer? Use this is valid_length
        if link_valid_length < valid_length:
            valid_length = link_valid_length

    # DO NOT! use valid length here, as valid_length and abs(valid_chunk)
    # may be different by this point if using link_exe_instance
    data_loss = data_length - abs(valid_chunk)

    
    if data_loss < 0:
        raise ValueError("Ahope needs fixing! Please contact a developer")
        
    # Loop over science segments and set up jobs
    for curr_seg in science_segs:
        # Is there enough data to analyse?
        curr_seg_length = abs(curr_seg)
        if curr_seg_length < data_length:
            continue
        # How many jobs do we need
        curr_seg_length = float(abs(curr_seg))
        num_jobs = int( math.ceil( \
                 (curr_seg_length - data_loss) / float(valid_length) ))
        # What is the incremental shift between jobs
        time_shift = (curr_seg_length - data_length) / float(num_jobs - 1)
        for job_num in range(num_jobs):
            # Get the science segment for this job
            # small factor of 0.0001 to avoid float round offs causing us to
            # miss a second at end of segments.
            shift_dur = curr_seg[0] + int(time_shift * job_num + 0.0001)
            job_data_seg = data_chunk.shift(shift_dur)
            # Sanity check that all data is used
            if job_num == 0:
               if job_data_seg[0] != curr_seg[0]:
                   errMsg = "Job is not using data from the start of the "
                   errMsg += "science segment. It should be using all data."
                   raise ValueError(errMsg)
            if job_num == (num_jobs - 1):
                if job_data_seg[1] != curr_seg[1]:
                    errMsg = "Job is not using data from the end of the "
                    errMsg += "science segment. It should be using all data."
                    raise ValueError(errMsg)
            job_valid_seg = valid_chunk.shift(shift_dur)
            # If we need to recalculate the valid times to avoid overlap
            if not allow_overlap:
                data_per_job = (curr_seg_length - data_loss) / float(num_jobs)
                lower_boundary = job_num*data_per_job + \
                                     valid_chunk[0] + curr_seg[0]
                upper_boundary = data_per_job + lower_boundary
                # NOTE: Convert to int after calculating both boundaries
                # small factor of 0.0001 to avoid float round offs causing us to
                # miss a second at end of segments.
                lower_boundary = int(lower_boundary)
                upper_boundary = int(upper_boundary + 0.0001)
                if lower_boundary < job_valid_seg[0] or \
                        upper_boundary > job_valid_seg[1]:
                    err_msg = ("Ahope is attempting to generate output "
                              "from a job at times where it is not valid.")
                    raise ValueError(err_msg)
                job_valid_seg = segments.segment([lower_boundary, 
                                                  upper_boundary])
                
            if parents:
                # Find the set of files with the best overlap
                curr_parent = parents.find_outputs_in_range(ifo, job_valid_seg)
                if not curr_parent:
                    err_string = ("No parent jobs found overlapping %d to %d." 
                                  %(job_valid_seg[0], job_valid_seg[1]))
                    err_string += "\nThis is a bad error! Contact a developer."
                    raise ValueError(err_string)
            else:
                curr_parent = [None]

            if datafind_outs:
                curr_dfouts = datafind_outs.find_all_output_in_range(ifo, 
                                                                  job_data_seg)
                if not curr_dfouts:
                    err_str = ("No datafind jobs found overlapping %d to %d."
                                %(job_data_seg[0],job_data_seg[1]))
                    err_str += "\nThis shouldn't happen. Contact a developer."
                    raise ValueError(err_str)

            for pnum, parent in enumerate(curr_parent):
                if pnum != 0:
                    tag = [str(pnum)]
                else:
                    tag = []
                # To ensure output file uniqueness I add a tag
                # We should generate unique names automatically, but it is a 
                # pain until we can set the output names for all executables              
                node = curr_exe_job.create_node(job_data_seg, job_valid_seg, 
                                                parent=parent,
                                                dfParents=curr_dfouts, 
                                                tags=tag)
                workflow.add_node(node)
                out_files += node.output_files
    return out_files

class PyCBCInspiralJob(Job):
    """
    The class used to create jobs for pycbc_inspiral executable.
    """
    def __init__(self, cp, exe_name, universe, ifo=None, out_dir=None, injection_file=None, tags=[]):
        Job.__init__(self, cp, exe_name, universe, ifo, out_dir, tags=tags)
        self.cp = cp
        self.set_memory(2000)
        self.injection_file = injection_file

        if self.get_opt('processing-scheme') == 'cuda':
            self.needs_gpu()

    def create_node(self, data_seg, valid_seg, parent=None, dfParents=None, tags=[]):
        node = Node(self)
        pad_data = int(self.get_opt('pad-data'))
        if pad_data is None:
            raise ValueError("The option pad-data is a required option of "
                             "%s. Please check the ini file." % self.exe_name)

        if not dfParents or len(dfParents) != 1:
            raise ValueError("%s must be supplied with a single cache file"
                              %(self.exe_name))

        # set remaining options flags   
        node.add_var_opt('gps-start-time', data_seg[0] + pad_data)
        node.add_var_opt('gps-end-time', data_seg[1] - pad_data)
        node.add_var_opt('trig-start-time', valid_seg[0])
        node.add_var_opt('trig-end-time', valid_seg[1])

        cache_file = dfParents[0]
        
        if self.injection_file is not None:
            node.add_input(self.injection_file, 'injection-file')

        # set the input and output files        
        node.make_and_add_output(valid_seg, '.xml.gz', 'output', tags=tags)
        node.add_input(cache_file, opt='frame-cache')
        node.add_input(parent, opt='bank-file')
        return node
        
    def get_valid_times(self):
        analysis_length = int(self.cp.get('ahope-inspiral', 'analysis-length'))
        pad_data = int(self.get_opt( 'pad-data'))
        start_pad = int(self.get_opt( 'segment-start-pad'))
        end_pad = int(self.get_opt('segment-end-pad'))

        data_length = analysis_length + pad_data * 2
        start = pad_data + start_pad
        end = data_length - pad_data - end_pad
        return data_length, segments.segment(start, end)

class PyCBCInspiralExec(Executable):
    """
    The class corresponding to the pycbc_inspiral executable. It can be used
    to create jobs and from that create nodes
    """
    def create_job(self, cp, ifo, out_dir=None, injection_file=None, tags=[]):
        return PyCBCInspiralJob(cp, self.exe_name, self.condor_universe,
                                ifo=ifo, out_dir=out_dir,
                                injection_file=injection_file, tags=tags)

class PyCBCTmpltbankJob(Job):
    """
    The class used to create jobs for pycbc_geom_nonspin_bank executable and
    any other executables using the same command line option groups.
    """
    def __init__(self, cp, exe_name, universe, ifo=None, out_dir=None,
                 tags=[]):
        Job.__init__(self, cp, exe_name, universe, ifo, out_dir, tags=tags)
        self.cp = cp
        self.set_memory(2000)

    def create_node(self, data_seg, valid_seg, parent=None, dfParents=None, tags=[]):
        node = Node(self)

        if not dfParents or len(dfParents) != 1:
            raise ValueError("%s must be supplied with a single cache file"
                              %(self.exe_name))

        pad_data = int(self.get_opt('pad-data'))
        if pad_data is None:
            raise ValueError("The option pad-data is a required option of "
                             "%s. Please check the ini file." % self.exe_name)

        # set the remaining option flags
        node.add_var_opt('gps-start-time', data_seg[0] + pad_data)
        node.add_var_opt('gps-end-time', data_seg[1] - pad_data)

        cache_file = dfParents[0]

        # set the input and output files      
        node.make_and_add_output(valid_seg, '.xml.gz', 'output-file', tags=tags)
        node.add_input(cache_file, opt='frame-cache')
        return node

    def create_nodata_node(self, valid_seg):
        """
        A simplified version of create_node that creates a node that does not
        need to read in data.
 
        Parameters
        -----------
        valid_seg : glue.segment
            The segment over which to declare the node valid. Usually this
            would be the duration of the analysis.

        Returns
        --------
        node : ahope.Node
            The instance corresponding to the created node.
        """
        node = Node(self)

        # Set the output file
        node.make_and_add_output(valid_seg, '.xml.gz', 'output-file')
        return node

    def get_valid_times(self):
        pad_data = int(self.get_opt( 'pad-data'))
        analysis_length = int(self.cp.get('ahope-tmpltbank', 'analysis-length'))
        
        #FIXME this should not be hard coded 
        data_length = analysis_length + pad_data * 2
        start = pad_data
        end = data_length - pad_data
        return data_length, segments.segment(start, end)

class PyCBCTmpltbankExec(Executable):
    """
    The class corresponding to pycbc_geom_nonspin_bank executable, and any
    other executables using the same command line option groups.
    """
    def create_job(self, cp, ifo, out_dir=None, tags=[]):
        return PyCBCTmpltbankJob(cp, self.exe_name, self.condor_universe,
                                 ifo=ifo, out_dir=out_dir, tags=tags)

class LigoLWCombineSegs(Job):
    """ 
    This class is used to create nodes for the ligolw_combine_segments 
    executable
    """
    def create_node(self, valid_seg, veto_files, segment_name):
        node = Node(self)
        node.add_var_opt('segment-name', segment_name)
        for fil in veto_files:
            node.add_input(fil, argument=True)   
        node.make_and_add_output(valid_seg, '.xml', 'output')      
        return node

class LigolwAddJob(Job):
    """
    The class used to create nodes for the ligolw_add executable.
    """
    def __init__(self, cp, exe_name, universe=None, ifo=None, out_dir=None, tags=[]):
        Job.__init__(self, cp, exe_name, universe, ifo, out_dir, tags=tags)
        self.set_memory(2000)

    def create_node(self, jobSegment, input_files, output=None):
        node = Node(self)

        # Very few options to ligolw_add, all input files are given as a long
        # argument list. If this becomes unwieldy we could dump all these files
        # to a cache file and read that in. ALL INPUT FILES MUST BE LISTED AS
        # INPUTS (with .add_input_file) IF THIS IS DONE THOUGH!
        for fil in input_files:
            node.add_input(fil, argument=True)

        # Currently we set the output file using the name of *all* active ifos,
        # even if one or more of these ifos is not active within jobSegment.
        # In ihope, these files were named with only participating ifos. Not
        # sure this is worth doing, can can be done with replacing self.ifo
        # here if desired
        if output:
            node.add_output(output, 'output')
        else:
            node.make_and_add_output(jobSegment, '.xml.gz', 'output')
        return node

class LigolwAddExec(Executable):
    """
    The class corresponding to the ligolw_add executable.
    """
    def __init__(self, exe_name):
        if exe_name != 'llwadd':
            raise ValueError('ligolw_add does not support setting '
                             'the exe_name to anything but "llwadd"')

        Executable.__init__(self, exe_name, 'vanilla')

    def create_job(self, cp, ifo, out_dir=None, tags=[]):
        return LigolwAddJob(cp, self.exe_name, self.condor_universe,
                            ifo=ifo, out_dir=out_dir, tags=tags)

class LigolwSSthincaJob(Job):
    """
    The class responsible for making jobs for ligolw_sstinca.
    """
    def __init__(self, cp, exe_name, universe, ifo=None, out_dir=None,
                 dqVetoName=None, tags=[]):
        Job.__init__(self, cp, exe_name, universe, ifo, out_dir, tags=tags)
        self.set_memory(2000)
        if dqVetoName:
            self.add_opt("vetoes-name", dqVetoName)

    def create_node(self, jobSegment, inputFile):
        node = Node(self)
        node.add_input(inputFile, argument=True)

        # Add the start/end times
        segString = "%f:%f" %(jobSegment[0], jobSegment[1]) 
        node.add_var_opt('coinc-end-time-segment', segString)

        # FIXME: This must match the *actual* output name!
        outFile = AhopeFile(self.ifo, self.exe_name, extension='.xml.gz',
                         segment=jobSegment,
                         directory=self.out_dir,
                         tags=self.tags)

        node.add_output(outFile)

        return node


class LigolwSSthincaExec(Executable):
    """
    The class corresponding to the ligolw_sstinca executable.
    """
    def __init__(self, exe_name):
        if exe_name != 'thinca':
            raise ValueError('ligolw_sstinca does not support setting '
                             'the exe_name to anything but "thinca"')

        Executable.__init__(self, exe_name, 'vanilla')

    def create_job(self, cp, ifo, out_dir=None, dqVetoName=None, tags=[]):
        return LigolwSSthincaJob(cp, self.exe_name, self.condor_universe,
                            ifo=ifo, out_dir=out_dir, 
                            dqVetoName=dqVetoName, tags=tags)
                            
class LalappsInspinjJob(Job):
    """
    The class used to create jobs for the lalapps_inspinj executable.
    """
    def create_node(self, segment):
        node = Node(self)
        
        if self.get_opt('write-compress') is not None:
            ext = '.xml.gz'
        else:
            ext = '.xml'
        
        node.add_var_opt('gps-start-time', segment[0])
        node.add_var_opt('gps-end-time', segment[1])    
        node.make_and_add_output(segment, '.xml', 'output')
        return node

class LalappsInspinjExec(Executable):
    """
    The class corresponding to the lalapps_inspinj executable.
    """
    def create_job(self, cp, out_dir=None, tags=[]):
        # FIXME: It is convention to name injection files with a 'HL' prefix
        # therefore I have hardcoded ifo=HL here. Maybe not a FIXME, but just
        # noting this.
        return LalappsInspinjJob(cp, self.exe_name, self.condor_universe,
                                 ifo='HL', out_dir=out_dir, tags=tags)
