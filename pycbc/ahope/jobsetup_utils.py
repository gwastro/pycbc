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
from glue import segments
from pycbc.ahope.ahope_utils import *
from pycbc.ahope.legacy_ihope import *

def select_tmpltbank_class(curr_exe):
    """
    This function returns an instance of the class that is appropriate for
    creating a template bank within ihope.
    
    Parameters
    ----------
    curr_exe : string
        The name of the AhopeExecutable that is being used.
    curr_section : string
        The name of the section storing options for this executble

    Returns
    --------
    Instanced class : exe_class
        An instance of the class that holds the utility functions appropriate
        for the given AhopeExecutable. This class **must** contain
        * exe_class.create_job()
        and the job returned by this **must** contain
        * job.get_valid_times(ifo, )
        * job.create_node()
    """
    # This is basically a list of if statements

    if curr_exe == 'lalapps_tmpltbank_ahope':
        exe_class = LegacyTmpltbankExecutable
    elif curr_exe == 'pycbc_geom_nonspinbank' or 'pycbc_aligned_stoch_bank':
        # These two codes have the same interface (shared option groups) and
        # can therefore used the same internal class
        exe_class = PyCBCTmpltbankExecutable
    elif curr_exe == 'pycbc_geom_nonspinbank':
        exe_class = PyCBCTmpltbankExecutable
    else:
        # Should we try some sort of default class??
        err_string = "No class exists for AhopeExecutable %s" %(curr_exe,)
        raise NotImplementedError(err_string)
    return exe_class

def select_matchedfilter_class(curr_exe):
    """
    This function returns an instance of the class that is appropriate for
    matched-filtering within ahope.
    
    Parameters
    ----------
    curr_exe : string
        The name of the AhopeExecutable that is being used.
    curr_section : string
        The name of the section storing options for this executble

    Returns
    --------
    Instanced class : exe_class
        An instance of the class that holds the utility functions appropriate
        for the given AhopeExecutable. This class **must** contain
        * exe_class.create_job()
        and the job returned by this **must** contain
        * job.get_valid_times(ifo, )
        * job.create_node()
    """
    # This is basically a list of if statements
    if curr_exe == 'lalapps_inspiral_ahope':
        exe_class = LegacyInspiralExecutable
    elif curr_exe == 'pycbc_inspiral':
        exe_class = PyCBCInspiralExecutable
    else:
        # Should we try some sort of default class??
        err_string = "No class exists for AhopeExecutable %s" %(curr_exe,)
        raise NotImplementedError(err_string)        
    return exe_class

def select_splitfilejob_instance(curr_exe, curr_section):
    """
    This function returns an instance of the class that is appropriate for
    splitting an output file up within ahope (for e.g. splitbank).
    
    Parameters
    ----------
    curr_exe : string
        The name of the AhopeExecutable that is being used.
    curr_section : string
        The name of the section storing options for this executble

    Returns
    --------
    Instanced class : exe_class
        An instance of the class that holds the utility functions appropriate
        for the given AhopeExecutable. This class **must** contain
        * exe_class.create_job()
        and the job returned by this **must** contain
        * job.create_node()
    """
    # This is basically a list of if statements
    if curr_exe == 'lalapps_splitbank':
        exe_class = LegacySplitBankExecutable(curr_section)
    # Some elif statements
    elif curr_exe == 'pycbc_splitbank':
        exe_class = PycbcSplitBankExecutable(curr_section)
    else:
        # Should we try some sort of default class??
        err_string = "No class exists for AhopeExecutable %s" %(curr_exe,)
        raise NotImplementedError(errString)

    return exe_class

def select_generic_executable(workflow, exe_tag):
    """
    This function returns an instance of the class that is appropriate for
    running the curr_exe. Curr_exe should not be a "specialized" job that fits
    into one of the select_XXX_instance functions above. IE. not a matched
    filter instance, or a template bank instance. Such specialized instances
    require extra setup. The only requirements here is that we can run
    create_job on the AhopeExecutable instance, and create_node on the resulting
    Job class.

    Parameters
    ----------
    workflow : ahope.Workflow instance
        The ahope workflow instance.

    exe_tag : string
        The name of the section storing options for this AhopeExecutable and the
        option giving the AhopeExecutable path in the [AhopeExecutables] section.


    Returns
    --------
    exe_class : ahope.AhopeExecutable sub-class
        The class that holds the utility functions appropriate
        for the given AhopeExecutable. This class this **must** contain
        * exe.create_node()
    """
    exe_path = workflow.cp.get("executables", exe_tag)
    exe_name = os.path.basename(exe_path)
    if exe_name == 'ligolw_add':
        exe_class = LigolwAddExecutable
    elif exe_name == 'ligolw_sstinca':
        exe_class = LigolwSSthincaExecutable
    elif exe_name == 'pycbc_sqlite_simplify':
        exe_class = PycbcSqliteSimplifyExecutable
    elif exe_name == 'ligolw_cbc_cluster_coincs':
        exe_class = SQLInOutExecutable
    elif exe_name == 'ligolw_dbinjfind':
        exe_class = SQLInOutExecutable
    elif exe_name == 'lalapps_inspinj':
        exe_class = LalappsInspinjExecutable
    elif exe_name == 'pycbc_timeslides':
        exe_class = PycbcTimeslidesExecutable
    elif exe_name == 'pycbc_compute_durations':
        exe_class = ComputeDurationsExecutable
    elif exe_name == 'pycbc_calculate_far':
        exe_class = SQLInOutExecutable
    else:
        # Should we try some sort of default class??
        err_string = "No class exists for AhopeExecutable %s" %(exe_name,)
        raise NotImplementedError(err_string)

    return exe_class

def sngl_ifo_job_setup(workflow, ifo, out_files, curr_exe_job, science_segs, 
                       datafind_outs, output_dir, parents=None, 
                       link_job_instance=None, allow_overlap=True,
                       compatibility_mode=True):
    """
    This function sets up a set of single ifo jobs. A basic overview of how this
    works is as follows:

    * (1) Identify the length of data that each job needs to read in, and what
      part of that data the job is valid for.
    * START LOOPING OVER SCIENCE SEGMENTS
    * (2) Identify how many jobs are needed (if any) to cover the given science
      segment and the time shift between jobs. If no jobs continue.
    * START LOPPING OVER JOBS
    * (3) Identify the time that the given job should produce valid output (ie.
      inspiral triggers) over.
    * (4) Identify the data range that the job will need to read in to produce
      the aforementioned valid output.
    * (5) Identify all parents/inputs of the job.
    * (6) Add the job to the workflow
    * END LOOPING OVER JOBS
    * END LOOPING OVER SCIENCE SEGMENTS

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
        Coordinate the valid times with another AhopeExecutable.
    allow_overlap : boolean (optional, kwarg, default = True)
        If this is set the times that jobs are valid for will be allowed to
        overlap. This may be desired for template banks which may have some
        overlap in the times they cover. This may not be desired for inspiral
        jobs, where you probably want triggers recorded by jobs to not overlap
        at all.
    compatibility_mode : boolean (optional,  kwarg, default = False)
        If given the jobs will be tiled in the same method as used in inspiral
        hipe. This requires that link_job_instance is also given. If not given
        ahope's methods are used.

    Returns
    --------
    out_files : ahope.AhopeFileList
        A list of the files that will be generated by this step in the
        ahope workflow.
    """

    if compatibility_mode and not link_job_instance:
        errMsg = "Compability mode requires a link_job_instance."
        raise ValueError(errMsg)

    cp = workflow.cp
    
    # Set up the condorJob class for the current AhopeExecutable
    data_length, valid_chunk = curr_exe_job.get_valid_times()
    
    # Begin by getting analysis start and end, and start and end of time
    # that the output file is valid for
    valid_length = abs(valid_chunk)

    data_chunk = segments.segment([0, data_length])
    job_tag = curr_exe_job.name.upper()
    
    ########### (1) ############
    # Get the times that can be analysed and needed data lengths
    data_length, valid_chunk, valid_length = identify_needed_data(curr_exe_job,\
                                           link_job_instance=link_job_instance)

    # DO NOT! use valid length here, as valid_length and abs(valid_chunk)
    # may be different if using link_exe_instance
    data_loss = data_length - abs(valid_chunk)

    
    # Loop over science segments and set up jobs
    for curr_seg in science_segs:
        ########### (2) ############
        # Initialize the class that identifies how many jobs are needed and the
        # shift between them.
        segmenter = JobSegmenter(data_length, valid_chunk, valid_length, 
                    curr_seg, data_loss, compatibility_mode=compatibility_mode)

        for job_num in range(segmenter.num_jobs):
            ############## (3) #############
            # Figure out over what times this job will be valid for

            job_valid_seg = segmenter.get_valid_times_for_job(job_num,
                                                   allow_overlap=allow_overlap)

            ############## (4) #############
            # Get the data that this job should read in

            job_data_seg = segmenter.get_data_times_for_job(job_num)
                
            ############# (5) ############
            # Identify parents/inputs to the job

            if parents:
                # Find the set of files with the best overlap
                curr_parent = parents.find_outputs_in_range(ifo, job_valid_seg,
                                                            useSplitLists=True)
                if not curr_parent:
                    err_string = ("No parent jobs found overlapping %d to %d." 
                                  %(job_valid_seg[0], job_valid_seg[1]))
                    err_string += "\nThis is a bad error! Contact a developer."
                    raise ValueError(err_string)
            else:
                curr_parent = [None]


            if datafind_outs:
                curr_dfouts = datafind_outs.find_all_output_in_range(ifo, 
                                              job_data_seg, useSplitLists=True)
                if not curr_dfouts:
                    err_str = ("No datafind jobs found overlapping %d to %d."
                                %(job_data_seg[0],job_data_seg[1]))
                    err_str += "\nThis shouldn't happen. Contact a developer."
                    raise ValueError(err_str)


            ############## (6) #############
            # Make node and add to workflow

            # Note if I have more than one curr_parent I need to make more than
            # one job. If there are no curr_parents it is set to [None] and I
            # make a single job. This catches the case of a split template bank
            # where I run a number of jobs to cover a single range of time.
            for pnum, parent in enumerate(curr_parent):
                if len(curr_parent) != 1:
                    tag = ["JOB%d" %(pnum,)]
                else:
                    tag = []
                # To ensure output file uniqueness I add a tag
                # We should generate unique names automatically, but it is a 
                # pain until we can set the output names for all AhopeExecutables              
                node = curr_exe_job.create_node(job_data_seg, job_valid_seg, 
                                                parent=parent,
                                                dfParents=curr_dfouts, 
                                                tags=tag)
                workflow.add_node(node)
                curr_out_files = node.output_files
                # FIXME: Here we remove PSD files if they are coming through.
                #        This should be done in a better way. On to-do list.
                curr_out_files = [i for i in curr_out_files if 'PSD_FILE'\
                                                                 not in i.tags]
                out_files += curr_out_files

    return out_files

def identify_needed_data(curr_exe_job, link_job_instance=None):
    """
    This function will identify the length of data that a specific executable
    needs to analyse and what part of that data is valid (ie. inspiral doesn't
    analyse the first or last 64+8s of data it reads in).
    In addition you can supply a second job instance to "link" to, which will
    ensure that the two jobs will have a one-to-one correspondence (ie. one
    template bank per one matched-filter job) and the corresponding jobs will
    be "valid" at the same times.

    Parameters
    -----------
    curr_exe_job : ahope.Job
        An instanced of the PyCBC Job class that has a get_valid times method.
    link_job_instance : Job instance (optional),
        Coordinate the valid times with another executable.

    Returns
    --------
    dataLength : float
        The amount of data (in seconds) that each instance of the job must read
        in.
    valid_chunk : glue.segment.segment
        The times within dataLength for which that jobs output **can** be
        valid (ie. for inspiral this is (72, dataLength-72) as, for a standard
        setup the inspiral job cannot look for triggers in the first 72 or
        last 72 seconds of data read in.)
    valid_length : float
        The maximum length of data each job can be valid for. If not using
        link_job_instance this is abs(valid_segment), but can be smaller than
        that if the linked job only analyses a small amount of data (for e.g.).
    """
    # Set up the condorJob class for the current executable
    data_length, valid_chunk = curr_exe_job.get_valid_times()

    # Begin by getting analysis start and end, and start and end of time
    # that the output file is valid for
    valid_length = abs(valid_chunk)

    data_chunk = segments.segment([0, data_length])
    job_tag = curr_exe_job.name.upper()

    if link_job_instance:
        # FIXME: Should we remove this, after testing is complete??
        # EURGHH! What we are trying to do here is, if this option is given,
        # line up the template bank and inspiral jobs so that there is one
        # template bank for each inspiral job. This is proving a little messy
        # and probably still isn't perfect.

        # What data does the linked exe use?
        link_data_length,link_valid_chunk = link_job_instance.get_valid_times()
        # What data is lost at start of both jobs? Take the maximum.
        start_data_loss = max(valid_chunk[0], link_valid_chunk[0])
        # What data is lost at end of both jobs? Take the maximum.
        end_data_loss = max(data_length - valid_chunk[1],\
                            link_data_length - link_valid_chunk[1])
        # Calculate valid_segments for both jobs based on the combined data
        # loss.
        valid_chunk = segments.segment(start_data_loss, \
                                       data_length - end_data_loss)
        link_valid_chunk = segments.segment(start_data_loss, \
                                       link_data_length - end_data_loss)

        # The maximum valid length should be the minimum of the two
        link_valid_length = abs(link_valid_chunk)

        # Which one is now longer? Use this is valid_length
        if link_valid_length < valid_length:
            valid_length = link_valid_length

    # DO NOT! use valid length here, as valid_length and abs(valid_chunk)
    # may be different by this point if using link_exe_instance
    data_loss = data_length - abs(valid_chunk)


    if data_loss < 0:
        raise ValueError("Ahope needs fixing! Please contact a developer")

    return data_length, valid_chunk, valid_length


class JobSegmenter(object):
    """
    This class is used when running sngl_ifo_job_setup to determine what times
    should be analysed be each job and what data is needed.
    """
    def __init__(self, data_length, valid_chunk, valid_length, curr_seg,
                 data_loss, compatibility_mode = False):
        """
        Initialize class.
        """
        self.curr_seg = curr_seg
        self.curr_seg_length = float(abs(curr_seg))
        self.data_loss = data_loss
        self.valid_chunk = valid_chunk
        self.valid_length = valid_length
        self.data_length = data_length
        self.data_chunk = segments.segment([0, self.data_length])

        if self.curr_seg_length < data_length:
            self.num_jobs = 0
            return
        # How many jobs do we need
        self.num_jobs = int( math.ceil( (self.curr_seg_length \
                                - self.data_loss) / float(self.valid_length) ))

        self.compatibility_mode = compatibility_mode
        if compatibility_mode and (self.valid_length != abs(self.valid_chunk)):
            errMsg = "In compatibility mode the template bank and matched-"
            errMsg += "filter jobs must read in the same amount of data."
            raise ValueError(errMsg)
        elif compatibility_mode:
            # What is the incremental shift between jobs
            self.job_time_shift = self.valid_length
        else:
            # What is the incremental shift between jobs
            self.job_time_shift = (self.curr_seg_length - self.data_length) / \
                                   float(self.num_jobs - 1)


    def get_valid_times_for_job(self, num_job, allow_overlap=True):
        """
        Get the times for which this job is valid.
        """
        if self.compatibility_mode:
            return self.get_valid_times_for_job_legacy(num_job)
        else:
           return self.get_valid_times_for_job_ahope(num_job, 
                                                   allow_overlap=allow_overlap)

    def get_valid_times_for_job_ahope(self, num_job, allow_overlap=True):
        """
        Get the times for which the job num_job will be valid, using ahope's
        method.
        """
        # small factor of 0.0001 to avoid float round offs causing us to
        # miss a second at end of segments.
        shift_dur = self.curr_seg[0] + int(self.job_time_shift * num_job\
                                           + 0.0001)
        job_valid_seg = self.valid_chunk.shift(shift_dur)
        # If we need to recalculate the valid times to avoid overlap
        if not allow_overlap:
            data_per_job = (self.curr_seg_length - self.data_loss) / \
                           float(self.num_jobs)
            lower_boundary = num_job*data_per_job + \
                                 self.valid_chunk[0] + self.curr_seg[0]
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
        return job_valid_seg

    def get_valid_times_for_job_legacy(self, num_job):
        """
        Get the times for which the job num_job will be valid, using the method
        use in inspiral hipe.
        """
        # All of this should be integers, so no rounding factors needed.
        shift_dur = self.curr_seg[0] + int(self.job_time_shift * num_job)
        job_valid_seg = self.valid_chunk.shift(shift_dur)

        # If this is the last job, push the end back
        if num_job == (self.num_jobs - 1):
            dataPushBack = self.data_length - self.valid_chunk[1]
            job_valid_seg = segments.segment(job_valid_seg[0],
                                               self.curr_seg[1] - dataPushBack)

        return job_valid_seg

    def get_data_times_for_job(self, num_job):
        """
        Get the data that this job will read in.
        """
        if self.compatibility_mode:
            return self.get_data_times_for_job_legacy(num_job)
        else:
           return self.get_data_times_for_job_ahope(num_job)


    def get_data_times_for_job_ahope(self, num_job):
        """
        Get the data that this job will need to read in.
        """
        # small factor of 0.0001 to avoid float round offs causing us to
        # miss a second at end of segments.
        shift_dur = self.curr_seg[0] + int(self.job_time_shift * num_job\
                                           + 0.0001)
        job_data_seg = self.data_chunk.shift(shift_dur)
        # Sanity check that all data is used
        if num_job == 0:
           if job_data_seg[0] != self.curr_seg[0]:
               errMsg = "Job is not using data from the start of the "
               errMsg += "science segment. It should be using all data."
               raise ValueError(errMsg)
        if num_job == (self.num_jobs - 1):
            if job_data_seg[1] != self.curr_seg[1]:
                errMsg = "Job is not using data from the end of the "
                errMsg += "science segment. It should be using all data."
                raise ValueError(errMsg)

        return job_data_seg

    def get_data_times_for_job_legacy(self, num_job):
        """
        Get the data that this job will need to read in.
        """
        # Should all be integers, so no rounding needed
        shift_dur = self.curr_seg[0] + int(self.job_time_shift * num_job)
        job_data_seg = self.data_chunk.shift(shift_dur)

        # If this is the last job, push the end back
        if num_job == (self.num_jobs - 1):
            dataPushBack = job_data_seg[1] - self.curr_seg[1]
            assert dataPushBack >= 0
            job_data_seg = segments.segment(job_data_seg[0] - dataPushBack,
                                                 self.curr_seg[1])
            assert (abs(job_data_seg) == self.data_length)

        # Sanity check that all data is used
        if num_job == 0:
           if job_data_seg[0] != self.curr_seg[0]:
               errMsg = "Job is not using data from the start of the "
               errMsg += "science segment. It should be using all data."
               raise ValueError(errMsg)
        if num_job == (self.num_jobs - 1):
            if job_data_seg[1] != self.curr_seg[1]:
                errMsg = "Job is not using data from the end of the "
                errMsg += "science segment. It should be using all data."
                raise ValueError(errMsg)

        return job_data_seg

class PyCBCInspiralExecutable(AhopeExecutable):
    """
    The class used to create jobs for pycbc_inspiral AhopeExecutable.
    """
    def __init__(self, cp, exe_name, ifo=None, out_dir=None, injection_file=None, tags=[]):
        AhopeExecutable.__init__(self, cp, exe_name, None, ifo, out_dir, tags=tags)
        self.cp = cp
        self.set_memory(2000)
        self.injection_file = injection_file

        if self.get_opt('processing-scheme') == 'cuda':
            self.needs_gpu()
         
        self.num_threads = 1  
        if self.get_opt('processing-scheme') is not None:
            stxt = self.get_opt('processing-scheme')
            if len(stxt.split(':')) > 1:
                self.num_threads = stxt.split(':')[1]

    def create_node(self, data_seg, valid_seg, parent=None, dfParents=None, tags=[]):
        node = AhopeNode(self)
        pad_data = int(self.get_opt('pad-data'))
        if pad_data is None:
            raise ValueError("The option pad-data is a required option of "
                             "%s. Please check the ini file." % self.name)

        if not dfParents:
            raise ValueError("%s must be supplied with data file(s)"
                              %(self.name))

        # set remaining options flags   
        node.add_opt('--gps-start-time', data_seg[0] + pad_data)
        node.add_opt('--gps-end-time', data_seg[1] - pad_data)
        node.add_opt('--trig-start-time', valid_seg[0])
        node.add_opt('--trig-end-time', valid_seg[1])

        node.add_profile('condor', 'request_cpus', self.num_threads)        

        if self.injection_file is not None:
            node.add_input_opt('--injection-file', self.injection_file)

        # set the input and output files        
        node.new_output_file_opt(valid_seg, '.xml.gz', '--output', tags=tags)
        node.add_input_list_opt('--frame-files', dfParents)
        node.add_input_opt('--bank-file', parent, )

        # FIXME: This hack is needed for pipedown compatibility. user-tag is
        #        no-op and is only needed for pipedown to know whether this is
        #        a "FULL_DATA" job or otherwise.
        outFile = os.path.basename(node.output_files[0].storage_path)
        userTag = outFile.split('-')[1]
        userTag = userTag.split('_')[1:]
        if userTag[0] == 'FULL' and userTag[1] == 'DATA':
            userTag = 'FULL_DATA'
        elif userTag[0] == 'PLAYGROUND':
            userTag = 'PLAYGROUND'
        elif userTag[0].endswith("INJ"):
            userTag = userTag[0]
        else:
            userTag = '_'.join(userTag)
        node.add_opt("--user-tag", userTag)

        return node
        
    def get_valid_times(self):
        # FIXME: IAN. I'm not happy about analysis_length being buried here.
        #        Maybe this should be something read in at the
        #        matchedfilter_utils level, and acted on *if* possible.
        analysis_length = int(self.cp.get('ahope-matchedfilter',
                                          'analysis-length'))
        pad_data = int(self.get_opt( 'pad-data'))
        start_pad = int(self.get_opt( 'segment-start-pad'))
        end_pad = int(self.get_opt('segment-end-pad'))

        data_length = analysis_length + pad_data * 2
        start = pad_data + start_pad
        end = data_length - pad_data - end_pad
        return data_length, segments.segment(start, end)

class PyCBCTmpltbankExecutable(AhopeExecutable):
    """
    The class used to create jobs for pycbc_geom_nonspin_bank AhopeExecutable and
    any other AhopeExecutables using the same command line option groups.
    """
    def __init__(self, cp, exe_name, ifo=None, out_dir=None,
                 tags=[], write_psd=False):
        AhopeExecutable.__init__(self, cp, exe_name, 'vanilla', ifo, out_dir, tags=tags)
        self.cp = cp
        self.set_memory(2000)
        self.write_psd = write_psd

    def create_node(self, data_seg, valid_seg, parent=None, dfParents=None, tags=[]):
        node = AhopeNode(self)

        if not dfParents:
            raise ValueError("%s must be supplied with data file(s)"
                              %(self.name))

        pad_data = int(self.get_opt('pad-data'))
        if pad_data is None:
            raise ValueError("The option pad-data is a required option of "
                             "%s. Please check the ini file." % self.name)

        # set the remaining option flags
        node.add_opt('--gps-start-time', data_seg[0] + pad_data)
        node.add_opt('--gps-end-time', data_seg[1] - pad_data)

        # set the input and output files      
        # Add the PSD file if needed
        if self.write_psd:
            node.make_and_add_output(valid_seg, 'txt', '--psd-output',
                                     tags=tags+['PSD_FILE'])
        node.new_output_file_opt(valid_seg, '.xml.gz', '--output-file', tags=tags)
        node.add_input_list_opt('--frame-files', dfParents)
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
        node = AhopeNode(self)

        # Set the output file
        # Add the PSD file if needed
        if self.write_psd:
            node.make_and_add_output(valid_seg, 'txt', 'psd-output', 
                                     tags=tags+['PSD_FILE'])

        node.new_output_file_opt(valid_seg, '.xml.gz', '--output-file')
        return node

    def get_valid_times(self):
        pad_data = int(self.get_opt( 'pad-data'))
        analysis_length = int(self.cp.get('ahope-tmpltbank', 'analysis-length'))
        
        #FIXME this should not be hard coded 
        data_length = analysis_length + pad_data * 2
        start = pad_data
        end = data_length - pad_data
        return data_length, segments.segment(start, end)

class LigoLWCombineSegsExecutable(AhopeExecutable):
    """ 
    This class is used to create nodes for the ligolw_combine_segments 
    AhopeExecutable
    """
    def create_node(self, valid_seg, veto_files, segment_name):
        node = AhopeNode(self)
        node.add_opt('--segment-name', segment_name)
        for fil in veto_files:
            node.add_input_arg(fil)   
        node.new_output_file_opt(valid_seg, '.xml', '--output')      
        return node

class LigolwAddExecutable(AhopeExecutable):
    """
    The class used to create nodes for the ligolw_add AhopeExecutable.
    """
    def __init__(self, cp, exe_name, universe=None, ifo=None, out_dir=None, tags=[]):
        AhopeExecutable.__init__(self, cp, exe_name, universe, ifo, out_dir, tags=tags)
        self.set_memory(2000)

    def create_node(self, jobSegment, input_files, output=None, tags=[]):
        node = AhopeNode(self)

        # Very few options to ligolw_add, all input files are given as a long
        # argument list. If this becomes unwieldy we could dump all these files
        # to a cache file and read that in. ALL INPUT FILES MUST BE LISTED AS
        # INPUTS (with .add_input_opt_file) IF THIS IS DONE THOUGH!
        for fil in input_files:
            node.add_input_arg(fil)

        # Currently we set the output file using the name of *all* active ifos,
        # even if one or more of these ifos is not active within jobSegment.
        # In ihope, these files were named with only participating ifos. Not
        # sure this is worth doing, can can be done with replacing self.ifo
        # here if desired
        if output:
            node.add_output_opt('--output', output)
        else:
            node.new_output_file_opt(jobSegment, '.xml.gz', '--output', tags=tags)
        return node

class LigolwSSthincaExecutable(AhopeExecutable):
    """
    The class responsible for making jobs for ligolw_sstinca.
    """
    def __init__(self, cp, exe_name, universe=None, ifo=None, out_dir=None,
                 dqVetoName=None, tags=[]):
        AhopeExecutable.__init__(self, cp, exe_name, universe, ifo, out_dir, tags=tags)
        self.set_memory(2000)
        if dqVetoName:
            self.add_opt("--vetoes-name", dqVetoName)

    def create_node(self, jobSegment, coincSegment, inputFile, tags=[]):
        node = AhopeNode(self)
        node.add_input_arg(inputFile)

        # Add the start/end times
        segString = ""
        if coincSegment[0]:
          segString += str(coincSegment[0])
        segString += ":"
        if coincSegment[1]:
          segString += str(coincSegment[1])

        node.add_opt('--coinc-end-time-segment', segString)

        # FIXME: This must match the *actual* output name!
        outFile = AhopeFile(self.ifo, self.name, jobSegment,
                         extension='.xml.gz', directory=self.out_dir,
                         tags=self.tags+tags)

        node._add_output(outFile)

        return node

class PycbcSqliteSimplifyExecutable(AhopeExecutable):
    """
    The class responsible for making jobs for pycbc_sqlite_simplify.
    """
    def __init__(self, cp, exe_name, universe=None, ifo=None, out_dir=None, tags=[]):
        AhopeExecutable.__init__(self, cp, exe_name, universe, ifo, out_dir, tags=tags)
        self.set_memory(2000)
        
    def create_node(self, job_segment, inputFiles, injFile=None, injString=None):
        node = AhopeNode(self)
        if injFile and not injString:
            raise ValueError("injString needed if injFile supplied.")
        for file in inputFiles:
            node.add_input_arg(file)
        if injFile:
            node.add_input_opt("--injection-file", injFile)
            node.add_opt("--simulation-tag", injString)
        node.new_output_file_opt(job_segment, '.sql', '--output-file',
                                 tags=self.tags) 
        return node

class SQLInOutExecutable(AhopeExecutable):
    """
    The class responsible for making jobs for SQL codes taking one input and
    one output.
    """
    def __init__(self, cp, exe_name, universe=None, ifo=None, out_dir=None, tags=[]):
        AhopeExecutable.__init__(self, cp, exe_name, universe, ifo, out_dir, tags=tags)

    def create_node(self, job_segment, inputFile):
        node = AhopeNode(self)
        node.add_input_opt('--input', inputFile)
        node.new_output_file_opt(job_segment, '.sql', '--output',
                                 tags=self.tags)
        return node

class ComputeDurationsExecutable(SQLInOutExecutable):
    """
    The class responsible for making jobs for pycbc_compute_durations.
    """
    def create_node(self, job_segment, input_file, summary_xml_file):
        node = AhopeNode(self)
        node.add_input_opt('--input', input_file)
        node.add_input_opt('--segment-file', summary_xml_file)
        node.new_output_file_opt(job_segment, '.sql', '--output',
                                 tags=self.tags)
        return node

class LalappsInspinjExecutable(AhopeExecutable):
    """
    The class used to create jobs for the lalapps_inspinj AhopeExecutable.
    """
    def create_node(self, segment):
        node = AhopeNode(self)
        
        if self.get_opt('write-compress') is not None:
            ext = '.xml.gz'
        else:
            ext = '.xml'
        
        node.add_opt('--gps-start-time', segment[0])
        node.add_opt('--gps-end-time', segment[1])    
        node.new_output_file_opt(segment, '.xml', '--output')
        return node

class PycbcTimeslidesExecutable(AhopeExecutable):
    """
    The class used to create jobs for the pycbc_timeslides AhopeExecutable.
    """
    def create_node(self, segment):
        node = AhopeNode(self)

        node.new_output_file_opt(segment, '.xml.gz', '--output-files')
        return node

class PycbcSplitBankExecutable(AhopeExecutable):
    """
    The class responsible for creating jobs for pycbc_splitbank.
    """
    def __init__(self, cp, exe_name, num_banks,
                 ifo=None, out_dir=None, tags=[], universe=None):
        AhopeExecutable.__init__(self, cp, exe_name, universe, ifo, out_dir, tags=tags)
        self.num_banks = int(num_banks)

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
        node = AhopeNode(self)
        node.add_input_opt('--bank-file', bank)

        # Get the output (taken from inspiral.py)
        out_files = AhopeFileList([])
        for i in range( 0, self.num_banks):
            curr_tag = 'bank%d' %(i)
            # FIXME: What should the tags actually be? The job.tags values are
            #        currently ignored.
            curr_tags = bank.tags + [curr_tag]
            job_tag = bank.description + "_" + self.name.upper()
            out_file = AhopeFile(bank.ifoList, job_tag, bank.segment,
                                 extension=".xml.gz", directory=self.out_dir,
                                 tags=curr_tags)
            out_files.append(out_file)
        node.add_output_list_opt('--output-filenames', out_files)
        return node
