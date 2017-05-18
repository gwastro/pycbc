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
and add jobs/nodes to a pycbc workflow. For details about pycbc.workflow see:
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/ahope.html
"""

import logging
import math, os
import lal
from pycbc_glue import segments
import Pegasus.DAX3 as dax
from pycbc.workflow.core import Executable, File, FileList, Node

def int_gps_time_to_str(t):
    """Takes an integer GPS time, either given as int or lal.LIGOTimeGPS, and
    converts it to a string. If a LIGOTimeGPS with nonzero decimal part is
    given, raises a ValueError."""

    if type(t) == int:
        return str(t)
    elif type(t) == lal.LIGOTimeGPS:
        if t.gpsNanoSeconds == 0:
            return str(t.gpsSeconds)
        else:
            raise ValueError('Need an integer GPS time, got %s' % str(t))

def select_tmpltbank_class(curr_exe):
    """ This function returns a class that is appropriate for setting up
    template bank jobs within workflow.

    Parameters
    ----------
    curr_exe : string
        The name of the executable to be used for generating template banks.

    Returns
    --------
    exe_class : Sub-class of pycbc.workflow.core.Executable that holds utility
        functions appropriate for the given executable.  Instances of the class
        ('jobs') **must** have methods
        * job.create_node()
        and
        * job.get_valid_times(ifo, )
    """
    exe_to_class_map = {
        'pycbc_geom_nonspinbank'  : PyCBCTmpltbankExecutable,
        'pycbc_aligned_stoch_bank': PyCBCTmpltbankExecutable
    }
    try:
        return exe_to_class_map[curr_exe]
    except KeyError:
        raise NotImplementedError(
            "No job class exists for executable %s, exiting" % curr_exe)

def select_matchedfilter_class(curr_exe):
    """ This function returns a class that is appropriate for setting up
    matched-filtering jobs within workflow.

    Parameters
    ----------
    curr_exe : string
        The name of the matched filter executable to be used.

    Returns
    --------
    exe_class : Sub-class of pycbc.workflow.core.Executable that holds utility
        functions appropriate for the given executable.  Instances of the class
        ('jobs') **must** have methods
        * job.create_node()
        and
        * job.get_valid_times(ifo, )
    """
    exe_to_class_map = {
        'pycbc_inspiral'          : PyCBCInspiralExecutable,
        'pycbc_inspiral_skymax'   : PyCBCInspiralExecutable
    }
    try:
        return exe_to_class_map[curr_exe]
    except KeyError:
        # also conceivable to introduce a default class??
        raise NotImplementedError(
            "No job class exists for executable %s, exiting" % curr_exe)

def select_generic_executable(workflow, exe_tag):
    """ Returns a class that is appropriate for setting up jobs to run executables
    having specific tags in the workflow config.
    Executables should not be "specialized" jobs fitting into one of the
    select_XXX_class functions above, i.e. not a matched filter or template
    bank job, which require extra setup.

    Parameters
    ----------
    workflow : pycbc.workflow.core.Workflow
        The Workflow instance.

    exe_tag : string
        The name of the config section storing options for this executable and
        the option giving the executable path in the [executables] section.

    Returns
    --------
    exe_class : Sub-class of pycbc.workflow.core.Executable that holds utility
        functions appropriate for the given executable.  Instances of the class
        ('jobs') **must** have a method job.create_node()
    """
    exe_path = workflow.cp.get("executables", exe_tag)
    exe_name = os.path.basename(exe_path)
    exe_to_class_map = {
        'ligolw_add'               : LigolwAddExecutable,
        'ligolw_cbc_sstinca'           : LigolwSSthincaExecutable,
        'pycbc_sqlite_simplify'    : PycbcSqliteSimplifyExecutable,
        'ligolw_cbc_cluster_coincs': SQLInOutExecutable,
        'ligolw_cbc_repop_coinc'   : SQLInOutExecutable,
        'repop_coinc_expfit'       : SQLInOutExecutable,
        'ligolw_cbc_dbinjfind'         : SQLInOutExecutable,
        'lalapps_inspinj'          : LalappsInspinjExecutable,
        'pycbc_dark_vs_bright_injections' : PycbcDarkVsBrightInjectionsExecutable,
        'pycbc_timeslides'         : PycbcTimeslidesExecutable,
        'pycbc_compute_durations'  : ComputeDurationsExecutable,
        'pycbc_calculate_far'      : PycbcCalculateFarExecutable,
        "pycbc_run_sqlite"         : SQLInOutExecutable,
        # FIXME: We may end up with more than one class for using ligolw_sqlite
        #        How to deal with this?
        "ligolw_sqlite"            : ExtractToXMLExecutable,
        "pycbc_inspinjfind"        : InspinjfindExecutable,
        "pycbc_pickle_horizon_distances" : PycbcPickleHorizonDistsExecutable,
        "pycbc_combine_likelihood" : PycbcCombineLikelihoodExecutable,
        "pycbc_gen_ranking_data"   : PycbcGenerateRankingDataExecutable,
        "pycbc_calculate_likelihood"     : PycbcCalculateLikelihoodExecutable,
        "gstlal_inspiral_marginalize_likelihood"      : GstlalMarginalizeLikelihoodExecutable,
        "pycbc_compute_far_from_snr_chisq_histograms" : GstlalFarfromsnrchisqhistExecutable,
        "gstlal_inspiral_plot_sensitivity"            : GstlalPlotSensitivity,
        "gstlal_inspiral_plot_background" : GstlalPlotBackground,
        "gstlal_inspiral_plotsummary"     : GstlalPlotSummary,
        "gstlal_inspiral_summary_page"    : GstlalSummaryPage,
        "pycbc_condition_strain"         : PycbcConditionStrainExecutable
    }
    try:
        return exe_to_class_map[exe_name]
    except KeyError:
        # Should we try some sort of default class??
        raise NotImplementedError(
            "No job class exists for executable %s, exiting" % exe_name)

def sngl_ifo_job_setup(workflow, ifo, out_files, curr_exe_job, science_segs,
                       datafind_outs, parents=None,
                       link_job_instance=None, allow_overlap=True,
                       compatibility_mode=True):
    """ This function sets up a set of single ifo jobs. A basic overview of how this
    works is as follows:

    * (1) Identify the length of data that each job needs to read in, and what
      part of that data the job is valid for.
    * START LOOPING OVER SCIENCE SEGMENTS
    * (2) Identify how many jobs are needed (if any) to cover the given science
      segment and the time shift between jobs. If no jobs continue.
    * START LOOPING OVER JOBS
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
    workflow: pycbc.workflow.core.Workflow
        An instance of the Workflow class that manages the constructed workflow.
    ifo : string
        The name of the ifo to set up the jobs for
    out_files : pycbc.workflow.core.FileList
        The FileList containing the list of jobs. Jobs will be appended
        to this list, and it does not need to be empty when supplied.
    curr_exe_job : Job
        An instanced of the Job class that has a get_valid times method.
    science_segs : glue.segments.segmentlist
        The list of times that the jobs should cover
    datafind_outs : pycbc.workflow.core.FileList
        The file list containing the datafind files.
    parents : pycbc.workflow.core.FileList (optional, kwarg, default=None)
        The FileList containing the list of jobs that are parents to
        the one being set up.
    link_job_instance : Job instance (optional),
        Coordinate the valid times with another Executable.
    allow_overlap : boolean (optional, kwarg, default = True)
        If this is set the times that jobs are valid for will be allowed to
        overlap. This may be desired for template banks which may have some
        overlap in the times they cover. This may not be desired for inspiral
        jobs, where you probably want triggers recorded by jobs to not overlap
        at all.
    compatibility_mode : boolean (optional,  kwarg, default = False)
        If given the jobs will be tiled in the same method as used in inspiral
        hipe. This requires that link_job_instance is also given. If not given
        workflow's methods are used.

    Returns
    --------
    out_files : pycbc.workflow.core.FileList
        A list of the files that will be generated by this step in the
        workflow.
    """
    if compatibility_mode and not link_job_instance:
        errMsg = "Compability mode requires a link_job_instance."
        raise ValueError(errMsg)

    ########### (1) ############
    # Get the times that can be analysed and needed data lengths
    data_length, valid_chunk, valid_length = identify_needed_data(curr_exe_job,
                                           link_job_instance=link_job_instance)

    # Loop over science segments and set up jobs
    for curr_seg in science_segs:
        ########### (2) ############
        # Initialize the class that identifies how many jobs are needed and the
        # shift between them.
        segmenter = JobSegmenter(data_length, valid_chunk, valid_length,
                                 curr_seg, curr_exe_job,
                                 compatibility_mode=compatibility_mode)

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
                # pain until we can set the output names for all Executables
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

def multi_ifo_coherent_job_setup(workflow, out_files, curr_exe_job,
                                 science_segs, datafind_outs, output_dir,
                                 parents=None, tags=None):
    """
    Method for setting up coherent inspiral jobs.
    """
    if tags is None:
        tags = []
    data_seg, job_valid_seg = curr_exe_job.get_valid_times()
    curr_out_files = FileList([])
    if 'IPN' in datafind_outs[-1].description \
            and 'bank_veto_bank' in datafind_outs[-2].description:
        ipn_sky_points = datafind_outs[-1]
        bank_veto = datafind_outs[-2]
        frame_files = datafind_outs[:-2]
    else:
        ipn_sky_points = None
        bank_veto = datafind_outs[-1]
        frame_files = datafind_outs[:-1]

    split_bank_counter = 0

    if curr_exe_job.injection_file is None:
        for split_bank in parents:
            tag = list(tags)
            tag.append(split_bank.tag_str)
            node = curr_exe_job.create_node(data_seg, job_valid_seg,
                    parent=split_bank, dfParents=frame_files,
                    bankVetoBank=bank_veto, ipn_file=ipn_sky_points, tags=tag)
            workflow.add_node(node)
            split_bank_counter += 1
            curr_out_files.extend(node.output_files)
    else:
        for inj_file in curr_exe_job.injection_file:
            for split_bank in parents:
                tag = list(tags)
                tag.append(inj_file.tag_str)
                tag.append(split_bank.tag_str)
                node = curr_exe_job.create_node(data_seg, job_valid_seg,
                        parent=split_bank, inj_file=inj_file, tags=tag,
                        dfParents=frame_files, bankVetoBank=bank_veto,
                        ipn_file=ipn_sky_points)
                workflow.add_node(node)
                split_bank_counter += 1
                curr_out_files.extend(node.output_files)

    # FIXME: Here we remove PSD files if they are coming
    #        through. This should be done in a better way. On
    #        to-do list.
    curr_out_files = [i for i in curr_out_files if 'PSD_FILE'\
                      not in i.tags]
    out_files += curr_out_files

    return out_files

def identify_needed_data(curr_exe_job, link_job_instance=None):
    """ This function will identify the length of data that a specific executable
    needs to analyse and what part of that data is valid (ie. inspiral doesn't
    analyse the first or last 64+8s of data it reads in).

    In addition you can supply a second job instance to "link" to, which will
    ensure that the two jobs will have a one-to-one correspondence (ie. one
    template bank per one matched-filter job) and the corresponding jobs will
    be "valid" at the same times.

    Parameters
    -----------
    curr_exe_job : Job
        An instance of the Job class that has a get_valid times method.
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
    data_lengths, valid_chunks = curr_exe_job.get_valid_times()

    # Begin by getting analysis start and end, and start and end of time
    # that the output file is valid for
    valid_lengths = [abs(valid_chunk) for valid_chunk in valid_chunks]

    if link_job_instance:
        # FIXME: Should we remove this, after testing is complete??
        # EURGHH! What we are trying to do here is, if this option is given,
        # line up the template bank and inspiral jobs so that there is one
        # template bank for each inspiral job. This is proving a little messy
        # and probably still isn't perfect.

        # What data does the linked exe use?
        link_data_length, link_valid_chunk = link_job_instance.get_valid_times()
        if len(link_data_length) > 1 or len(valid_lengths) > 1:
            raise ValueError('Linking job instances for tiling is not supported'
                             ' between jobs that allow variable tile size')

        # What data is lost at start of both jobs? Take the maximum.
        start_data_loss = max(valid_chunks[0][0], link_valid_chunk[0][0])
        # What data is lost at end of both jobs? Take the maximum.
        end_data_loss = max(data_lengths[0] - valid_chunks[0][1],\
                            link_data_length[0] - link_valid_chunk[0][1])
        # Calculate valid_segments for both jobs based on the combined data
        # loss.

        valid_chunks[0] = segments.segment(start_data_loss, \
                                       data_lengths[0] - end_data_loss)
        link_valid_chunk = segments.segment(start_data_loss, \
                                       link_data_length[0] - end_data_loss)

        # The maximum valid length should be the minimum of the two
        link_valid_length = abs(link_valid_chunk)

        # Which one is now longer? Use this is valid_length
        if link_valid_length < valid_lengths[0]:
            valid_lengths[0] = link_valid_length

    return data_lengths, valid_chunks, valid_lengths


class JobSegmenter(object):
    """ This class is used when running sngl_ifo_job_setup to determine what times
    should be analysed be each job and what data is needed.
    """

    def __init__(self, data_lengths, valid_chunks, valid_lengths, curr_seg,
                 curr_exe_class, compatibility_mode = False):
        """ Initialize class. """
        self.exe_class = curr_exe_class
        self.curr_seg = curr_seg
        self.curr_seg_length = float(abs(curr_seg))

        self.data_length, self.valid_chunk, self.valid_length = \
                      self.pick_tile_size(self.curr_seg_length, data_lengths,
                                                  valid_chunks, valid_lengths)

        self.data_chunk = segments.segment([0, self.data_length])
        self.data_loss = self.data_length - abs(self.valid_chunk)

        if self.data_loss < 0:
            raise ValueError("pycbc.workflow.jobsetup needs fixing! Please contact a developer")

        if self.curr_seg_length < self.data_length:
            self.num_jobs = 0
            return

        # How many jobs do we need
        self.num_jobs = int( math.ceil( (self.curr_seg_length \
                                - self.data_loss) / float(self.valid_length) ))

        self.compatibility_mode = compatibility_mode
        if compatibility_mode and (self.valid_length != abs(self.valid_chunk)):
            errMsg = "In compatibility mode the template bank and matched-"
            errMsg += "filter jobs must read in the same amount of data."
            print(self.valid_length, self.valid_chunk)
            raise ValueError(errMsg)
        elif compatibility_mode and len(data_lengths) > 1:
            raise ValueError("Cannot enable compatibility mode tiling with "
                             "variable tile size")
        elif compatibility_mode:
            # What is the incremental shift between jobs
            self.job_time_shift = self.valid_length
        elif self.curr_seg_length == self.data_length:
            # If the segment length is identical to the data length then I
            # will have exactly 1 job!
            self.job_time_shift = 0
        else:
            # What is the incremental shift between jobs
            self.job_time_shift = (self.curr_seg_length - self.data_length) / \
                                   float(self.num_jobs - 1)

    def pick_tile_size(self, seg_size, data_lengths, valid_chunks, valid_lengths):
        """ Choose job tiles size based on science segment length """

        if len(valid_lengths) == 1:
            return data_lengths[0], valid_chunks[0], valid_lengths[0]
        else:
            # Pick the tile size that is closest to 1/3 of the science segment
            target_size = seg_size / 3
            pick, pick_diff = 0, abs(valid_lengths[0] - target_size)
            for i, size in enumerate(valid_lengths):
                if abs(size - target_size) < pick_diff:
                    pick, pick_diff  = i, abs(size - target_size)
            return data_lengths[pick], valid_chunks[pick], valid_lengths[pick]

    def get_valid_times_for_job(self, num_job, allow_overlap=True):
        """ Get the times for which this job is valid. """
        if self.compatibility_mode:
            return self.get_valid_times_for_job_legacy(num_job)
        else:
            return self.get_valid_times_for_job_workflow(num_job,
                                                   allow_overlap=allow_overlap)

    def get_valid_times_for_job_workflow(self, num_job, allow_overlap=True):
        """ Get the times for which the job num_job will be valid, using workflow's
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
                err_msg = ("Workflow is attempting to generate output "
                          "from a job at times where it is not valid.")
                raise ValueError(err_msg)
            job_valid_seg = segments.segment([lower_boundary,
                                              upper_boundary])
        return job_valid_seg

    def get_valid_times_for_job_legacy(self, num_job):
        """ Get the times for which the job num_job will be valid, using the method
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
        """ Get the data that this job will read in. """
        if self.compatibility_mode:
            job_data_seg =  self.get_data_times_for_job_legacy(num_job)
        else:
            job_data_seg =  self.get_data_times_for_job_workflow(num_job)

        # Sanity check that all data is used
        if num_job == 0:
            if job_data_seg[0] != self.curr_seg[0]:
                err= "Job is not using data from the start of the "
                err += "science segment. It should be using all data."
                raise ValueError(err)
        if num_job == (self.num_jobs - 1):
            if job_data_seg[1] != self.curr_seg[1]:
                err = "Job is not using data from the end of the "
                err += "science segment. It should be using all data."
                raise ValueError(err)

        if hasattr(self.exe_class, 'zero_pad_data_extend'):
            job_data_seg = self.exe_class.zero_pad_data_extend(job_data_seg,
                                                                 self.curr_seg)

        return job_data_seg

    def get_data_times_for_job_workflow(self, num_job):
        """ Get the data that this job will need to read in. """
        # small factor of 0.0001 to avoid float round offs causing us to
        # miss a second at end of segments.
        shift_dur = self.curr_seg[0] + int(self.job_time_shift * num_job\
                                           + 0.0001)
        job_data_seg = self.data_chunk.shift(shift_dur)
        return job_data_seg

    def get_data_times_for_job_legacy(self, num_job):
        """ Get the data that this job will need to read in. """
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
        return job_data_seg

class PyCBCInspiralExecutable(Executable):
    """ The class used to create jobs for pycbc_inspiral Executable. """

    current_retention_level = Executable.ALL_TRIGGERS
    file_input_options = ['--gating-file']

    def __init__(self, cp, exe_name, ifo=None, out_dir=None,
                 injection_file=None, tags=None):
        if tags is None:
            tags = []
        super(PyCBCInspiralExecutable, self).__init__(cp, exe_name, None, ifo,
                                                      out_dir, tags=tags)
        self.cp = cp
        self.set_memory(2000)
        self.injection_file = injection_file
        self.ext = '.hdf'

        self.num_threads = 1
        if self.get_opt('processing-scheme') is not None:
            stxt = self.get_opt('processing-scheme')
            if len(stxt.split(':')) > 1:
                self.num_threads = stxt.split(':')[1]

    def create_node(self, data_seg, valid_seg, parent=None, dfParents=None, tags=None):
        if tags is None:
            tags = []
        node = Node(self)
        if not self.has_opt('pad-data'):
            raise ValueError("The option pad-data is a required option of "
                             "%s. Please check the ini file." % self.name)
        pad_data = int(self.get_opt('pad-data'))

        if not dfParents:
            raise ValueError("%s must be supplied with data file(s)"
                              %(self.name))

        # set remaining options flags
        node.add_opt('--gps-start-time',
                     int_gps_time_to_str(data_seg[0] + pad_data))
        node.add_opt('--gps-end-time',
                     int_gps_time_to_str(data_seg[1] - pad_data))
        node.add_opt('--trig-start-time', int_gps_time_to_str(valid_seg[0]))
        node.add_opt('--trig-end-time', int_gps_time_to_str(valid_seg[1]))
        node.add_profile('condor', 'request_cpus', self.num_threads)

        if self.injection_file is not None:
            node.add_input_opt('--injection-file', self.injection_file)

        # set the input and output files
        fil = node.new_output_file_opt(valid_seg, self.ext, '--output', tags=tags,
                         store_file=self.retain_files, use_tmp_subdirs=True)
        fil.add_metadata('data_seg', data_seg)
        node.add_input_list_opt('--frame-files', dfParents)
        node.add_input_opt('--bank-file', parent, )

        # FIXME: This hack is needed for pipedown compatibility. user-tag is
        #        no-op and is only needed for pipedown to know whether this is
        #        a "FULL_DATA" job or otherwise. Alex wants to burn this code
        #        with fire.
        if node.output_files[0].storage_path is not None:
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
        """ Determine possible dimensions of needed input and valid output
        """

        if self.cp.has_option('workflow-matchedfilter',
                                                      'min-analysis-segments'):
            min_analysis_segs = int(self.cp.get('workflow-matchedfilter',
                                                'min-analysis-segments'))
        else:
            min_analysis_segs = 0

        if self.cp.has_option('workflow-matchedfilter',
                                                      'max-analysis-segments'):
            max_analysis_segs = int(self.cp.get('workflow-matchedfilter',
                                                'max-analysis-segments'))
        else:
            # Choose ridiculously large default value
            max_analysis_segs = 1000

        if self.cp.has_option('workflow-matchedfilter', 'min-analysis-length'):
            min_analysis_length = int(self.cp.get('workflow-matchedfilter',
                                                  'min-analysis-length'))
        else:
            min_analysis_length = 0

        if self.cp.has_option('workflow-matchedfilter', 'max-analysis-length'):
            max_analysis_length = int(self.cp.get('workflow-matchedfilter',
                                                  'max-analysis-length'))
        else:
            # Choose a ridiculously large default value
            max_analysis_length = 100000

        segment_length = int(self.get_opt('segment-length'))
        pad_data = 0
        if self.has_opt('pad-data'):
            pad_data += int(self.get_opt('pad-data'))

        # NOTE: Currently the tapered data is ignored as it is short and
        #       will lie within the segment start/end pad. This means that
        #       the tapered data *will* be used for PSD estimation (but this
        #       effect should be small). It will also be in the data segments
        #       used for SNR generation (when in the middle of a data segment
        #       where zero-padding is not being used) but the templates should
        #       not be long enough to use this data assuming segment start/end
        #       pad take normal values. When using zero-padding this data will
        #       be used for SNR generation.

        #if self.has_opt('taper-data'):
        #    pad_data += int(self.get_opt( 'taper-data' ))
        if self.has_opt('allow-zero-padding'):
            self.zero_padding=True
        else:
            self.zero_padding=False

        start_pad = int(self.get_opt( 'segment-start-pad'))
        end_pad = int(self.get_opt('segment-end-pad'))

        seg_ranges = range(min_analysis_segs, max_analysis_segs + 1)
        data_lengths = []
        valid_regions = []
        for nsegs in seg_ranges:
            analysis_length = (segment_length - start_pad - end_pad) * nsegs
            if not self.zero_padding:
                data_length = analysis_length + pad_data * 2 \
                              + start_pad + end_pad
                start = pad_data + start_pad
                end = data_length - pad_data - end_pad
            else:
                data_length = analysis_length + pad_data * 2
                start = pad_data
                end = data_length - pad_data
            if data_length > max_analysis_length: continue
            if data_length < min_analysis_length: continue
            data_lengths += [data_length]
            valid_regions += [segments.segment(start, end)]
        # If min_analysis_length is given, ensure that it is added as an option
        # for job analysis length.
        if min_analysis_length:
            data_length = min_analysis_length
            if not self.zero_padding:
                start = pad_data + start_pad
                end = data_length - pad_data - end_pad
            else:
                start = pad_data
                end = data_length - pad_data
            if end > start:
                data_lengths += [data_length]
                valid_regions += [segments.segment(start, end)]

        return data_lengths, valid_regions

    def zero_pad_data_extend(self, job_data_seg, curr_seg):
        """When using zero padding, *all* data is analysable, but the setup
        functions must include the padding data where it is available so that
        we are not zero-padding in the middle of science segments. This
        function takes a job_data_seg, that is chosen for a particular node
        and extends it with segment-start-pad and segment-end-pad if that
        data is available.
        """
        if self.zero_padding is False:
            return job_data_seg
        else:
            start_pad = int(self.get_opt( 'segment-start-pad'))
            end_pad = int(self.get_opt('segment-end-pad'))
            new_data_start = max(curr_seg[0], job_data_seg[0] - start_pad)
            new_data_end = min(curr_seg[1], job_data_seg[1] + end_pad)
            new_data_seg = segments.segment([new_data_start, new_data_end])
            return new_data_seg


class PyCBCTmpltbankExecutable(Executable):
    """ The class used to create jobs for pycbc_geom_nonspin_bank Executable and
    any other Executables using the same command line option groups.
    """

    current_retention_level = Executable.MERGED_TRIGGERS
    def __init__(self, cp, exe_name, ifo=None, out_dir=None,
                 tags=None, write_psd=False, psd_files=None):
        if tags is None:
            tags = []
        super(PyCBCTmpltbankExecutable, self).__init__(cp, exe_name, 'vanilla', ifo, out_dir, tags=tags)
        self.cp = cp
        self.set_memory(2000)
        self.write_psd = write_psd
        self.psd_files = psd_files

    def create_node(self, data_seg, valid_seg, parent=None, dfParents=None, tags=None):
        if tags is None:
            tags = []
        node = Node(self)

        if not dfParents:
            raise ValueError("%s must be supplied with data file(s)"
                              % self.name)

        pad_data = int(self.get_opt('pad-data'))
        if pad_data is None:
            raise ValueError("The option pad-data is a required option of "
                             "%s. Please check the ini file." % self.name)

        # set the remaining option flags
        node.add_opt('--gps-start-time',
                     int_gps_time_to_str(data_seg[0] + pad_data))
        node.add_opt('--gps-end-time',
                     int_gps_time_to_str(data_seg[1] - pad_data))

        # set the input and output files
        # Add the PSD file if needed
        if self.write_psd:
            node.new_output_file_opt(valid_seg, '.txt', '--psd-output',
                                     tags=tags+['PSD_FILE'], store_file=self.retain_files)
        node.new_output_file_opt(valid_seg, '.xml.gz', '--output-file',
                                 tags=tags, store_file=self.retain_files)
        node.add_input_list_opt('--frame-files', dfParents)
        return node

    def create_nodata_node(self, valid_seg, tags=None):
        """ A simplified version of create_node that creates a node that does not
        need to read in data.

        Parameters
        -----------
        valid_seg : glue.segment
            The segment over which to declare the node valid. Usually this
            would be the duration of the analysis.

        Returns
        --------
        node : pycbc.workflow.core.Node
            The instance corresponding to the created node.
        """
        if tags is None:
            tags = []
        node = Node(self)

        # Set the output file
        # Add the PSD file if needed
        if self.write_psd:
            node.new_output_file_opt(valid_seg, '.txt', '--psd-output',
                                     tags=tags+['PSD_FILE'], store_file=self.retain_files)

        node.new_output_file_opt(valid_seg, '.xml.gz', '--output-file',
                                 store_file=self.retain_files)

        if self.psd_files is not None:
            should_add = False

            # If any of the ifos for this job are in the set
            # of ifos for which a static psd was provided.
            for ifo in self.ifo_list:
                for psd_file in self.psd_files:
                    if ifo in psd_file.ifo_list:
                        should_add = True

            if should_add:
                node.add_input_opt('--psd-file', psd_file)

        return node

    def get_valid_times(self):
        pad_data = int(self.get_opt( 'pad-data'))
        analysis_length = int(self.cp.get('workflow-tmpltbank', 'analysis-length'))

        #FIXME this should not be hard coded
        data_length = analysis_length + pad_data * 2
        start = pad_data
        end = data_length - pad_data
        return [data_length], [segments.segment(start, end)]

class LigoLWCombineSegsExecutable(Executable):
    """ This class is used to create nodes for the ligolw_combine_segments
    Executable
    """

    # Always want to keep the segments
    current_retention_level = Executable.FINAL_RESULT
    def create_node(self, valid_seg, veto_files, segment_name):
        node = Node(self)
        node.add_opt('--segment-name', segment_name)
        for fil in veto_files:
            node.add_input_arg(fil)
        node.new_output_file_opt(valid_seg, '.xml', '--output',
                                 store_file=self.retain_files)
        return node

class LigolwAddExecutable(Executable):
    """ The class used to create nodes for the ligolw_add Executable. """

    current_retention_level = Executable.INTERMEDIATE_PRODUCT
    def __init__(self, *args, **kwargs):
        super(LigolwAddExecutable, self).__init__(*args, **kwargs)
        self.set_memory(2000)

    def create_node(self, jobSegment, input_files, output=None,
                    use_tmp_subdirs=True, tags=None):
        if tags is None:
            tags = []
        node = Node(self)

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
            node.new_output_file_opt(jobSegment, '.xml.gz', '--output',
                                    tags=tags, store_file=self.retain_files,
                                    use_tmp_subdirs=use_tmp_subdirs)
        return node

class LigolwSSthincaExecutable(Executable):
    """ The class responsible for making jobs for ligolw_sstinca. """

    current_retention_level = Executable.MERGED_TRIGGERS
    def __init__(self, cp, exe_name, universe=None, ifo=None, out_dir=None,
                 dqVetoName=None, tags=None):
        if tags is None:
            tags = []
        super(LigolwSSthincaExecutable, self).__init__(cp, exe_name, universe, ifo, out_dir, tags=tags)
        self.set_memory(2000)
        if dqVetoName:
            self.add_opt("--vetoes-name", dqVetoName)

    def create_node(self, jobSegment, coincSegment, inputFile, tags=None,
                    write_likelihood=False):
        if tags is None:
            tags = []
        node = Node(self)
        node.add_input_arg(inputFile)

        # Add the start/end times
        segString = ""
        if coincSegment[0]:
            segString += str(coincSegment[0])
        segString += ":"
        if coincSegment[1]:
            segString += str(coincSegment[1])

        node.add_opt('--coinc-end-time-segment', segString)

        node.new_output_file_opt(jobSegment, '.xml.gz', '--output-file',
                        tags=self.tags+tags, store_file=self.retain_files,
                        use_tmp_subdirs=True)

        if write_likelihood:
            node.new_output_file_opt(jobSegment, '.xml.gz',
                        '--likelihood-output-file',
                        tags=['DIST_STATS']+self.tags+tags,
                        store_file=self.retain_files)
        return node

class PycbcSqliteSimplifyExecutable(Executable):
    """ The class responsible for making jobs for pycbc_sqlite_simplify. """

    current_retention_level = Executable.INTERMEDIATE_PRODUCT
    def __init__(self, cp, exe_name, universe=None, ifo=None, out_dir=None, tags=None):
        if tags is None:
            tags = []
        super(PycbcSqliteSimplifyExecutable, self).__init__(cp, exe_name, universe, ifo, out_dir, tags=tags)
        self.set_memory(2000)

    def create_node(self, job_segment, inputFiles, injFile=None,
                    injString=None, workflow=None):
        node = Node(self)
        if injFile and not injString:
            raise ValueError("injString needed if injFile supplied.")
        # Need to check if I am dealing with a single split job
        extra_tags = []
        if len(inputFiles) == 1:
            for tag in inputFiles[0].tags:
                if tag.startswith('JOB'):
                    extra_tags.append(tag)
        # If the number of input files is large, split this up
        num_inputs = len(inputFiles)
        if num_inputs > 20 and workflow is not None:
            reduced_inputs = FileList([])
            count = 0
            testing = 0
            curr_node = Node(self)
            # hack the retention level ..
            if self.global_retention_threshold > \
                                 Executable.INTERMEDIATE_PRODUCT:
                curr_store_file = False
            else:
                curr_store_file = True
            for i, file in enumerate(inputFiles):
                curr_node.add_input_arg(file)
                if ( (not (i % 20)) and (i != 0) ) or ((i+1) == num_inputs):
                    count_tag = ['SPLIT%d' %(count)]
                    curr_node.new_output_file_opt(job_segment, '.sqlite',
                                    '--output-file',
                                     tags=self.tags+extra_tags+count_tag,
                                     store_file=curr_store_file)
                    workflow.add_node(curr_node)
                    reduced_inputs.append(curr_node.output_file)
                    testing += len(curr_node._inputs)
                    curr_node = Node(self)
                    count+=1
            inputFiles = reduced_inputs
        elif num_inputs > 20:
            err_msg = "Please provide the workflow keyword to the "
            err_msg += "pycbc_sqlite_simplify create node function if "
            err_msg += "supplying a large number of inputs. This allows us "
            err_msg += "to parallelize the operation and speedup the workflow."
            logging.warn(err_msg)

        for file in inputFiles:
            node.add_input_arg(file)
        if injFile:
            node.add_input_opt("--injection-file", injFile)
            node.add_opt("--simulation-tag", injString)
        node.new_output_file_opt(job_segment, '.sqlite', '--output-file',
                        tags=self.tags+extra_tags, store_file=self.retain_files)
        return node

class SQLInOutExecutable(Executable):
    """
    The class responsible for making jobs for SQL codes taking one input and
    one output.
    """
    current_retention_level=Executable.ALL_TRIGGERS
    def __init__(self, cp, exe_name, universe=None, ifo=None, out_dir=None, tags=None):
        if tags is None:
            tags = []
        super(SQLInOutExecutable, self).__init__(cp, exe_name, universe, ifo, out_dir, tags=tags)

    def create_node(self, job_segment, input_file):
        # Need to follow through the split job tag if present
        extra_tags = []
        for tag in input_file.tags:
            if tag.startswith('JOB'):
                extra_tags.append(tag)

        node = Node(self)
        node.add_input_opt('--input', input_file)
        node.new_output_file_opt(job_segment, '.sqlite', '--output',
                        tags=self.tags+extra_tags, store_file=self.retain_files)
        return node

class PycbcCalculateFarExecutable(SQLInOutExecutable):
    """
    The class responsible for making jobs for the FAR calculation code. This
    only raises the default retention level
    """
    current_retention_level=Executable.FINAL_RESULT

class ExtractToXMLExecutable(Executable):
    """
    This class is responsible for running ligolw_sqlite jobs that will take an
    SQL file and dump it back to XML.
    """
    current_retention_level = Executable.INTERMEDIATE_PRODUCT
    def __init__(self, cp, exe_name, universe=None, ifo=None, out_dir=None,
                 tags=None):
        if tags is None:
            tags = []
        Executable.__init__(self, cp, exe_name, universe, ifo, out_dir,
                                  tags=tags)
    def create_node(self, job_segment, input_file):
        node = Node(self)
        node.add_input_opt('--database', input_file)
        node.new_output_file_opt(job_segment, '.xml', '--extract',
                                 tags=self.tags, store_file=self.retain_files)
        return node

class ComputeDurationsExecutable(SQLInOutExecutable):
    """
    The class responsible for making jobs for pycbc_compute_durations.
    """
    current_retention_level = Executable.INTERMEDIATE_PRODUCT
    def create_node(self, job_segment, input_file, summary_xml_file):
        node = Node(self)
        node.add_input_opt('--input', input_file)
        node.add_input_opt('--segment-file', summary_xml_file)
        node.new_output_file_opt(job_segment, '.sqlite', '--output',
                                 tags=self.tags, store_file=self.retain_files)
        return node

class InspinjfindExecutable(Executable):
    """
    The class responsible for running jobs with pycbc_inspinjfind
    """
    current_retention_level = Executable.INTERMEDIATE_PRODUCT
    def __init__(self, cp, exe_name, universe=None, ifo=None, out_dir=None,
                 tags=None):
        if tags is None:
            tags = []
        Executable.__init__(self, cp, exe_name, universe, ifo, out_dir,
                                  tags=tags)
    def create_node(self, job_segment, input_file):
        node = Node(self)
        node.add_input_opt('--input-file', input_file)
        node.new_output_file_opt(job_segment, '.xml', '--output-file',
                                 tags=self.tags, store_file=self.retain_files)
        return node

class PycbcSplitInspinjExecutable(Executable):
    """
    The class responsible for running the pycbc_split_inspinj executable
    """
    current_retention_level = Executable.INTERMEDIATE_PRODUCT
    def __init__(self, cp, exe_name, num_splits, universe=None, ifo=None,
                 out_dir=None):
        super(PycbcSplitInspinjExecutable, self).__init__(cp, exe_name,
                universe, ifo, out_dir, tags=[])
        self.num_splits = int(num_splits)
    def create_node(self, parent, tags=None):
        if tags is None:
            tags = []
        node = Node(self)

        node.add_input_opt('--input-file', parent)

        if parent.name.endswith("gz"):
            ext = ".xml.gz"
        else:
            ext = ".xml"

        out_files = FileList([])
        for i in range(self.num_splits):
            curr_tag = 'split%d' % i
            curr_tags = parent.tags + [curr_tag]
            job_tag = parent.description + "_" + self.name.upper()
            out_file = File(parent.ifo_list, job_tag, parent.segment,
                            extension=ext, directory=self.out_dir,
                            tags=curr_tags, store_file=self.retain_files)
            out_files.append(out_file)

        node.add_output_list_opt('--output-files', out_files)
        return node

class PycbcPickleHorizonDistsExecutable(Executable):
    """
    The class responsible for running the pycbc_pickle_horizon_distances
    executable which is part 1 of 4 of the gstlal_inspiral_calc_likelihood
    functionality
    """
    # FIXME: This class will soon be removed when gstlal post-proc is updated
    current_retention_level = Executable.INTERMEDIATE_PRODUCT
    def __init__(self, cp, exe_name, universe=None, ifo=None, out_dir=None,
                 tags=None):
        if tags is None:
            tags = []
        Executable.__init__(self, cp, exe_name, universe, ifo, out_dir,
                                  tags=tags)
    def create_node(self, job_segment, trigger_files):
        node = Node(self)
        for file in trigger_files:
            node.add_input_arg(file)
        node.new_output_file_opt(job_segment, '.pickle', '--output-file',
                                 tags=self.tags, store_file=self.retain_files)
        return node

class PycbcCombineLikelihoodExecutable(Executable):
    """
    The class responsible for running the pycbc_combine_likelihood
    executable which is part 2 of 4 of the gstlal_inspiral_calc_likelihood
    functionality
    """
    current_retention_level = Executable.INTERMEDIATE_PRODUCT
    def __init__(self, cp, exe_name, universe=None, ifo=None, out_dir=None,
                 tags=None):
        if tags is None:
            tags = []
        Executable.__init__(self, cp, exe_name, universe, ifo, out_dir,
                                  tags=tags)
    def create_node(self, job_segment, likelihood_files, horizon_dist_file):
        node = Node(self)
        node.add_input_list_opt('--likelihood-urls', likelihood_files)
        node.add_input_opt('--horizon-dist-file', horizon_dist_file)
        node.new_output_file_opt(job_segment, '.xml.gz', '--output-file',
                                 tags=self.tags, store_file=self.retain_files)
        return node

class PycbcGenerateRankingDataExecutable(Executable):
    """
    The class responsible for running the pycbc_gen_ranking_data
    executable which is part 3 of 4 of the gstlal_inspiral_calc_likelihood
    functionality
    """
    current_retention_level = Executable.INTERMEDIATE_PRODUCT
    def __init__(self, cp, exe_name, universe=None, ifo=None, out_dir=None,
                 tags=None):
        if tags is None:
            tags = []
        Executable.__init__(self, cp, exe_name, universe, ifo, out_dir,
                                  tags=tags)
        self.set_num_cpus(4)
        self.set_memory('8000')
    def create_node(self, job_segment, likelihood_file, horizon_dist_file):
        node = Node(self)
        node.add_input_opt('--likelihood-file', likelihood_file)
        node.add_input_opt('--horizon-dist-file', horizon_dist_file)
        node.new_output_file_opt(job_segment, '.xml.gz', '--output-file',
                                 tags=self.tags, store_file=self.retain_files)
        return node

class PycbcCalculateLikelihoodExecutable(Executable):
    """
    The class responsible for running the pycbc_calculate_likelihood
    executable which is part 4 of 4 of the gstlal_inspiral_calc_likelihood
    functionality
    """
    current_retention_level = Executable.FINAL_RESULT
    def __init__(self, cp, exe_name, universe=None, ifo=None, out_dir=None,
                 tags=None):
        if tags is None:
            tags = []
        Executable.__init__(self, cp, exe_name, universe, ifo, out_dir,
                                  tags=tags)
    def create_node(self, job_segment, trigger_file, likelihood_file,
                    horizon_dist_file):
        node = Node(self)
        # Need to follow through the split job tag if present
        extra_tags = []
        for tag in trigger_file.tags:
            if tag.startswith('JOB'):
                extra_tags.append(tag)

        node.add_input_opt('--trigger-file', trigger_file)
        node.add_input_opt('--horizon-dist-file', horizon_dist_file)
        node.add_input_opt('--likelihood-file', likelihood_file)
        node.new_output_file_opt(job_segment, '.sqlite', '--output-file',
                                 tags=self.tags + extra_tags,
                                 store_file=self.retain_files)
        return node

class GstlalMarginalizeLikelihoodExecutable(Executable):
    """
    The class responsible for running the gstlal marginalize_likelihood jobs
    """
    current_retention_level = Executable.FINAL_RESULT
    def __init__(self, cp, exe_name, universe=None, ifo=None, out_dir=None,
                 tags=None):
        if tags is None:
            tags = []
        Executable.__init__(self, cp, exe_name, universe, ifo, out_dir,
                                  tags=tags)
    def create_node(self, job_segment, input_file):
        node = Node(self)
        node.add_input_arg(input_file)
        node.new_output_file_opt(job_segment, '.xml.gz', '--output',
                                 tags=self.tags, store_file=self.retain_files)
        return node

class GstlalFarfromsnrchisqhistExecutable(Executable):
    """
    The class responsible for running the gstlal far from chisq hist jobs
    """
    current_retention_level = Executable.FINAL_RESULT
    def __init__(self, cp, exe_name, universe=None, ifo=None, out_dir=None,
                 tags=None):
        if tags is None:
            tags = []
        Executable.__init__(self, cp, exe_name, universe, ifo, out_dir,
                                  tags=tags)
    def create_node(self, job_segment, non_inj_db, marg_input_file,
                   inj_database=None, write_background_bins=False):
        node = Node(self)
        node.add_input_opt("--non-injection-db", non_inj_db)
        if inj_database is not None:
            node.add_input_opt("--input-database", inj_database)
        node.add_input_opt("--background-bins-file", marg_input_file)
        node.new_output_file_opt(job_segment, '.sqlite', '--output-database',
                                 tags=self.tags, store_file=self.retain_files)
        # FIXME: Not supported yet
        if write_background_bins:
            node.new_output_file_opt(job_segment, '.xml.gz',
                                  "--background-bins-out-file",
                                  tags=["POSTMARG"] + self.tags,
                                  store_file=self.retain_files)
        return node

class GstlalPlotSensitivity(Executable):
    """
    The class responsible for running gstlal_plot_sensitivity
    """
    current_retention_level = Executable.FINAL_RESULT
    def __init__(self, cp, exe_name, universe=None, ifo=None, out_dir=None,
                 tags=None):
        if tags is None:
            tags = []
        Executable.__init__(self, cp, exe_name, universe, ifo, out_dir,
                                  tags=tags)
        self.set_memory('4000')
    def create_node(self, non_inj_db, injection_dbs):
        node = Node(self)
        node.add_input_opt("--zero-lag-database", non_inj_db)
        for file in injection_dbs:
            node.add_input_arg(file)
        return node

class GstlalPlotSummary(Executable):
    """
    The class responsible for running gstlal_plot_summary
    """
    current_retention_level = Executable.FINAL_RESULT
    def __init__(self, cp, exe_name, universe=None, ifo=None, out_dir=None,
                 tags=None):
        if tags is None:
            tags = []
        Executable.__init__(self, cp, exe_name, universe, ifo, out_dir,
                                  tags=tags)
        self.set_memory('4000')
    def create_node(self, non_inj_db, injection_dbs):
        node = Node(self)
        node.add_input_arg(non_inj_db)
        for file in injection_dbs:
            node.add_input_arg(file)
        return node

class GstlalPlotBackground(Executable):
    """
    The class responsible for running gstlal_plot_background
    """
    current_retention_level = Executable.FINAL_RESULT
    def __init__(self, cp, exe_name, universe=None, ifo=None, out_dir=None,
                 tags=None):
        if tags is None:
            tags = []
        Executable.__init__(self, cp, exe_name, universe, ifo, out_dir,
                                  tags=tags)
        self.set_memory('4000')
    def create_node(self, non_inj_db, likelihood_file):
        node = Node(self)
        node.add_input_opt("--database", non_inj_db)
        node.add_input_arg(likelihood_file)
        return node

class GstlalSummaryPage(Executable):
    """
    The class responsible for running gstlal_inspiral_summary_page
    """
    current_retention_level = Executable.FINAL_RESULT
    def __init__(self, cp, exe_name, universe=None, ifo=None, out_dir=None,
                 tags=None):
        if tags is None:
            tags = []
        Executable.__init__(self, cp, exe_name, universe, ifo, out_dir,
                                  tags=tags)
    def create_and_add_node(self, workflow, parent_nodes):
        node = Node(self)
        # FIXME: As the parents of this job (the plotting scripts) do not track
        # the location of their output files, we must set explicit parent-child
        # relations the "old-fashioned" way here. Possibly this is how the AEI
        # DB stuff could work?
        workflow.add_node(node)
        for parent in parent_nodes:
            dep = dax.Dependency(parent=parent._dax_node, child=node._dax_node)
            workflow._adag.addDependency(dep)
        return node

class LalappsInspinjExecutable(Executable):
    """
    The class used to create jobs for the lalapps_inspinj Executable.
    """
    current_retention_level = Executable.FINAL_RESULT
    def create_node(self, segment, exttrig_file=None, tags=None):
        if tags is None:
            tags = []
        node = Node(self)

        curr_tags = self.tags + tags
        # This allows the desired number of injections to be given explicitly
        # in the config file. Used for coh_PTF as segment length is unknown
        # before run time.
        if self.get_opt('write-compress') is not None:
            ext = '.xml.gz'
        else:
            ext = '.xml'

        # Check if these injections are using trigger information to choose
        # sky positions for the simulated signals
        if (self.get_opt('l-distr') == 'exttrig' and exttrig_file is not None \
                and 'trigger' in exttrig_file.description):
            # Use an XML file containing trigger information
            triggered = True
            node.add_input_opt('--exttrig-file', exttrig_file)
        elif (self.get_opt('l-distr') == 'ipn' and exttrig_file is not None \
                and 'IPN' in exttrig_file.description):
            # Use an IPN sky points file
            triggered = True
            node.add_input_opt('--ipn-file', exttrig_file)
        elif (self.get_opt('l-distr') != 'exttrig') \
                and (self.get_opt('l-distr') != 'ipn' and not \
                     self.has_opt('ipn-file')):
            # Use no trigger information for generating injections
            triggered = False
        else:
            err_msg = "The argument 'l-distr' passed to the "
            err_msg += "%s job has the value " % self.tagged_name
            err_msg += "'%s' but you have not " % self.get_opt('l-distr')
            err_msg += "provided the corresponding ExtTrig or IPN file. "
            err_msg += "Please check your configuration files and try again."
            raise ValueError(err_msg)

        if triggered:
            num_injs = int(self.cp.get_opt_tags('workflow-injections',
                                                'num-injs', curr_tags))
            inj_tspace = float(segment[1] - segment[0]) / num_injs
            node.add_opt('--time-interval', inj_tspace)
            node.add_opt('--time-step', inj_tspace)

        node.new_output_file_opt(segment, ext, '--output',
                                 store_file=self.retain_files)

        node.add_opt('--gps-start-time', int_gps_time_to_str(segment[0]))
        node.add_opt('--gps-end-time', int_gps_time_to_str(segment[1]))
        return node

class PycbcDarkVsBrightInjectionsExecutable(Executable):
    """
    The clase used to create jobs for the pycbc_dark_vs_bright_injections Executable.
    """
    current_retention_level = Executable.FINAL_RESULT
    def __init__(self, cp, exe_name, universe=None, ifos=None, out_dir=None,
                 tags=None):
        if tags is None:
            tags = []
        Executable.__init__(self, cp, exe_name, universe, ifos, out_dir,
                            tags=tags)
        self.cp = cp
        self.out_dir = out_dir
        self.exe_name = exe_name

    def create_node(self, parent, segment, tags=None):
        if tags is None:
            tags = []
        node = Node(self)
        if not parent:
            raise ValueError("Must provide an input file.")

        node = Node(self)
        # Standard injection file produced by lalapps_inspinj
        # becomes the input here
        node.add_input_opt('-i', parent)
        if self.has_opt('write-compress'):
            ext = '.xml.gz'
        else:
            ext = '.xml'
        # The output are two files:
        # 1) the list of potentially EM bright injections
        tag=['POTENTIALLY_BRIGHT']
        node.new_output_file_opt(segment, ext, '--output-bright',
                                 store_file=self.retain_files, tags=tag)
        # 2) the list of EM dim injections
        tag=['DIM_ONLY']
        node.new_output_file_opt(segment,
                                 ext, '--output-dim',
                                 store_file=self.retain_files, tags=tag)
        return node

class LigolwCBCJitterSkylocExecutable(Executable):
    """
    The class used to create jobs for the ligolw_cbc_skyloc_jitter executable.
    """
    current_retention_level = Executable.MERGED_TRIGGERS
    def __init__(self, cp, exe_name, universe=None, ifos=None, out_dir=None,
                 tags=None):
        if tags is None:
            tags = []
        Executable.__init__(self, cp, exe_name, universe, ifos, out_dir,
                            tags=tags)
        self.cp = cp
        self.out_dir = out_dir
        self.exe_name = exe_name

    def create_node(self, parent, segment, tags=None):
        if tags is None:
            tags = []
        if not parent:
            raise ValueError("Must provide an input file.")

        node = Node(self)
        node.add_input_opt('--input-file', parent)
        output_file = File(parent.ifo_list, self.exe_name,
                           segment, extension='.xml', store_file=True,
                           directory=self.out_dir, tags=tags)
        node.add_output_opt('--output-file', output_file)

        return node

class LigolwCBCAlignTotalSpinExecutable(Executable):
    """
    The class used to create jobs for the ligolw_cbc_skyloc_jitter executable.
    """
    current_retention_level = Executable.MERGED_TRIGGERS
    def __init__(self, cp, exe_name, universe=None, ifos=None, out_dir=None,
                 tags=None):
        if tags is None:
            tags = []
        Executable.__init__(self, cp, exe_name, universe, ifos, out_dir,
                            tags=tags)
        self.cp = cp
        self.out_dir = out_dir
        self.exe_name = exe_name

    def create_node(self, parent, segment, tags=None):
        if tags is None:
            tags = []
        if not parent:
            raise ValueError("Must provide an input file.")

        node = Node(self)
        output_file = File(parent.ifo_list, self.exe_name, segment,
                           extension='.xml', store_file=self.retain_files,
                           directory=self.out_dir, tags=tags)
        node.add_output_opt('--output-file', output_file)
        node.add_input_arg(parent)
        return node

class PycbcTimeslidesExecutable(Executable):
    """
    The class used to create jobs for the pycbc_timeslides Executable.
    """
    current_retention_level = Executable.FINAL_RESULT
    def create_node(self, segment):
        node = Node(self)

        node.new_output_file_opt(segment, '.xml.gz', '--output-files',
                                 store_file=self.retain_files)
        return node

class PycbcSplitBankExecutable(Executable):
    """ The class responsible for creating jobs for pycbc_hdf5_splitbank. """

    extension = '.hdf'
    current_retention_level = Executable.ALL_TRIGGERS
    def __init__(self, cp, exe_name, num_banks,
                 ifo=None, out_dir=None, universe=None):
        super(PycbcSplitBankExecutable, self).__init__(cp, exe_name, universe,
                ifo, out_dir, tags=[])
        self.num_banks = int(num_banks)

    def create_node(self, bank, tags=None):
        """
        Set up a CondorDagmanNode class to run splitbank code

        Parameters
        ----------
        bank : pycbc.workflow.core.File
            The File containing the template bank to be split

        Returns
        --------
        node : pycbc.workflow.core.Node
            The node to run the job
        """
        if tags is None:
            tags = []
        node = Node(self)
        node.add_input_opt('--bank-file', bank)

        # Get the output (taken from inspiral.py)
        out_files = FileList([])
        for i in range( 0, self.num_banks):
            curr_tag = 'bank%d' %(i)
            # FIXME: What should the tags actually be? The job.tags values are
            #        currently ignored.
            curr_tags = bank.tags + [curr_tag] + tags
            job_tag = bank.description + "_" + self.name.upper()
            out_file = File(bank.ifo_list, job_tag, bank.segment,
                            extension=self.extension, directory=self.out_dir,
                            tags=curr_tags, store_file=self.retain_files)
            out_files.append(out_file)
        node.add_output_list_opt('--output-filenames', out_files)
        return node

class PycbcSplitBankXmlExecutable(PycbcSplitBankExecutable):
    """ Subclass resonsible for creating jobs for pycbc_splitbank. """

    extension='.xml.gz'

class PycbcConditionStrainExecutable(Executable):
    """ The class responsible for creating jobs for pycbc_condition_strain. """

    current_retention_level = Executable.ALL_TRIGGERS
    def __init__(self, cp, exe_name, ifo=None, out_dir=None, universe=None,
            tags=None):
        super(PycbcConditionStrainExecutable, self).__init__(cp, exe_name, universe,
              ifo, out_dir, tags)

    def create_node(self, input_files, tags=None):
        if tags is None:
            tags = []
        node = Node(self)
        start_time = self.cp.get("workflow", "start-time")
        end_time = self.cp.get("workflow", "end-time")
        node.add_opt('--gps-start-time', start_time)
        node.add_opt('--gps-end-time', end_time)
        node.add_input_list_opt('--frame-files', input_files)

        out_file = File(self.ifo, "gated",
                        segments.segment(int(start_time), int(end_time)),
                        directory=self.out_dir, store_file=self.retain_files,
                        extension=input_files[0].name.split('.', 1)[-1],
                        tags=tags)
        node.add_output_opt('--output-strain-file', out_file)

        out_gates_file = File(self.ifo, "output_gates",
                              segments.segment(int(start_time), int(end_time)),
                              directory=self.out_dir, extension='txt',
                              store_file=self.retain_files, tags=tags)
        node.add_output_opt('--output-gates-file', out_gates_file)

        return node, out_file

