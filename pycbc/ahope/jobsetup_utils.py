import math
from glue import segments
from pycbc.ahope import AhopeFile
from pycbc.ahope.legacy_ihope import *

def select_tmpltbankjob_instance(curr_exe, curr_section):
    """This function returns an instance of the class that is appropriate for
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
        * exe_class.get_valid_times(ifo, )
        * exe_class.create_job()
        * exe_class.create_node()
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
    """This function returns an instance of the class that is appropriate for
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
        * exe_class.get_valid_times()
        * exe_class.create_condorjob()
        * exe_class.create_condornode()
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
    """This function returns an instance of the class that is appropriate for
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
        * exe_class.create_node()
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

def sngl_ifo_job_setup(workflow, ifo, out_files, exe_instance, science_segs, 
                       datafind_outs, output_dir, parents=None, 
                       link_exe_instance=False, allow_overlap=True):
    """
    This function sets up a set of single ifo jobs.

    Parameters
    -----------
    workflow: Workflow
        An instanced class that manages the constructed workflow.
    ifo : string
        The name of the ifo to set up the jobs for
    out_files : AhopeOutFileList or AhopeOutGroupList
        The AhopeOutFileList containing the list of jobs. Jobs will be appended
        to this list, and it does not need to be empty when supplied.
    exe_instance : Instanced class
        An instanced class that contains the functions needed to set up things
        that are specific to the executable being run.
    science_segs : segments.segmentlist
        The list of times that the jobs should cover
    datafind_outs : AhopeFileList
        The file list containing the datafind files.
    output_dir : path string
        The directory where data products will be placed.
    parents : AhopeFileList (optional, kwarg, default=None)
        The AhopeOutFileList containing the list of jobs that are parents to
        the one being set up.
    link_exe_instance : Executable instance (optional),
        Coordinate the valid times with another executable.
    allow_overlap : boolean (optional, kwarg, default = True)
        If this is set the times that jobs are valid for will be allowed to
        overlap. This may be desired for template banks which may have some
        overlap in the times they cover. This may not be desired for inspiral
        jobs, where you probably want triggers recorded by jobs to not overlap
        at all.
    """
    cp = workflow.cp
    
    # Begin by getting analysis start and end, and start and end of time
    # that the output file is valid for
    data_length, valid_chunk = exe_instance.get_valid_times(cp, ifo)

    data_chunk = segments.segment([0, data_length])
    job_tag = exe_instance.exe_name.upper()
    
    if link_exe_instance:
        _, link_valid_chunk = link_exe_instance.get_valid_times(cp, ifo)
        valid_chunk_start = max(valid_chunk[0], link_valid_chunk[0])
        valid_chunk_end = min(valid_chunk[1], link_valid_chunk[1])
        valid_chunk = segments.segment([valid_chunk_start, valid_chunk_end])

    # Set up the condorJob class for the current executable
    curr_exe_job = exe_instance.create_job(cp, ifo, output_dir)
    
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
                 (curr_seg_length - data_loss) / float(abs(valid_chunk)) ))
        # What is the incremental shift between jobs
        time_shift = (curr_seg_length - data_length) / float(num_jobs - 1)
        for job_num in range(num_jobs):
            # Get the science segment for this job
            shift_dur = curr_seg[0] + int(time_shift * job_num)
            job_data_seg = data_chunk.shift(shift_dur)
            job_valid_seg = valid_chunk.shift(shift_dur)
            # If we need to recalculate the valid times to avoid overlap
            if not allow_overlap:
                data_per_job = (curr_seg_length - data_loss) / float(num_jobs)
                lower_boundary = int(job_num*data_per_job +
                                     valid_chunk[0] + curr_seg[0])
                upper_boundary = int(data_per_job + lower_boundary)
                if lower_boundary < job_valid_seg[0] or \
                        upper_boundary > job_valid_seg[1]:
                    err_msg = ("Ahope is attempting to generate output "
                              "from a job at times where it is not valid.")
                    raise ValueError(err_msg)
                job_valid_seg = segments.segment([lower_boundary, 
                                                  upper_boundary])
                
            if parents:
                curr_parent = parents.find_output(ifo, job_valid_seg)
                if not curr_parent:
                    err_string = ("No parent jobs found overlapping %d to %d." 
                                  %(job_valid_seg[0], job_valid_seg[1]))
                    err_string += "\nThis is a bad error! Contact a developer."
                    raise ValueError(err_string)
            else:
                curr_parent = None

            if datafind_outs:
                curr_dfouts = datafind_outs.find_all_output_in_range(ifo, 
                                                                  job_data_seg)
                if not curr_dfouts:
                    err_str = ("No datafind jobs found overlapping %d to %d."
                                %(job_data_seg[0],job_data_seg[1]))
                    err_str += "\nThis shouldn't happen. Contact a developer."
                    raise ValueError(err_str)

            node = curr_exe_job.create_node(job_data_seg, job_valid_seg, 
                                            parent=curr_parent,
                                            dfParents=curr_dfouts)
            workflow.add_node(node)
            out_files += node.output_files
    return out_files

