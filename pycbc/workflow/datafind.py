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
This module is responsible for querying a datafind server to determine the
availability of the data that the code is attempting to run on. It also
performs a number of tests and can act on these as described below. Full
documentation for this function can be found here:
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/ahope/datafind.html
"""

import copy
import logging
import os
import urllib.parse

import igwn_segments as segments
from gwdatafind import find_urls as find_frame_urls
from igwn_ligolw import ligolw, utils

from pycbc.io.ligolw import LIGOLWContentHandler
from pycbc.workflow.core import File, FileList, SegFile, make_analysis_dir

# NOTE urllib is weird. For some reason it only allows known schemes and will
# give *wrong* results, rather then failing, if you use something like gsiftp
# We can add schemes explicitly, as below, but be careful with this!
# (urllib is used indirectly through lal.Cache objects)
urllib.parse.uses_relative.append("osdf")
urllib.parse.uses_netloc.append("osdf")

logger = logging.getLogger(__name__)


def setup_datafind_workflow(
    workflow, science_segs, output_dir, seg_file=None, tags=None
):
    """
    Setup datafind section of the workflow. This section is responsible for
    generating, or setting up the workflow to generate, a list of files that
    record the location of the frame files needed to perform the analysis.
    There could be multiple options here, the datafind jobs could be done at
    run time or could be put into a dag. The subsequent jobs will know
    what was done here from the OutFileList containing the datafind jobs
    (and the Dagman nodes if appropriate.
    For now the only implemented option is to generate the datafind files at
    runtime. This module can also check if the frameFiles actually exist, check
    whether the obtained segments line up with the original ones and update the
    science segments to reflect missing data files.

    Parameters
    ----------
    workflow: pycbc.workflow.core.Workflow
        The workflow class that stores the jobs that will be run.
    science_segs : Dictionary of ifo keyed igwn_segments.segmentlist instances
        This contains the times that the workflow is expected to analyse.
    output_dir : path
        All output files written by datafind processes will be written to this
        directory.
    seg_file : SegFile, optional (default=None)
        The file returned by get_science_segments containing the science
        segments and the associated segment_summary. This will
        be used for the segment_summary test and is required if, and only if,
        performing that test.
    tags : list of string, optional (default=None)
        Use this to specify tags. This can be used if this module is being
        called more than once to give call specific configuration (by setting
        options in [workflow-datafind-${TAG}] rather than [workflow-datafind]).
        This is also used to tag the Files returned by the class to uniqueify
        the Files and uniqueify the actual filename.
        FIXME: Filenames may not be unique with current codes!

    Returns
    -------
    datafindOuts : OutGroupList
        List of all the datafind output files for use later in the pipeline.
    sci_avlble_file : SegFile
        SegFile containing the analysable time after checks in the datafind
        module are applied to the input segment list. For production runs this
        is expected to be equal to the input segment list.
    science_segs : Dictionary of ifo keyed igwn_segments.segmentlist instances
        This contains the times that the workflow is expected to analyse. If
        the updateSegmentTimes kwarg is given this will be updated to reflect
        any instances of missing data.
    sci_avlble_name : string
        The name with which the analysable time is stored in the
        sci_avlble_file.

    """
    if tags is None:
        tags = []
    logger.info("Entering datafind module")
    make_analysis_dir(output_dir)
    cp = workflow.cp

    # Parse for options in ini file
    datafind_method = cp.get_opt_tags("workflow-datafind", "datafind-method", tags)

    if cp.has_option_tags("workflow-datafind", "datafind-check-segment-gaps", tags):
        check_segment_gaps = cp.get_opt_tags(
            "workflow-datafind", "datafind-check-segment-gaps", tags
        )
    else:
        check_segment_gaps = "no_test"
    if cp.has_option_tags("workflow-datafind", "datafind-check-frames-exist", tags):
        check_frames_exist = cp.get_opt_tags(
            "workflow-datafind", "datafind-check-frames-exist", tags
        )
    else:
        check_frames_exist = "no_test"
    if cp.has_option_tags("workflow-datafind", "datafind-check-segment-summary", tags):
        check_segment_summary = cp.get_opt_tags(
            "workflow-datafind", "datafind-check-segment-summary", tags
        )
    else:
        check_segment_summary = "no_test"

    logger.info("Starting datafind with setup_datafind_runtime_generated")
    if datafind_method == "AT_RUNTIME_MULTIPLE_CACHES":
        datafindcaches, datafindouts = setup_datafind_runtime_cache_multi_calls_perifo(
            cp, science_segs, output_dir, tags=tags
        )
    elif datafind_method == "AT_RUNTIME_SINGLE_CACHES":
        datafindcaches, datafindouts = setup_datafind_runtime_cache_single_call_perifo(
            cp, science_segs, output_dir, tags=tags
        )
    elif datafind_method == "AT_RUNTIME_MULTIPLE_FRAMES":
        datafindcaches, datafindouts = setup_datafind_runtime_frames_multi_calls_perifo(
            cp, science_segs, output_dir, tags=tags
        )
    elif datafind_method == "AT_RUNTIME_SINGLE_FRAMES":
        datafindcaches, datafindouts = setup_datafind_runtime_frames_single_call_perifo(
            cp, science_segs, output_dir, tags=tags
        )
    elif datafind_method == "AT_RUNTIME_FAKE_DATA":
        pass
    elif datafind_method == "FROM_PREGENERATED_LCF_FILES":
        ifos = science_segs.keys()
        datafindcaches, datafindouts = setup_datafind_from_pregenerated_lcf_files(
            cp, ifos, output_dir, tags=tags
        )
    else:
        msg = """Entry datafind-method in [workflow-datafind] does not have "
              expected value. Valid values are
              AT_RUNTIME_MULTIPLE_FRAMES, AT_RUNTIME_SINGLE_FRAMES
              AT_RUNTIME_MULTIPLE_CACHES, AT_RUNTIME_SINGLE_CACHES,
              FROM_PREGENERATED_LCF_FILES, or AT_RUNTIME_FAKE_DATA.
              Consult the documentation for more info."""
        raise ValueError(msg)

    using_backup_server = False
    if (
        datafind_method == "AT_RUNTIME_MULTIPLE_FRAMES"
        or datafind_method == "AT_RUNTIME_SINGLE_FRAMES"
    ):
        if cp.has_option_tags(
            "workflow-datafind", "datafind-backup-datafind-server", tags
        ):
            using_backup_server = True
            backup_server = cp.get_opt_tags(
                "workflow-datafind", "datafind-backup-datafind-server", tags
            )
            cp_new = copy.deepcopy(cp)
            cp_new.set(
                "workflow-datafind", "datafind-ligo-datafind-server", backup_server
            )
            cp_new.set("datafind", "urltype", "gsiftp")
            backup_datafindcaches, backup_datafindouts = (
                setup_datafind_runtime_frames_single_call_perifo(
                    cp_new, science_segs, output_dir, tags=tags
                )
            )
            backup_datafindouts = datafind_keep_unique_backups(
                backup_datafindouts, datafindouts
            )
            datafindcaches.extend(backup_datafindcaches)
            datafindouts.extend(backup_datafindouts)

    logger.info("setup_datafind_runtime_generated completed")
    # If we don't have frame files covering all times we can update the science
    # segments.
    if check_segment_gaps in ["warn", "update_times", "raise_error"]:
        logger.info("Checking science segments against datafind output....")
        new_science_segs = get_science_segs_from_datafind_outs(datafindcaches)
        logger.info("New segments calculated from data find output.....")
        missing_data = False
        for ifo in science_segs.keys():
            # If no science segments in input then do nothing
            if not science_segs[ifo]:
                msg = f"No science segments are present for ifo {ifo}, "
                msg += "the segment metadata indicates there is no analyzable"
                msg += " strain data between the selected GPS start and end "
                msg += "times."
                logger.warning(msg)
                continue
            if ifo not in new_science_segs:
                msg = "No data frames were found corresponding to the science "
                msg += f"segments for ifo {ifo}"
                logger.error(msg)
                missing_data = True
                if check_segment_gaps == "update_times":
                    science_segs[ifo] = segments.segmentlist()
                continue
            missing = science_segs[ifo] - new_science_segs[ifo]
            if abs(missing):
                msg = f"From ifo {ifo} we are missing frames covering:"
                msg += "\n{}".format("\n".join(map(str, missing)))
                missing_data = True
                logger.error(msg)
                if check_segment_gaps == "update_times":
                    # Remove missing time, so that we can carry on if desired
                    logger.info("Updating science segments for ifo %s.", ifo)
                    science_segs[ifo] = science_segs[ifo] - missing

        if check_segment_gaps == "raise_error" and missing_data:
            raise ValueError("Workflow cannot find needed data, exiting.")
        logger.info("Done checking, any discrepancies are reported above.")
    elif check_segment_gaps == "no_test":
        pass
    else:
        err_msg = "check_segment_gaps kwarg must take a value from 'no_test', "
        err_msg += "'warn', 'update_times' or 'raise_error'."
        raise ValueError(err_msg)

    # Do all of the frame files that were returned actually exist?
    if check_frames_exist in ["warn", "update_times", "raise_error"]:
        logger.info("Verifying that all frames exist on disk.")
        missing_fr_segs, missing_frames = get_missing_segs_from_frame_file_cache(
            datafindcaches
        )
        missing_flag = False
        for ifo in missing_frames.keys():
            # If no data in the input then do nothing
            if not science_segs[ifo]:
                continue
            # If using a backup server, does the frame exist remotely?
            if using_backup_server:
                # WARNING: This will be slow, but hopefully it will not occur
                #          for too many frames. This could be optimized if
                #          it becomes necessary.
                new_list = []
                for frame in missing_frames[ifo]:
                    for dfout in datafindouts:
                        dfout_pfns = list(dfout.pfns)
                        dfout_urls = [a.url for a in dfout_pfns]
                        if frame.url in dfout_urls:
                            pfn = dfout_pfns[dfout_urls.index(frame.url)]
                            dfout.removePFN(pfn)
                            if len(dfout.pfns) == 0:
                                new_list.append(frame)
                            else:
                                msg = f"Frame {frame.url} not found locally. "
                                msg += (
                                    "Replacing with remote url(s) "
                                    f"{[a.url for a in dfout.pfns]!s}."
                                )
                                logger.info(msg)
                            break
                    else:
                        new_list.append(frame)
                missing_frames[ifo] = new_list
            if missing_frames[ifo]:
                msg = f"From ifo {ifo} we are missing the following frames:"
                msg += "\n".join([a.url for a in missing_frames[ifo]])
                missing_flag = True
                logger.error(msg)
            if check_frames_exist == "update_times":
                # Remove missing times, so that we can carry on if desired
                logger.info("Updating science times for ifo %s.", ifo)
                science_segs[ifo] = science_segs[ifo] - missing_fr_segs[ifo]

        if check_frames_exist == "raise_error" and missing_flag:
            raise ValueError("Workflow cannot find all frames, exiting.")
        logger.info("Finished checking frames.")
    elif check_frames_exist == "no_test":
        pass
    else:
        err_msg = "check_frames_exist kwarg must take a value from 'no_test', "
        err_msg += "'warn', 'update_times' or 'raise_error'."
        raise ValueError(err_msg)

    # Check if there are cases where frames exist, but no entry in the segment
    # summary table are present.
    if check_segment_summary in ["warn", "raise_error"]:
        logger.info("Checking the segment summary table against frames.")
        df_science_segs = get_science_segs_from_datafind_outs(datafindcaches)
        missing_flag = False
        # NOTE: Should this be overrideable in the config file?
        sci_seg_name = "SCIENCE"
        if seg_file is None:
            err_msg = "You must provide the science segments SegFile object "
            err_msg += "if using the datafind-check-segment-summary option."
            raise ValueError(err_msg)
        if seg_file.seg_summ_dict is None:
            err_msg = "The provided science segments SegFile object must "
            err_msg += "contain a valid segment_summary table if using the "
            err_msg += "datafind-check-segment-summary option."
            raise ValueError(err_msg)
        seg_summary_times = seg_file.seg_summ_dict
        for ifo in df_science_segs.keys():
            curr_seg_summ_times = seg_summary_times[ifo + ":" + sci_seg_name]
            missing = df_science_segs[ifo] & seg_file.valid_segments
            missing.coalesce()
            missing = missing - curr_seg_summ_times
            missing.coalesce()
            science_but_not_frame = science_segs[ifo] - df_science_segs[ifo]
            science_but_not_frame.coalesce()
            missing2 = science_segs[ifo] - science_but_not_frame
            missing2.coalesce()
            missing2 = missing2 - curr_seg_summ_times
            missing2.coalesce()
            if abs(missing):
                msg = f"From ifo {ifo} the following times have frames, "
                msg += "but are not covered in the segment summary table."
                msg += "\n{}".format("\n".join(map(str, missing)))
                logger.error(msg)
                missing_flag = True
            if abs(missing2):
                msg = f"From ifo {ifo} the following times have frames, "
                msg += "are science, and are not covered in the segment "
                msg += "summary table."
                msg += "\n{}".format("\n".join(map(str, missing2)))
                logger.error(msg)
                missing_flag = True
        if check_segment_summary == "raise_error" and missing_flag:
            err_msg = "Segment_summary discrepancy detected, exiting."
            raise ValueError(err_msg)
    elif check_segment_summary == "no_test":
        pass
    else:
        err_msg = "check_segment_summary kwarg must take a value from 'no_test', "
        err_msg += "'warn', or 'raise_error'."
        raise ValueError(err_msg)

    # Now need to create the file for SCIENCE_AVAILABLE
    sci_avlble_dict = segments.segmentlistdict()
    # NOTE: Should this be overrideable in the config file?
    sci_avlble_name = "SCIENCE_AVAILABLE"
    for ifo in science_segs.keys():
        sci_avlble_dict[ifo + ":" + sci_avlble_name] = science_segs[ifo]

    sci_avlble_file = SegFile.from_segment_list_dict(
        "SCIENCE_AVAILABLE",
        sci_avlble_dict,
        ifo_list=science_segs.keys(),
        valid_segment=workflow.analysis_time,
        extension=".xml",
        tags=tags,
        directory=output_dir,
    )

    logger.info("Leaving datafind module")
    if datafind_method == "AT_RUNTIME_FAKE_DATA":
        datafindouts = None
    else:
        datafindouts = FileList(datafindouts)

    return datafindouts, sci_avlble_file, science_segs, sci_avlble_name


def setup_datafind_runtime_cache_multi_calls_perifo(
    cp, science_segs, output_dir, tags=None
):
    """
    This function uses the `gwdatafind` library to obtain the location of all
    the frame files that will be needed to cover the analysis of the data
    given in science_segs. This function will not check if the returned frames
    cover the whole time requested, such sanity checks are done in the
    pycbc.workflow.setup_datafind_workflow entry function. As opposed to
    setup_datafind_runtime_single_call_perifo this call will one call to the
    datafind server for every science segment. This function will return a list
    of output files that correspond to the cache .lcf files that are produced,
    which list the locations of all frame files. This will cause problems with
    pegasus, which expects to know about all input files (ie. the frame files
    themselves.)

    Parameters
    ----------
    cp : ConfigParser.ConfigParser instance
        This contains a representation of the information stored within the
        workflow configuration files
    science_segs : Dictionary of ifo keyed igwn_segments.segmentlist instances
        This contains the times that the workflow is expected to analyse.
    output_dir : path
        All output files written by datafind processes will be written to this
        directory.
    tags : list of strings, optional (default=None)
        Use this to specify tags. This can be used if this module is being
        called more than once to give call specific configuration (by setting
        options in [workflow-datafind-${TAG}] rather than [workflow-datafind]).
        This is also used to tag the Files returned by the class to uniqueify
        the Files and uniqueify the actual filename.
        FIXME: Filenames may not be unique with current codes!

    Returns
    -------
    datafindcaches : list of glue.lal.Cache instances
       The glue.lal.Cache representations of the various calls to the datafind
       server and the returned frame files.
    datafindOuts : pycbc.workflow.core.FileList
        List of all the datafind output files for use later in the pipeline.

    """
    if tags is None:
        tags = []

    # Now ready to loop over the input segments
    datafindouts = []
    datafindcaches = []
    logger.info("Querying datafind server for all science segments.")
    for ifo, science_segs_ifo in science_segs.items():
        observatory = ifo[0].upper()
        frame_type = cp.get_opt_tags(
            "workflow-datafind", f"datafind-{ifo.lower()}-frame-type", tags
        )
        for seg in science_segs_ifo:
            msg = f"Finding data between {int(seg[0])} and {int(seg[1])} "
            msg += f"for ifo {ifo}"
            logger.info(msg)
            # WARNING: For now the workflow will expect times to be in integer seconds
            start_time = int(seg[0])
            end_time = int(seg[1])

            # Sometimes the connection can drop, so try a backup here
            try:
                cache, cache_file = run_datafind_instance(
                    cp,
                    output_dir,
                    observatory,
                    frame_type,
                    start_time,
                    end_time,
                    ifo,
                    tags=tags,
                )
            except:
                cache, cache_file = run_datafind_instance(
                    cp,
                    output_dir,
                    observatory,
                    frame_type,
                    start_time,
                    end_time,
                    ifo,
                    tags=tags,
                )
            datafindouts.append(cache_file)
            datafindcaches.append(cache)
    return datafindcaches, datafindouts


def setup_datafind_runtime_cache_single_call_perifo(
    cp, science_segs, output_dir, tags=None
):
    """
    This function uses the `gwdatafind` library to obtain the location of all
    the frame files that will be needed to cover the analysis of the data
    given in science_segs. This function will not check if the returned frames
    cover the whole time requested, such sanity checks are done in the
    pycbc.workflow.setup_datafind_workflow entry function. As opposed to
    setup_datafind_runtime_generated this call will only run one call to
    datafind per ifo, spanning the whole time. This function will return a list
    of output files that correspond to the cache .lcf files that are produced,
    which list the locations of all frame files. This will cause problems with
    pegasus, which expects to know about all input files (ie. the frame files
    themselves.)

    Parameters
    ----------
    cp : ConfigParser.ConfigParser instance
        This contains a representation of the information stored within the
        workflow configuration files
    science_segs : Dictionary of ifo keyed igwn_segments.segmentlist instances
        This contains the times that the workflow is expected to analyse.
    output_dir : path
        All output files written by datafind processes will be written to this
        directory.
    tags : list of strings, optional (default=None)
        Use this to specify tags. This can be used if this module is being
        called more than once to give call specific configuration (by setting
        options in [workflow-datafind-${TAG}] rather than [workflow-datafind]).
        This is also used to tag the Files returned by the class to uniqueify
        the Files and uniqueify the actual filename.
        FIXME: Filenames may not be unique with current codes!

    Returns
    -------
    datafindcaches : list of glue.lal.Cache instances
       The glue.lal.Cache representations of the various calls to the datafind
       server and the returned frame files.
    datafindOuts : pycbc.workflow.core.FileList
        List of all the datafind output files for use later in the pipeline.

    """
    if tags is None:
        tags = []

    # We want to ignore gaps as the detectors go up and down and calling this
    # way will give gaps. See the setup_datafind_runtime_generated function
    # for datafind calls that only query for data that will exist
    cp.set("datafind", "on_gaps", "ignore")

    # Now ready to loop over the input segments
    datafindouts = []
    datafindcaches = []
    logger.info("Querying datafind server for all science segments.")
    for ifo, science_segs_ifo in science_segs.items():
        observatory = ifo[0].upper()
        checked_times = segments.segmentlist([])
        frame_types = cp.get_opt_tags(
            "workflow-datafind", f"datafind-{ifo.lower()}-frame-type", tags
        )
        # Check if this is one type, or time varying
        frame_types = frame_types.replace(" ", "").strip().split(",")
        for ftype in frame_types:
            # Check the times, default to full time initially
            # This REQUIRES a coalesced segment list to work
            start = int(science_segs_ifo[0][0])
            end = int(science_segs_ifo[-1][1])
            # Then check for limits. We're expecting something like:
            # value[start:end], so need to extract value, start and end
            if "[" in ftype:
                # This gets start and end out
                bopt = ftype.split("[")[1].split("]")[0]
                newstart, newend = bopt.split(":")
                # Then check if the times are within science time
                start = max(int(newstart), start)
                end = min(int(newend), end)
                if end <= start:
                    continue
                # This extracts value
                ftype = ftype.split("[")[0]
            curr_times = segments.segment(start, end)
            # The times here must be distinct. We cannot have two different
            # frame files at the same time from the same ifo.
            if checked_times.intersects_segment(curr_times):
                err_msg = "Different frame types cannot overlap in time."
                raise ValueError(err_msg)
            checked_times.append(curr_times)

            # Ask datafind where the frames are
            try:
                cache, cache_file = run_datafind_instance(
                    cp, output_dir, observatory, ftype, start, end, ifo, tags=tags
                )
            except:
                cache, cache_file = run_datafind_instance(
                    cp, output_dir, observatory, ftype, start, end, ifo, tags=tags
                )

            datafindouts.append(cache_file)
            datafindcaches.append(cache)
    return datafindcaches, datafindouts


def setup_datafind_runtime_frames_single_call_perifo(
    cp, science_segs, output_dir, tags=None
):
    """
    This function uses the `gwdatafind` library to obtain the location of all
    the frame files that will be needed to cover the analysis of the data
    given in science_segs. This function will not check if the returned frames
    cover the whole time requested, such sanity checks are done in the
    pycbc.workflow.setup_datafind_workflow entry function. As opposed to
    setup_datafind_runtime_generated this call will only run one call to
    datafind per ifo, spanning the whole time. This function will return a list
    of files corresponding to the individual frames returned by the datafind
    query. This will allow pegasus to more easily identify all the files used
    as input, but may cause problems for codes that need to take frame cache
    files as input.

    Parameters
    ----------
    cp : ConfigParser.ConfigParser instance
        This contains a representation of the information stored within the
        workflow configuration files
    science_segs : Dictionary of ifo keyed igwn_segments.segmentlist instances
        This contains the times that the workflow is expected to analyse.
    output_dir : path
        All output files written by datafind processes will be written to this
        directory.
    tags : list of strings, optional (default=None)
        Use this to specify tags. This can be used if this module is being
        called more than once to give call specific configuration (by setting
        options in [workflow-datafind-${TAG}] rather than [workflow-datafind]).
        This is also used to tag the Files returned by the class to uniqueify
        the Files and uniqueify the actual filename.
        FIXME: Filenames may not be unique with current codes!

    Returns
    -------
    datafindcaches : list of glue.lal.Cache instances
       The glue.lal.Cache representations of the various calls to the datafind
       server and the returned frame files.
    datafindOuts : pycbc.workflow.core.FileList
        List of all the datafind output files for use later in the pipeline.

    """
    datafindcaches, _ = setup_datafind_runtime_cache_single_call_perifo(
        cp, science_segs, output_dir, tags=tags
    )

    datafindouts = convert_cachelist_to_filelist(datafindcaches)

    return datafindcaches, datafindouts


def setup_datafind_runtime_frames_multi_calls_perifo(
    cp, science_segs, output_dir, tags=None
):
    """
    This function uses the `gwdatafind` library to obtain the location of all
    the frame files that will be needed to cover the analysis of the data
    given in science_segs. This function will not check if the returned frames
    cover the whole time requested, such sanity checks are done in the
    pycbc.workflow.setup_datafind_workflow entry function. As opposed to
    setup_datafind_runtime_single_call_perifo this call will one call to the
    datafind server for every science segment. This function will return a list
    of files corresponding to the individual frames returned by the datafind
    query. This will allow pegasus to more easily identify all the files used
    as input, but may cause problems for codes that need to take frame cache
    files as input.

    Parameters
    ----------
    cp : ConfigParser.ConfigParser instance
        This contains a representation of the information stored within the
        workflow configuration files
    science_segs : Dictionary of ifo keyed igwn_segments.segmentlist instances
        This contains the times that the workflow is expected to analyse.
    output_dir : path
        All output files written by datafind processes will be written to this
        directory.
    tags : list of strings, optional (default=None)
        Use this to specify tags. This can be used if this module is being
        called more than once to give call specific configuration (by setting
        options in [workflow-datafind-${TAG}] rather than [workflow-datafind]).
        This is also used to tag the Files returned by the class to uniqueify
        the Files and uniqueify the actual filename.
        FIXME: Filenames may not be unique with current codes!

    Returns
    -------
    datafindcaches : list of glue.lal.Cache instances
       The glue.lal.Cache representations of the various calls to the datafind
       server and the returned frame files.
    datafindOuts : pycbc.workflow.core.FileList
        List of all the datafind output files for use later in the pipeline.

    """
    datafindcaches, _ = setup_datafind_runtime_cache_multi_calls_perifo(
        cp, science_segs, output_dir, tags=tags
    )

    datafindouts = convert_cachelist_to_filelist(datafindcaches)

    return datafindcaches, datafindouts


def setup_datafind_from_pregenerated_lcf_files(cp, ifos, output_dir, tags=None):
    """
    This function is used if you want to run with pregenerated lcf frame
    cache files.

    Parameters
    ----------
    cp : ConfigParser.ConfigParser instance
        This contains a representation of the information stored within the
        workflow configuration files
    ifos : list of ifo strings
        List of ifos to get pregenerated files for.
    output_dir : path
        All output files written by datafind processes will be written to this
        directory. Currently this sub-module writes no output.
    tags : list of strings, optional (default=None)
        Use this to specify tags. This can be used if this module is being
        called more than once to give call specific configuration (by setting
        options in [workflow-datafind-${TAG}] rather than [workflow-datafind]).
        This is also used to tag the Files returned by the class to uniqueify
        the Files and uniqueify the actual filename.

    Returns
    -------
    datafindcaches : list of glue.lal.Cache instances
       The glue.lal.Cache representations of the various calls to the datafind
       server and the returned frame files.
    datafindOuts : pycbc.workflow.core.FileList
        List of all the datafind output files for use later in the pipeline.

    """
    from glue import lal

    if tags is None:
        tags = []

    datafindcaches = []
    for ifo in ifos:
        search_string = f"datafind-pregenerated-cache-file-{ifo.lower()}"
        frame_cache_file_name = cp.get_opt_tags(
            "workflow-datafind", search_string, tags=tags
        )
        curr_cache = lal.Cache.fromfilenames(
            [frame_cache_file_name], coltype=lal.LIGOTimeGPS
        )
        curr_cache.ifo = ifo
        datafindcaches.append(curr_cache)
    datafindouts = convert_cachelist_to_filelist(datafindcaches)

    return datafindcaches, datafindouts


def convert_cachelist_to_filelist(datafindcache_list):
    """
    Take as input a list of glue.lal.Cache objects and return a pycbc FileList
    containing all frames within those caches.

    Parameters
    ----------
    datafindcache_list : list of glue.lal.Cache objects
        The list of cache files to convert.

    Returns
    -------
    datafind_filelist : FileList of frame File objects
        The list of frame files.

    """
    prev_file = None
    prev_name = None
    this_name = None

    datafind_filelist = FileList([])

    for cache in datafindcache_list:
        # sort the cache into time sequential order
        cache.sort()
        curr_ifo = cache.ifo
        for frame in cache:
            # Pegasus doesn't like "localhost" in URLs.
            frame.url = frame.url.replace("file://localhost", "file://")
            # Not sure why it happens in OSDF URLs!!
            # May need to remove use of Cache objects
            frame.url = frame.url.replace("osdf://localhost/", "osdf:///")

            # Create one File() object for each unique frame file that we
            # get back in the cache.
            if prev_file:
                prev_name = os.path.basename(prev_file.cache_entry.url)
                this_name = os.path.basename(frame.url)

            if (prev_file is None) or (prev_name != this_name):
                curr_file = File(
                    curr_ifo,
                    frame.description,
                    frame.segment,
                    file_url=frame.url,
                    use_tmp_subdirs=True,
                )
                datafind_filelist.append(curr_file)
                prev_file = curr_file

            # Populate the PFNs for the File() we just created
            cvmfs_urls = ("file:///cvmfs/", "osdf://")
            if frame.url.startswith(cvmfs_urls):
                # Frame is on CVMFS/OSDF, so let all sites read it directly.
                curr_file.add_pfn(frame.url, site="all")
            elif frame.url.startswith("file://"):
                # Frame not on CVMFS, so may need transferring.
                # Be careful here! If all your frames files are on site
                # = local and you try to run on OSG, it will likely
                # overwhelm the condor file transfer process!
                for site in ["local", "condorpool_shared"]:
                    curr_file.add_pfn(frame.url, site=site)
            else:
                # Frame is at some unknown URL. Pegasus will decide how to deal
                # with this, but will likely transfer to local site first, and
                # from there transfer to remote sites as needed.
                curr_file.add_pfn(frame.url, site="notlocal")

    return datafind_filelist


def get_science_segs_from_datafind_outs(datafindcaches):
    """
    This function will calculate the science segments that are covered in
    the OutGroupList containing the frame files returned by various
    calls to the datafind server. This can then be used to check whether this
    list covers what it is expected to cover.

    Parameters
    ----------
    datafindcaches : OutGroupList
        List of all the datafind output files.

    Returns
    -------
    new_science_segs : Dictionary of ifo keyed igwn_segments.segmentlist instances
        The times covered by the frames found in datafindOuts.

    """
    new_science_segs = {}
    for cache in datafindcaches:
        if len(cache) > 0:
            group_segs = segments.segmentlist(e.segment for e in cache).coalesce()
            ifo = cache.ifo
            if ifo not in new_science_segs:
                new_science_segs[ifo] = group_segs
            else:
                new_science_segs[ifo].extend(group_segs)
                new_science_segs[ifo].coalesce()
    return new_science_segs


def get_missing_segs_from_frame_file_cache(datafindcaches):
    """
    This function will use os.path.isfile to determine if all the frame files
    returned by the local datafind server actually exist on the disk. This can
    then be used to update the science times if needed.

    Parameters
    ----------
    datafindcaches : OutGroupList
        List of all the datafind output files.

    Returns
    -------
    missing_frame_segs : Dict. of ifo keyed igwn_segments.segmentlist instances
        The times corresponding to missing frames found in datafindOuts.
    missing_frames: Dict. of ifo keyed lal.Cache instances
        The list of missing frames

    """
    from glue import lal

    missing_frame_segs = {}
    missing_frames = {}
    for cache in datafindcaches:
        if len(cache) > 0:
            # Don't bother if these are not file:// urls, assume all urls in
            # one cache file must be the same type
            if not cache[0].scheme == "file":
                warn_msg = f"We have {cache[0].scheme} entries in the "
                warn_msg += "cache file. I do not check if these exist."
                logger.warning(warn_msg)
                continue
            _, curr_missing_frames = cache.checkfilesexist(on_missing="warn")
            missing_segs = segments.segmentlist(
                e.segment for e in curr_missing_frames
            ).coalesce()
            ifo = cache.ifo
            if ifo not in missing_frame_segs:
                missing_frame_segs[ifo] = missing_segs
                missing_frames[ifo] = lal.Cache(curr_missing_frames)
            else:
                missing_frame_segs[ifo].extend(missing_segs)
                # NOTE: This .coalesce probably isn't needed as the segments
                # should be disjoint. If speed becomes an issue maybe remove it?
                missing_frame_segs[ifo].coalesce()
                missing_frames[ifo].extend(curr_missing_frames)
    return missing_frame_segs, missing_frames


def get_segment_summary_times(science_file, segment_name):
    """
    This function will find the times for which the segment_summary is set
    for the flag given by segment_name.

    Parameters
    ----------
    science_file : SegFile
        The segment file that we want to use to determine this.
    segment_name : string
        The DQ flag to search for times in the segment_summary table.

    Returns
    -------
    summ_seg_list : igwn_segments.segmentlist
        The times that are covered in the segment summary table.

    """
    # Parse the segment_name
    segment_name = segment_name.split(":")
    if len(segment_name) not in [2, 3]:
        raise ValueError(f"Invalid channel name {segment_name}.")
    ifo = segment_name[0]
    channel = segment_name[1]
    version = ""
    if len(segment_name) == 3:
        version = int(segment_name[2])

    # Load the filename
    xmldoc = utils.load_filename(
        science_file.cache_entry.path,
        compress="auto",
        contenthandler=LIGOLWContentHandler,
    )

    # Get the segment_def_id for the segment_name
    segment_def_table = ligolw.Table.get_table(xmldoc, "segment_definer")
    for entry in segment_def_table:
        if (entry.ifos == ifo) and (entry.name == channel):
            if len(segment_name) == 2 or (entry.version == version):
                seg_def_id = entry.segment_def_id
                break
    else:
        raise ValueError(
            f"Cannot find channel {segment_name} in segment_definer table."
        )

    # Get the segmentlist corresponding to this segment_name in segment_summary
    segment_summ_table = ligolw.Table.get_table(xmldoc, "segment_summary")
    summ_seg_list = segments.segmentlist([])
    for entry in segment_summ_table:
        if entry.segment_def_id == seg_def_id:
            segment = segments.segment(entry.start_time, entry.end_time)
            summ_seg_list.append(segment)
    summ_seg_list.coalesce()

    return summ_seg_list


def run_datafind_instance(
    cp, output_dir, observatory, frame_type, start_time, end_time, ifo, tags=None
):
    """
    This function will query the datafind server once to find frames between
    the specified times for the specified frame type and observatory.

    Parameters
    ----------
    cp : ConfigParser instance
        Source for any kwargs that should be sent to the datafind module
    output_dir : Output cache files will be written here. We also write the
        commands for reproducing what is done in this function to this
        directory.
    observatory : string
        The observatory to query frames for. Ex. 'H', 'L' or 'V'.  NB: not
        'H1', 'L1', 'V1' which denote interferometers.
    frame_type : string
        The frame type to query for.
    start_time : int
        Integer start time to query the datafind server for frames.
    end_time : int
        Integer end time to query the datafind server for frames.
    ifo : string
        The interferometer to use for naming output. Ex. 'H1', 'L1', 'V1'.
        Maybe this could be merged with the observatory string, but this
        could cause issues if running on old 'H2' and 'H1' data.
    tags : list of string, optional (default=None)
        Use this to specify tags. This can be used if this module is being
        called more than once to give call specific configuration (by setting
        options in [workflow-datafind-${TAG}] rather than [workflow-datafind]).
        This is also used to tag the Files returned by the class to uniqueify
        the Files and uniquify the actual filename.
        FIXME: Filenames may not be unique with current codes!

    Returns
    -------
    df_cache : glue.lal.Cache instance
       The glue.lal.Cache representation of the call to the datafind
       server and the returned frame files.
    cacheFile : pycbc.workflow.core.File
        Cache file listing all of the datafind output files for use later in the pipeline.

    """
    from glue import lal

    if tags is None:
        tags = []

    # Determine if we should override the default datafind server
    if cp.has_option_tags("workflow-datafind", "datafind-ligo-datafind-server", tags):
        datafind_server = cp.get_opt_tags(
            "workflow-datafind", "datafind-ligo-datafind-server", tags
        )
    else:
        datafind_server = None

    seg = segments.segment([start_time, end_time])
    # Take the datafind kwargs from config (usually urltype=file is
    # given).
    df_kwargs = {}
    # By default ignore missing frames, this case is dealt with outside of here
    df_kwargs["on_gaps"] = "ignore"
    if cp.has_section("datafind"):
        for item, value in cp.items("datafind"):
            df_kwargs[item] = value
    for tag in tags:
        if cp.has_section(f"datafind-{tag}"):
            for item, value in cp.items(f"datafind-{tag}"):
                df_kwargs[item] = value

    # It is useful to print the corresponding command to the logs
    # directory to check if this was expected.
    log_datafind_command(
        observatory,
        frame_type,
        start_time,
        end_time,
        os.path.join(output_dir, "logs"),
        **df_kwargs,
    )
    logger.debug("Asking datafind server for frames.")
    df_cache = lal.Cache.from_urls(
        find_frame_urls(
            observatory,
            frame_type,
            start_time,
            end_time,
            host=datafind_server,
            **df_kwargs,
        ),
    )
    logger.debug("Frames returned")
    # workflow format output file
    cache_file = File(
        ifo, "DATAFIND", seg, extension="lcf", directory=output_dir, tags=tags
    )
    cache_file.add_pfn(cache_file.cache_entry.path, site="local")

    df_cache.ifo = ifo
    # Dump output to file
    fp = open(cache_file.storage_path, "w")
    # FIXME: CANNOT use df_cache.tofile because it will print 815901601.00000
    #        as a gps time which is incompatible with the lal cache format
    #        (and the C codes) which demand an integer.
    # df_cache.tofile(fp)
    for entry in df_cache:
        start = str(int(entry.segment[0]))
        duration = str(int(abs(entry.segment)))
        print(
            f"{entry.observatory} {entry.description} {start} {duration} {entry.url}",
            file=fp,
        )
        entry.segment = segments.segment(int(entry.segment[0]), int(entry.segment[1]))

    fp.close()
    return df_cache, cache_file


def log_datafind_command(
    observatory, frame_type, start_time, end_time, output_dir, **df_kwargs
):
    """
    This command will print an equivalent gw_data_find command to disk that
    can be used to debug why the internal datafind module is not working.
    """
    # FIXME: This does not accurately reproduce the call as assuming the
    # kwargs will be the same is wrong, so some things need to be converted
    # "properly" to the command line equivalent.
    gw_command = [
        "gw_data_find",
        "--observatory",
        observatory,
        "--type",
        frame_type,
        "--gps-start-time",
        str(start_time),
        "--gps-end-time",
        str(end_time),
    ]

    for name, value in df_kwargs.items():
        if name == "match":
            gw_command.append("--match")
            gw_command.append(str(value))
        elif name == "urltype":
            gw_command.append("--url-type")
            gw_command.append(str(value))
        elif name == "on_gaps":
            pass
        else:
            err_msg = f"Unknown datafind kwarg given: {name}. "
            err_msg += "This argument is stripped in the logged .sh command."
            logger.warning(err_msg)

    file_name = (
        f"{observatory}-{frame_type}-{int(start_time)}-{int(end_time - start_time)}.sh"
    )
    file_path = os.path.join(output_dir, file_name)
    fp = open(file_path, "w")
    fp.write(" ".join(gw_command))
    fp.close()


def datafind_keep_unique_backups(backup_outs, orig_outs):
    """
    This function will take a list of backup datafind files, presumably
    obtained by querying a remote datafind server, e.g. CIT, and compares
    these against a list of original datafind files, presumably obtained by
    querying the local datafind server. Only the datafind files in the backup
    list that do not appear in the original list are returned. This allows us
    to use only files that are missing from the local cluster.

    Parameters
    ----------
    backup_outs : FileList
        List of datafind files from the remote datafind server.
    orig_outs : FileList
        List of datafind files from the local datafind server.

    Returns
    -------
    FileList
        List of datafind files in backup_outs and not in orig_outs.

    """
    # NOTE: This function is not optimized and could be made considerably
    #       quicker if speed becomes in issue. With 4s frame files this might
    #       be slow, but for >1000s files I don't foresee any issue, so I keep
    #       this simple.
    return_list = FileList([])
    # We compare the LFNs to determine uniqueness
    # Is there a way to associate two paths with one LFN??
    orig_names = [f.name for f in orig_outs]
    for file in backup_outs:
        if file.name not in orig_names:
            return_list.append(file)
        else:
            index_num = orig_names.index(file.name)
            orig_out = orig_outs[index_num]
            pfns = list(file.pfns)
            # This shouldn't happen, but catch if it does
            assert len(pfns) == 1
            orig_out.add_pfn(pfns[0].url, site="notlocal")

    return return_list
