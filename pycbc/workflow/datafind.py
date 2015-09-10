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

import os, copy
import urlparse
import logging
from glue import segments, lal
from glue.ligolw import utils, table, lsctables, ligolw
from pycbc.workflow.core import OutSegFile, File, FileList, make_analysis_dir
from pycbc.frame import datafind_connection

class ContentHandler(ligolw.LIGOLWContentHandler):
        pass

lsctables.use_in(ContentHandler)

def setup_datafind_workflow(workflow, scienceSegs,  outputDir, segFilesList,
                            tag=None):
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
    scienceSegs : Dictionary of ifo keyed glue.segment.segmentlist instances
        This contains the times that the workflow is expected to analyse.
    outputDir : path
        All output files written by datafind processes will be written to this
        directory.
    segFilesList : List of the files returned by segment_utils
        This contains representations of the various segment files that were
        constructed at the segment generation stage of the workflow. This will
        be used for the segment_summary test, or if any of the other tests are
        given "update_times" (and can be given a value of None otherwise).
    tag : string, optional (default=None)
        Use this to specify a tag. This can be used if this module is being
        called more than once to give call specific configuration (by setting
        options in [workflow-datafind-${TAG}] rather than [workflow-datafind]). 
        This is also used to tag the Files returned by the class to uniqueify
        the Files and uniqueify the actual filename.
        FIXME: Filenames may not be unique with current codes!

    Returns
    --------
    datafindOuts : OutGroupList
        List of all the datafind output files for use later in the pipeline.
    scienceSegs : Dictionary of ifo keyed glue.segment.segmentlist instances
        This contains the times that the workflow is expected to analyse. If 
        the updateSegmentTimes kwarg is given this will be updated to reflect 
        any instances of missing data.
    """
    logging.info("Entering datafind module")
    make_analysis_dir(outputDir)
    cp = workflow.cp

    # Parse for options in ini file
    datafindMethod = cp.get_opt_tags("workflow-datafind",
                                     "datafind-method", [tag])

    if cp.has_option_tags("workflow-datafind",
                          "datafind-check-segment-gaps", [tag]):
        checkSegmentGaps = cp.get_opt_tags("workflow-datafind", 
                                          "datafind-check-segment-gaps", [tag])
    else:
        checkSegmentGaps = "no_test"
    if cp.has_option_tags("workflow-datafind",
                          "datafind-check-frames-exist", [tag]):
        checkFramesExist = cp.get_opt_tags("workflow-datafind",
                                          "datafind-check-frames-exist", [tag])
    else:
        checkFramesExist = "no_test"
    if cp.has_option_tags("workflow-datafind",
                          "datafind-check-segment-summary", [tag]):
        checkSegmentSummary = cp.get_opt_tags("workflow-datafind",
                                       "datafind-check-segment-summary", [tag])
    else:
        checkSegmentSummary = "no_test"
    
    logging.info("Starting datafind with setup_datafind_runtime_generated")
    if datafindMethod == "AT_RUNTIME_MULTIPLE_CACHES":
        datafindcaches, datafindouts = \
            setup_datafind_runtime_cache_multi_calls_perifo(cp, scienceSegs,
                                                            outputDir, tag=tag)
    elif datafindMethod == "AT_RUNTIME_SINGLE_CACHES":
        datafindcaches, datafindouts = \
            setup_datafind_runtime_cache_single_call_perifo(cp, scienceSegs, 
                                                            outputDir, tag=tag)
    elif datafindMethod == "AT_RUNTIME_MULTIPLE_FRAMES":
        datafindcaches, datafindouts = \
            setup_datafind_runtime_frames_multi_calls_perifo(cp, scienceSegs,
                                                            outputDir, tag=tag)
    elif datafindMethod == "AT_RUNTIME_SINGLE_FRAMES":
        datafindcaches, datafindouts = \
            setup_datafind_runtime_frames_single_call_perifo(cp, scienceSegs,
                                                            outputDir, tag=tag)

    elif datafindMethod == "FROM_PREGENERATED_LCF_FILES":
        ifos = scienceSegs.keys()
        datafindcaches, datafindouts = \
            setup_datafind_from_pregenerated_lcf_files(cp, ifos,
                                                       outputDir, tag=tag)
    else:
        msg = "Entry datafind-method in [workflow-datafind] does not have "
        msg += "expected value. Valid values are "
        msg += "AT_RUNTIME_MULTIPLE_FRAMES, AT_RUNTIME_SINGLE_FRAMES "
        msg += "AT_RUNTIME_MULTIPLE_CACHES or AT_RUNTIME_SINGLE_CACHES. "
        msg += "Consult the documentation for more info."
        raise ValueError(msg)

    using_backup_server = False
    if datafindMethod == "AT_RUNTIME_MULTIPLE_FRAMES" or \
                                  datafindMethod == "AT_RUNTIME_SINGLE_FRAMES":
        if cp.has_option_tags("workflow-datafind",
                          "datafind-backup-datafind-server", [tag]):
            using_backup_server = True
            backup_server = cp.get_opt_tags("workflow-datafind",
                                      "datafind-backup-datafind-server", [tag])
            cp_new = copy.deepcopy(cp)
            cp_new.set("workflow-datafind",
                                "datafind-ligo-datafind-server", backup_server)
            cp_new.set('datafind', 'urltype', 'gsiftp')
            backup_datafindcaches, backup_datafindouts =\
                setup_datafind_runtime_frames_single_call_perifo(cp_new,
                                               scienceSegs, outputDir, tag=tag)
            backup_datafindouts = datafind_keep_unique_backups(\
                                             backup_datafindouts, datafindouts)
            datafindcaches.extend(backup_datafindcaches)
            datafindouts.extend(backup_datafindouts)

    logging.info("setup_datafind_runtime_generated completed")
    # If we don't have frame files covering all times we can update the science
    # segments.
    if checkSegmentGaps in ['warn','update_times','raise_error']:
        logging.info("Checking science segments against datafind output....")
        newScienceSegs = get_science_segs_from_datafind_outs(datafindcaches)
        logging.info("Datafind segments calculated.....")
        missingData = False
        msg = "Any errors directly following this message refer to times that"
        msg += " the segment server says are science, but datafind cannot find"
        msg += "frames for:"
        logging.info(msg)
        for ifo in scienceSegs.keys():
            # If no data in the input then do nothing
            if not scienceSegs[ifo]:
                msg = "No input science segments for ifo %s " %(ifo)
                msg += "so, surprisingly, no data has been found. "
                msg += "Was this expected?"
                logging.warning(msg)
                continue
            if not newScienceSegs.has_key(ifo):
                msg = "IFO %s's science segments " %(ifo)
                msg += "are completely missing."
                logging.error(msg)
                missingData = True
                if checkSegmentGaps == 'update_times':
                    scienceSegs[ifo] = segments.segmentlist()
                continue
            missing = scienceSegs[ifo] - newScienceSegs[ifo]
            if abs(missing):
                msg = "From ifo %s we are missing frames covering:" %(ifo)
                msg += "\n%s" % "\n".join(map(str, missing))
                missingData = True
                logging.error(msg)
                if checkSegmentGaps == 'update_times':
                    # Remove missing time, so that we can carry on if desired
                    logging.info("Updating science times for ifo %s." %(ifo))
                    scienceSegs[ifo] = scienceSegs[ifo] - missing

        if checkSegmentGaps == 'raise_error' and missingData:
            raise ValueError("Workflow cannot find needed data, exiting.")
        logging.info("Done checking, any discrepancies are reported above.")
    elif checkSegmentGaps == 'no_test':
        # Do nothing
        pass
    else:
        errMsg = "checkSegmentGaps kwArg must take a value from 'no_test', "
        errMsg += "'warn', 'update_times' or 'raise_error'."
        raise ValueError(errMsg)

    # Do all of the frame files that were returned actually exist?
    if checkFramesExist in ['warn','update_times','raise_error']:
        logging.info("Verifying that all frames exist on disk.")
        missingFrSegs, missingFrames = \
                          get_missing_segs_from_frame_file_cache(datafindcaches)
        missingFlag = False
        for ifo in missingFrames.keys():
            # If no data in the input then do nothing
            if not scienceSegs[ifo]:
                continue
            # If using a backup server, does the frame exist remotely?
            if using_backup_server:
                # WARNING: This will be slow, but hopefully it will not occur
                #          for too many frames. This could be optimized if
                #          it becomes necessary.
                new_list = []
                for frame in missingFrames[ifo]:
                    for dfout in datafindouts:
                        dfout_pfns = list(dfout.pfns)
                        dfout_urls = [a.url for a in dfout_pfns]
                        if frame.url in dfout_urls:
                            pfn = dfout_pfns[dfout_urls.index(frame.url)]
                            dfout.removePFN(pfn)
                            if len(dfout.pfns) == 0:
                                new_list.append(frame)
                            else:
                                msg = "Frame %s not found locally. "\
                                                                  %(frame.url,)
                                msg += "Replacing with remote url(s) "
                                msg += "%s." \
                                           %(str([a.url for a in dfout.pfns]),)
                                logging.info(msg)
                            break
                    else:
                        new_list.append(frame)
                missingFrames[ifo] = new_list
            if missingFrames[ifo]:
                msg = "From ifo %s we are missing the following frames:" %(ifo)
                msg +='\n'.join([a.url for a in missingFrames[ifo]])
                missingFlag = True
                logging.error(msg)
            if checkFramesExist == 'update_times':
                # Remove missing times, so that we can carry on if desired
                logging.info("Updating science times for ifo %s." %(ifo))
                scienceSegs[ifo] = scienceSegs[ifo] - missingFrSegs[ifo]
                
        if checkFramesExist == 'raise_error' and missingFlag:
            raise ValueError("Workflow cannot find all frames, exiting.")
        logging.info("Finished checking frames.")
    elif checkFramesExist == 'no_test':
        # Do nothing
        pass
    else:
        errMsg = "checkFramesExist kwArg must take a value from 'no_test', "
        errMsg += "'warn', 'update_times' or 'raise_error'."
        raise ValueError(errMsg)

    # Check if there are cases where frames exist, but no entry in the segment
    # summary table are present.
    if checkSegmentSummary in ['warn', 'raise_error']:
        logging.info("Checking the segment summary table against frames.")
        dfScienceSegs = get_science_segs_from_datafind_outs(datafindcaches)
        missingFlag = False
        for ifo in dfScienceSegs.keys():
            scienceFile = segFilesList.find_output_with_ifo(ifo)
            scienceFile = scienceFile.find_output_with_tag('SCIENCE')
            if not len(scienceFile) == 1:
                errMsg = "Did not find exactly 1 science file."
                raise ValueError(errMsg)
            scienceFile = scienceFile[0]

            scienceChannel = cp.get('workflow-segments',\
                                'segments-%s-science-name'%(ifo.lower()))
            segSummaryTimes = get_segment_summary_times(scienceFile,
                                                        scienceChannel)
            missing = dfScienceSegs[ifo] - segSummaryTimes
            scienceButNotFrame = scienceSegs[ifo] - dfScienceSegs[ifo]
            missing2 = scienceSegs[ifo] - scienceButNotFrame
            missing2 = missing2 - segSummaryTimes
            if abs(missing):
                msg = "From ifo %s the following times have frames, " %(ifo)
                msg += "but are not covered in the segment summary table."
                msg += "\n%s" % "\n".join(map(str, missing))
                logging.error(msg)
                missingFlag = True
            if abs(missing2):
                msg = "From ifo %s the following times have frames, " %(ifo)
                msg += "are science, and are not covered in the segment "
                msg += "summary table."
                msg += "\n%s" % "\n".join(map(str, missing2))
                logging.error(msg)
                missingFlag = True
        if checkSegmentSummary == 'raise_error' and missingFlag:
            errMsg = "Segment_summary discrepancy detected, exiting."
            raise ValueError(errMsg)
    elif checkSegmentSummary == 'no_test':
        # Do nothing
        pass
    else:
        errMsg = "checkSegmentSummary kwArg must take a value from 'no_test', "
        errMsg += "'warn', or 'raise_error'."
        raise ValueError(errMsg)

    # Now need to create the file for SCIENCE_AVAILABLE
    for ifo in scienceSegs.keys():
        availableSegsFile = os.path.abspath(os.path.join(outputDir, 
                           "%s-SCIENCE_AVAILABLE_SEGMENTS.xml" %(ifo.upper()) ))
        currUrl = urlparse.urlunparse(['file', 'localhost', availableSegsFile,
                          None, None, None])
        if tag:
            currTags = [tag, 'SCIENCE_AVAILABLE']
        else:
            currTags = ['SCIENCE_AVAILABLE']
        currFile = OutSegFile(ifo, 'SEGMENTS', workflow.analysis_time,
                            currUrl, segment_list=scienceSegs[ifo], tags = currTags)
        currFile.PFN(availableSegsFile, site='local')
        segFilesList.append(currFile)
        currFile.toSegmentXml()
   

    logging.info("Leaving datafind module")
    return FileList(datafindouts), scienceSegs
    

def setup_datafind_runtime_cache_multi_calls_perifo(cp, scienceSegs, 
                                                    outputDir, tag=None):
    """
    This function uses the glue.datafind library to obtain the location of all
    the frame files that will be needed to cover the analysis of the data
    given in scienceSegs. This function will not check if the returned frames
    cover the whole time requested, such sanity checks are done in the
    pycbc.workflow.setup_datafind_workflow entry function. As opposed to
    setup_datafind_runtime_single_call_perifo this call will one call to the
    datafind server for every science segment. This function will return a list
    of output files that correspond to the cache .lcf files that are produced,
    which list the locations of all frame files. This will cause problems with
    pegasus, which expects to know about all input files (ie. the frame files
    themselves.)

    Parameters
    -----------
    cp : ConfigParser.ConfigParser instance
        This contains a representation of the information stored within the
        workflow configuration files
    scienceSegs : Dictionary of ifo keyed glue.segment.segmentlist instances
        This contains the times that the workflow is expected to analyse.
    outputDir : path
        All output files written by datafind processes will be written to this
        directory.
    tag : string, optional (default=None)
        Use this to specify a tag. This can be used if this module is being
        called more than once to give call specific configuration (by setting
        options in [workflow-datafind-${TAG}] rather than [workflow-datafind]). This
        is also used to tag the Files returned by the class to uniqueify
        the Files and uniqueify the actual filename.
        FIXME: Filenames may not be unique with current codes!

    Returns
    --------
    datafindcaches : list of glue.lal.Cache instances
       The glue.lal.Cache representations of the various calls to the datafind
       server and the returned frame files.
    datafindOuts : pycbc.workflow.core.FileList
        List of all the datafind output files for use later in the pipeline.

    """
    # First job is to do setup for the datafind jobs
    # First get the server name
    logging.info("Setting up connection to datafind server.")
    connection = setup_datafind_server_connection(cp, tag=tag)

    # Now ready to loop over the input segments
    datafindouts = []
    datafindcaches = []
    ifos = scienceSegs.keys()
    logging.info("Querying datafind server for all science segments.")
    for ifo, scienceSegsIfo in scienceSegs.items():
        observatory = ifo[0].upper()
        frameType = cp.get_opt_tags("workflow-datafind", 
                                    "datafind-%s-frame-type" % (ifo.lower()), [tag])
        for seg in scienceSegsIfo:
            msg = "Finding data between %d and %d " %(seg[0],seg[1])
            msg += "for ifo %s" %(ifo)
            logging.debug(msg)
            # WARNING: For now the workflow will expect times to be in integer seconds
            startTime = int(seg[0])
            endTime = int(seg[1])

            # Sometimes the connection can drop, so try a backup here
            try:
                cache, cache_file = run_datafind_instance(cp, outputDir,
                                           connection, observatory, frameType,
                                           startTime, endTime, ifo, tag=tag)
            except:
                connection = setup_datafind_server_connection(cp, tag=tag)
                cache, cache_file = run_datafind_instance(cp, outputDir,
                                           connection, observatory, frameType,
                                           startTime, endTime, ifo, tag=tag)
            datafindouts.append(cache_file)
            datafindcaches.append(cache)
    return datafindcaches, datafindouts

def setup_datafind_runtime_cache_single_call_perifo(cp, scienceSegs, outputDir,
                                              tag=None):
    """
    This function uses the glue.datafind library to obtain the location of all
    the frame files that will be needed to cover the analysis of the data
    given in scienceSegs. This function will not check if the returned frames
    cover the whole time requested, such sanity checks are done in the
    pycbc.workflow.setup_datafind_workflow entry function. As opposed to 
    setup_datafind_runtime_generated this call will only run one call to
    datafind per ifo, spanning the whole time. This function will return a list
    of output files that correspond to the cache .lcf files that are produced,
    which list the locations of all frame files. This will cause problems with
    pegasus, which expects to know about all input files (ie. the frame files
    themselves.)

    Parameters
    -----------
    cp : ConfigParser.ConfigParser instance
        This contains a representation of the information stored within the
        workflow configuration files
    scienceSegs : Dictionary of ifo keyed glue.segment.segmentlist instances
        This contains the times that the workflow is expected to analyse.
    outputDir : path
        All output files written by datafind processes will be written to this
        directory.
    tag : string, optional (default=None)
        Use this to specify a tag. This can be used if this module is being
        called more than once to give call specific configuration (by setting
        options in [workflow-datafind-${TAG}] rather than [workflow-datafind]). This
        is also used to tag the Files returned by the class to uniqueify
        the Files and uniqueify the actual filename.
        FIXME: Filenames may not be unique with current codes!

    Returns
    --------
    datafindcaches : list of glue.lal.Cache instances
       The glue.lal.Cache representations of the various calls to the datafind
       server and the returned frame files.
    datafindOuts : pycbc.workflow.core.FileList
        List of all the datafind output files for use later in the pipeline.

    """
    # First job is to do setup for the datafind jobs
    # First get the server name
    logging.info("Setting up connection to datafind server.")
    connection = setup_datafind_server_connection(cp, tag=tag)

    # We want to ignore gaps as the detectors go up and down and calling this
    # way will give gaps. See the setup_datafind_runtime_generated function
    # for datafind calls that only query for data that will exist
    cp.set("datafind","on_gaps","ignore")

    # Now ready to loop over the input segments
    datafindouts = []
    datafindcaches = []
    ifos = scienceSegs.keys()
    logging.info("Querying datafind server for all science segments.")
    for ifo, scienceSegsIfo in scienceSegs.items():
        observatory = ifo[0].upper()
        frameType = cp.get_opt_tags("workflow-datafind",
                                    "datafind-%s-frame-type" % (ifo.lower()), [tag])
        # This REQUIRES a coalesced segment list to work
        startTime = int(scienceSegsIfo[0][0])
        endTime = int(scienceSegsIfo[-1][1])
        try:
            cache, cache_file = run_datafind_instance(cp, outputDir, connection,
                                       observatory, frameType, startTime,
                                       endTime, ifo, tag=tag)
        except:
            connection = setup_datafind_server_connection(cp, tag=tag)
            cache, cache_file = run_datafind_instance(cp, outputDir, connection,
                                       observatory, frameType, startTime,
                                       endTime, ifo, tag=tag)

        datafindouts.append(cache_file)
        datafindcaches.append(cache)
    return datafindcaches, datafindouts

def setup_datafind_runtime_frames_single_call_perifo(cp, scienceSegs,
                                              outputDir, tag=None):
    """
    This function uses the glue.datafind library to obtain the location of all
    the frame files that will be needed to cover the analysis of the data
    given in scienceSegs. This function will not check if the returned frames
    cover the whole time requested, such sanity checks are done in the
    pycbc.workflow.setup_datafind_workflow entry function. As opposed to 
    setup_datafind_runtime_generated this call will only run one call to
    datafind per ifo, spanning the whole time. This function will return a list
    of files corresponding to the individual frames returned by the datafind
    query. This will allow pegasus to more easily identify all the files used
    as input, but may cause problems for codes that need to take frame cache
    files as input.

    Parameters
    -----------
    cp : ConfigParser.ConfigParser instance
        This contains a representation of the information stored within the
        workflow configuration files
    scienceSegs : Dictionary of ifo keyed glue.segment.segmentlist instances
        This contains the times that the workflow is expected to analyse.
    outputDir : path
        All output files written by datafind processes will be written to this
        directory.
    tag : string, optional (default=None)
        Use this to specify a tag. This can be used if this module is being
        called more than once to give call specific configuration (by setting
        options in [workflow-datafind-${TAG}] rather than [workflow-datafind]). This
        is also used to tag the Files returned by the class to uniqueify
        the Files and uniqueify the actual filename.
        FIXME: Filenames may not be unique with current codes!

    Returns
    --------
    datafindcaches : list of glue.lal.Cache instances
       The glue.lal.Cache representations of the various calls to the datafind
       server and the returned frame files.
    datafindOuts : pycbc.workflow.core.FileList
        List of all the datafind output files for use later in the pipeline.

    """
    datafindcaches, _ = \
        setup_datafind_runtime_cache_single_call_perifo(cp, scienceSegs,
                                                        outputDir, tag=tag)

    datafindouts = convert_cachelist_to_filelist(datafindcaches)

    return datafindcaches, datafindouts

def setup_datafind_runtime_frames_multi_calls_perifo(cp, scienceSegs,
                                                     outputDir, tag=None):
    """
    This function uses the glue.datafind library to obtain the location of all
    the frame files that will be needed to cover the analysis of the data
    given in scienceSegs. This function will not check if the returned frames
    cover the whole time requested, such sanity checks are done in the
    pycbc.workflow.setup_datafind_workflow entry function. As opposed to
    setup_datafind_runtime_single_call_perifo this call will one call to the
    datafind server for every science segment. This function will return a list
    of files corresponding to the individual frames returned by the datafind
    query. This will allow pegasus to more easily identify all the files used
    as input, but may cause problems for codes that need to take frame cache
    files as input.

    Parameters
    -----------
    cp : ConfigParser.ConfigParser instance
        This contains a representation of the information stored within the
        workflow configuration files
    scienceSegs : Dictionary of ifo keyed glue.segment.segmentlist instances
        This contains the times that the workflow is expected to analyse.
    outputDir : path
        All output files written by datafind processes will be written to this
        directory.
    tag : string, optional (default=None)
        Use this to specify a tag. This can be used if this module is being
        called more than once to give call specific configuration (by setting
        options in [workflow-datafind-${TAG}] rather than [workflow-datafind]). This
        is also used to tag the Files returned by the class to uniqueify
        the Files and uniqueify the actual filename.
        FIXME: Filenames may not be unique with current codes!

    Returns
    --------
    datafindcaches : list of glue.lal.Cache instances
       The glue.lal.Cache representations of the various calls to the datafind
       server and the returned frame files.
    datafindOuts : pycbc.workflow.core.FileList
        List of all the datafind output files for use later in the pipeline.

    """
    datafindcaches, _ = \
        setup_datafind_runtime_cache_multi_calls_perifo(cp, scienceSegs,
                                                        outputDir, tag=tag)

    datafindouts = convert_cachelist_to_filelist(datafindcaches)

    return datafindcaches, datafindouts

def setup_datafind_from_pregenerated_lcf_files(cp, ifos, outputDir, tag=None):
    """
    This function is used if you want to run with pregenerated lcf frame
    cache files. 

    Parameters
    -----------
    cp : ConfigParser.ConfigParser instance
        This contains a representation of the information stored within the
        workflow configuration files
    ifos : list of ifo strings
        List of ifos to get pregenerated files for.
    outputDir : path
        All output files written by datafind processes will be written to this
        directory. Currently this sub-module writes no output.
    tag : string, optional (default=None)
        Use this to specify a tag. This can be used if this module is being
        called more than once to give call specific configuration (by setting
        options in [workflow-datafind-${TAG}] rather than [workflow-datafind]).
        This is also used to tag the Files returned by the class to uniqueify
        the Files and uniqueify the actual filename.

    Returns
    --------
    datafindcaches : list of glue.lal.Cache instances
       The glue.lal.Cache representations of the various calls to the datafind
       server and the returned frame files.
    datafindOuts : pycbc.workflow.core.FileList
        List of all the datafind output files for use later in the pipeline.
    """
    datafindcaches = []
    for ifo in ifos:
        search_string = "datafind-pregenerated-cache-file-%s" %(ifo.lower(),)
        frame_cache_file_name = cp.get_opt_tags("workflow-datafind",
                                                search_string, tags=[tag])
        curr_cache = lal.Cache.fromfilenames([frame_cache_file_name],
                                             coltype=lal.LIGOTimeGPS)
        curr_cache.ifo = ifo
        datafindcaches.append(curr_cache)
    datafindouts = convert_cachelist_to_filelist(datafindcaches)

    return datafindcaches, datafindouts

def convert_cachelist_to_filelist(datafindcache_list):
    """
    Take as input a list of glue.lal.Cache objects and return a pycbc FileList
    containing all frames within those caches.
   
    Parameters
    -----------
    datafindcache_list : list of glue.lal.Cache objects
        The list of cache files to convert.
  
    Returns
    --------
    datafind_filelist : FileList of frame File objects
        The list of frame files.
    """ 
    datafind_filelist = FileList([])
    prev_file = None
    for cache in datafindcache_list:
        curr_ifo = cache.ifo
        for frame in cache:
            # Don't add a new workflow file entry for this frame if
            # if is a duplicate. These are assumed to be returned in time
            # order
            if prev_file and prev_file.cache_entry.url == frame.url:
                continue

            # Pegasus doesn't like "localhost" in URLs.
            frame.url = frame.url.replace('file://localhost','file://')

            currFile = File(curr_ifo, frame.description,
                    frame.segment, file_url=frame.url, use_tmp_subdirs=True)
            if frame.url.startswith('file://'):
                currFile.PFN(frame.url, site='local')
            else:
                currFile.PFN(frame.url, site='notlocal')
            datafind_filelist.append(currFile)
            prev_file = currFile
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
    --------
    newScienceSegs : Dictionary of ifo keyed glue.segment.segmentlist instances
        The times covered by the frames found in datafindOuts.
    """
    newScienceSegs = {}
    for cache in datafindcaches:
        if len(cache) > 0:
            groupSegs = segments.segmentlist(e.segment for e in cache).coalesce()
            ifo = cache.ifo
            if not newScienceSegs.has_key(ifo):
                newScienceSegs[ifo] = groupSegs
            else:
                newScienceSegs[ifo].extend(groupSegs)
                newScienceSegs[ifo].coalesce()
    return newScienceSegs

def get_missing_segs_from_frame_file_cache(datafindcaches):
    """
    This function will use os.path.isfile to determine if all the frame files
    returned by the local datafind server actually exist on the disk. This can
    then be used to update the science times if needed.
   
    Parameters
    -----------
    datafindcaches : OutGroupList
        List of all the datafind output files.

    Returns
    --------
    missingFrameSegs : Dict. of ifo keyed glue.segment.segmentlist instances
        The times corresponding to missing frames found in datafindOuts.
    missingFrames: Dict. of ifo keyed lal.Cache instances
        The list of missing frames
    """
    missingFrameSegs = {}
    missingFrames = {}
    for cache in datafindcaches:
        if len(cache) > 0:
            # Don't bother if these are not file:// urls, assume all urls in
            # one cache file must be the same type
            if not cache[0].scheme == 'file':
                warn_msg = "We have %s entries in the " %(cache[0].scheme,)
                warn_msg += "cache file. I do not check if these exist."
                logging.info(warn_msg)
                continue
            _, currMissingFrames = cache.checkfilesexist(on_missing="warn")
            missingSegs = segments.segmentlist(e.segment \
                                         for e in currMissingFrames).coalesce()
            ifo = cache.ifo
            if not missingFrameSegs.has_key(ifo):
                missingFrameSegs[ifo] = missingSegs
                missingFrames[ifo] = lal.Cache(currMissingFrames)
            else:
                missingFrameSegs[ifo].extend(missingSegs)
                # NOTE: This .coalesce probably isn't needed as the segments
                # should be disjoint. If speed becomes an issue maybe remove it?
                missingFrameSegs[ifo].coalesce()
                missingFrames[ifo].extend(currMissingFrames)
    return missingFrameSegs, missingFrames

def setup_datafind_server_connection(cp, tag=None):
    """
    This function is resposible for setting up the connection with the datafind
    server.

    Parameters
    -----------
    cp : pycbc.workflow.configuration.WorkflowConfigParser
        The memory representation of the ConfigParser
    Returns
    --------
    connection
        The open connection to the datafind server.
    """
    if cp.has_option_tags("workflow-datafind",
                          "datafind-ligo-datafind-server", [tag]):
        datafind_server = cp.get_opt_tags("workflow-datafind",
                                        "datafind-ligo-datafind-server", [tag])
    else:
        datafind_server = None
        
    return datafind_connection(datafind_server)

def get_segment_summary_times(scienceFile, segmentName):
    """
    This function will find the times for which the segment_summary is set
    for the flag given by segmentName.

    Parameters
    -----------
    scienceFile : SegFile
        The segment file that we want to use to determine this.
    segmentName : string
        The DQ flag to search for times in the segment_summary table.

    Returns
    ---------
    summSegList : glue.segments.segmentlist
        The times that are covered in the segment summary table.
    """
    # Parse the segmentName
    segmentName = segmentName.split(':')
    if not len(segmentName) in [2,3]:
        raise ValueError("Invalid channel name %s." %(segmentName))
    ifo = segmentName[0]
    channel = segmentName[1]
    version = ''
    if len(segmentName) == 3:
        version = int(segmentName[2])

    # Load the filename
    xmldoc = utils.load_filename(scienceFile.cache_entry.path,
                             gz=scienceFile.cache_entry.path.endswith("gz"),
                             contenthandler=ContentHandler)

    # Get the segment_def_id for the segmentName
    segmentDefTable = table.get_table(xmldoc, "segment_definer")
    for entry in segmentDefTable:
        if (entry.ifos == ifo) and (entry.name == channel):
            if len(segmentName) == 2 or (entry.version==version):
                segDefID = entry.segment_def_id
                break
    else:
        raise ValueError("Cannot find channel %s in segment_definer table."\
                         %(segmentName))

    # Get the segmentlist corresponding to this segmentName in segment_summary
    segmentSummTable = table.get_table(xmldoc, "segment_summary")
    summSegList = segments.segmentlist([])
    for entry in segmentSummTable:
        if entry.segment_def_id == segDefID:
            segment = segments.segment(entry.start_time, entry.end_time)
            summSegList.append(segment)
    summSegList.coalesce()
   
    return summSegList

def run_datafind_instance(cp, outputDir, connection, observatory, frameType,
                          startTime, endTime, ifo, tag=None):
    """
    This function will query the datafind server once to find frames between
    the specified times for the specified frame type and observatory.

    Parameters
    ----------
    cp : ConfigParser instance
        This is used to find any kwArgs that should be sent to the datafind
        mod
    outputDir : Output cache files will be written here. We also write the
        commands for reproducing what is done in this function to this
        directory.
    connection : datafind connection object
        Initialized through the glue.datafind module, this is the open
        connection to the datafind server.
    observatory : string
        The observatory to query frames for. This is 'H', 'L' or 'V' and not
        the usual 'H1', 'L1', 'V1' ... because.
    frameType : string
        The frame type to query for.
    startTime : int
        Integer start time to query the datafind server for frames.
    endTime : int
        Integer end time to query the datafind server for frames.
    ifo : string
        The interferometer to use for naming output. This is 'H1', 'L1', 'V1',
        etc. Maybe this could be merged with the observatory string, but this
        could cause issues if running on old 'H2' and 'H1' data.
    tag : string, optional (default=None)
        Use this to specify a tag. This can be used if this module is being
        called more than once to give call specific configuration (by setting
        options in [workflow-datafind-${TAG}] rather than [workflow-datafind]). This
        is also used to tag the Files returned by the class to uniqueify
        the Files and uniqueify the actual filename.
        FIXME: Filenames may not be unique with current codes!


    Returns
    --------
    dfCache : glue.lal.Cache instance
       The glue.lal.Cache representation of the call to the datafind
       server and the returned frame files.
    cacheFile : pycbc.workflow.core.File
        Cache file listing all of the datafind output files for use later in the pipeline.

    """
    seg = segments.segment([startTime, endTime])
    if tag:
        currTags = [tag]
    else:
        currTags = []
    # Take the datafind KWargs from config (usually urltype=file is
    # given).
    dfKwargs = {}
    # By default ignore missing frames, this case is dealt with outside of here
    dfKwargs['on_gaps'] = 'ignore'
    if cp.has_section("datafind"):
        for item, value in cp.items("datafind"):
            dfKwargs[item] = value
    if tag:
        if cp.has_section('datafind-%s' %(tag)):
            for item, value in cp.items("datafind-%s" %(tag)):
                dfKwargs[item] = value

    # It is useful to print the corresponding command to the logs
    # directory to check if this was expected.
    log_datafind_command(observatory, frameType, startTime, endTime,
                         os.path.join(outputDir,'logs'), **dfKwargs)
    logging.debug("Asking datafind server for frames.")
    dfCache = connection.find_frame_urls(observatory, frameType, 
                                        startTime, endTime, **dfKwargs)
    logging.debug("Frames returned")
    # workflow format output file
    cache_file = File(ifo, 'DATAFIND', seg, extension='lcf',
                      directory=outputDir, tags=currTags)
    cache_file.PFN(cache_file.cache_entry.path, site='local')
    
    dfCache.ifo = ifo
    # Dump output to file
    fP = open(cache_file.storage_path, "w")
    # FIXME: CANNOT use dfCache.tofile because it will print 815901601.00000
    #        as a gps time which is incompatible with the lal cache format
    #        (and the C codes) which demand an integer.
    #dfCache.tofile(fP)
    for entry in dfCache:
        start = str(int(entry.segment[0]))
        duration = str(int(abs(entry.segment)))
        print >> fP, "%s %s %s %s %s" \
            %(entry.observatory, entry.description, start, duration, entry.url)
        entry.segment = segments.segment(int(entry.segment[0]), int(entry.segment[1]))

    fP.close()
    return dfCache, cache_file


def log_datafind_command(observatory, frameType, startTime, endTime,
                         outputDir, **dfKwargs):
    """
    This command will print an equivalent gw_data_find command to disk that
    can be used to debug why the internal datafind module is not working.
    """
    # FIXME: This does not accurately reproduce the call as assuming the
    # kwargs will be the same is wrong, so some things need to be converted
    # "properly" to the command line equivalent.
    gw_command = ['gw_data_find', '--observatory', observatory,
                  '--type', frameType,
                  '--gps-start-time', str(startTime),
                  '--gps-end-time', str(endTime)]

    for name, value in dfKwargs.items():
        if name == 'match':
            gw_command.append("--match")
            gw_command.append(str(value))
        elif name == 'urltype':
            gw_command.append("--url-type")
            gw_command.append(str(value))
        elif name == 'on_gaps':
            pass
        else:
            errMsg = "Unknown datafind kwarg given: %s. " %(name)
            errMsg+= "This argument is stripped in the logged .sh command."
            logging.warn(errMsg)
  
    fileName = "%s-%s-%d-%d.sh" \
               %(observatory, frameType, startTime, endTime-startTime)
    filePath = os.path.join(outputDir, fileName)
    fP = open(filePath, 'w')
    fP.write(' '.join(gw_command))
    fP.close()

def datafind_keep_unique_backups(backup_outs, orig_outs):
    """This function will take a list of backup datafind files, presumably
    obtained by querying a remote datafind server, e.g. CIT, and compares
    these against a list of original datafind files, presumably obtained by
    querying the local datafind server. Only the datafind files in the backup
    list that do not appear in the original list are returned. This allows us
    to use only files that are missing from the local cluster.

    Parameters
    -----------
    backup_outs : FileList
        List of datafind files from the remote datafind server.
    orig_outs : FileList
        List of datafind files from the local datafind server.

    Returns
    --------
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
            assert(len(pfns) == 1)
            orig_out.PFN(pfns[0].url, site='notlocal')

    return return_list
