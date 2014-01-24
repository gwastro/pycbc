import os,sys,optparse
import urlparse,urllib
import logging
from glue import datafind
from glue import lal
from glue import segments,segmentsUtils,git_version
from glue.ligolw import utils, table, lsctables, ligolw
from pycbc.ahope import AhopeFile, AhopeFileList

def setup_datafind_workflow(workflow, scienceSegs,  outputDir, 
                            checkSegmentGaps='no_test',
                            checkFramesExist='no_test', 
                            checkSegmentSummary='no_test', segFileDict=None):
    """
    Setup datafind section of ahope workflow. This section is responsible for
    generating, or setting up the workflow to generate, a list of files that
    record the location of the frame files needed to perform the analysis.
    There could be multiple options here, the datafind jobs could be done at
    run time or could be put into a dag. The subsequent jobs will know
    what was done here from the AhopeOutFileList containing the datafind jobs
    (and the Dagman nodes if appropriate.
    For now the only implemented option is to generate the datafind files at
    runtime. This module can also check if the frameFiles actually exist, check
    whether the obtained segments line up with the original ones and update the
    science segments to reflect missing data files.

    Parameters
    ----------
    workflow: Workflow
        The ahope workflow class that stores the jobs that will be run.
    scienceSegs : Dictionary of ifo keyed glue.segment.segmentlist instances
        This contains the times that ahope is expected to analyse.
    outputDir : path
        All output files written by datafind processes will be written to this
        directory.
    checkSegmentGaps : string (optional, default='no_test')
        If this option takes any value other than 'no_test' ahope will check
        that the local datafind server has returned frames covering all of the
        listed science times. Its behaviour is then as follows
        * 'no_test': Do not perform this test. Any discrepancies will cause
          later failures.
        * 'warn': Perform the test, print warnings covering any discrepancies
          but do nothing about them. Discrepancies will cause failures later in
          the workflow.
        * 'update_times': Perform the test, print warnings covering any
          discrepancies and update the input science times to remove times that
          are not present on the host cluster.
        * 'raise_error': Perform the test. If any discrepancies occur, raise a
          ValueError.
    checkFramesExist : string (optional, default='no_test')
        If this options takes any value other than 'no_test' ahope will check
        that the frames returned by the local datafind server are accessible
        from the machine that is running ahope. Its behaviour is then as follows
        * 'no_test': Do not perform this test. Any discrepancies will cause
          later failures.
        * 'warn': Perform the test, print warnings covering any discrepancies
          but do nothing about them. Discrepancies will cause failures later in
          the workflow.
        * 'update_times': Perform the test, print warnings covering any
          discrepancies and update the input science times to remove times that
          are not present on the host cluster.
        * 'raise_error': Perform the test. If any discrepancies occur, raise a
          ValueError.
    checkSegmentSummary: string (optional, default='no_test')
        If this option takes any value other than 'no_test' ahope will check
        that all frames returned by datafind are covered by the segment_summary
        table (for the science flag). Its behaviour is then as follows:
        * 'no_test': Do not perform this test. 
        * 'warn': Perform the test, print warnings covering any discrepancies
          but do nothing about them.
        * 'raise_error': Perform the test. If any discrepancies occur, raise a
          ValueError.

    

    Returns
    --------
    datafindOuts : AhopeOutGroupList
        List of all the datafind output files for use later in the pipeline.
    scienceSegs : Dictionary of ifo keyed glue.segment.segmentlist instances
        This contains the times that ahope is expected to analyse. If the 
        updateSegmentTimes kwarg is given this will be updated to reflect any
        instances of missing data.
    """
    logging.info("Entering datafind module")
    cp = workflow.cp
    
    logging.info("Starting datafind with setup_datafind_runtime_generated")
    if cp.get("ahope-datafind","datafind-method") == "AT_RUNTIME_MULTIPLE":
        datafindcaches, datafindouts = setup_datafind_runtime_generated(cp, scienceSegs,
                                                                    outputDir)
    elif cp.get("ahope-datafind","datafind-method") == "AT_RUNTIME_SINGLE":
        datafindcaches, datafindouts = setup_datafind_runtime_single_call_perifo(cp, 
                                                       scienceSegs, outputDir)
    else:
        msg = "Entry datafind-method in [ahope-datafind] does not have "
        msg += "expected value. Valid values are AT_RUNTIME_MULTIPLE, "
        msg += "AT_RUNTIME_SINGLE.."
        raise ValueError(msg)

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
                msg = "From ifo %s we are missing segments:" %(ifo)
                msg += "\n%s" % "\n".join(map(str, missing))
                missingData = True
                logging.error(msg)
                if checkSegmentGaps == 'update_times':
                    # Remove missing time, so that we can carry on if desired
                    logging.info("Updating science times for ifo %s." %(ifo))
                    scienceSegs[ifo] = scienceSegs[ifo] - missing
        if checkSegmentGaps == 'raise_error' and missingData:
            raise ValueError("Ahope cannot find needed data, exiting.")
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
        for ifo in scienceSegs.keys():
            # If no data in the input then do nothing
            if not scienceSegs[ifo]:
                continue
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
            raise ValueError("Ahope cannot find all frames, exiting.")
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
        if not cp.get("ahope-datafind", "datafind-method")\
                == "AT_RUNTIME_SINGLE":
            logging.info("This is much more meaningful if using "
                         "datafind-method = AT_RUNTIME_SINGLE.")
        dfScienceSegs = get_science_segs_from_datafind_outs(datafindcaches)
        missingFlag = False
        for ifo in dfScienceSegs.keys():
            scienceFile = segFileDict[ifo]['SCIENCE']
            scienceChannel = cp.get('ahope-segments',\
                                'segments-%s-science-name'%(ifo.lower()))
            segSummaryTimes = get_segment_summary_times(scienceFile,
                                                        scienceChannel)
            missing = dfScienceSegs[ifo] - segSummaryTimes
            if abs(missing):
                msg = "From ifo %s the following times have frames, " %(ifo)
                msg += "but are not covered in the segment summary table."
                msg += "\n%s" % "\n".join(map(str, missing))
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


    logging.info("Leaving datafind module")
    return AhopeFileList(datafindouts), scienceSegs
    

def setup_datafind_runtime_generated(cp, scienceSegs, outputDir):
    """
    This function uses the glue.datafind library to obtain the location of all
    the frame files that will be needed to cover the analysis of the data
    given in scienceSegs. This function will not check if the returned frames
    cover the whole time requested, such sanity checks are done in the
    pycbc.ahope.setup_datafind_workflow entry function. As opposed to
    setup_datafind_runtime_single_call_perifo this call will one call to the
    datafind server for every science segment.

    Parameters
    -----------
    cp : ConfigParser.ConfigParser instance
        This contains a representation of the information stored within the
        ahope configuration files
    scienceSegs : Dictionary of ifo keyed glue.segment.segmentlist instances
        This contains the times that ahope is expected to analyse.
    outputDir : path
        All output files written by datafind processes will be written to this
        directory.

    Returns
    --------
    datafindOuts : AhopeOutGroupList
        List of all the datafind output files for use later in the pipeline.

    """
    # First job is to do setup for the datafind jobs
    # First get the server name
    logging.info("Setting up connection to datafind server.")
    connection = setup_datafind_server_connection(cp)

    # Now ready to loop over the input segments
    datafindouts = []
    datafindcaches = []
    ifos = scienceSegs.keys()
    jobTag = "DATAFIND"
    logging.info("Querying datafind server for all science segments.")
    for ifo, scienceSegsIfo in scienceSegs.items():
        observatory = ifo[0].upper()
        frameType = cp.get("ahope-datafind", "datafind-%s-frame-type"%(ifo))
        for seg in scienceSegsIfo:
            msg = "Finding data between %d and %d " %(seg[0],seg[1])
            msg += "for ifo %s" %(ifo)
            logging.debug(msg)
            # WARNING: For now ahope will expect times to be in integer seconds
            startTime = int(seg[0])
            endTime = int(seg[1])

            # Sometimes the connection can drop, so try a backup here
            try:
                cache, cache_file = run_datafind_instance(cp, outputDir, connection,
                                           observatory, frameType, startTime,
                                           endTime, ifo, jobTag)
            except:
                connection = setup_datafind_server_connection(cp)
                cache, cache_file = run_datafind_instance(cp, outputDir, connection,
                                           observatory, frameType, startTime,
                                           endTime, ifo, jobTag)
            datafindouts.append(cache_file)
            datafindcaches.append(cache)
    return datafindcaches, datafindouts

def setup_datafind_runtime_single_call_perifo(cp, scienceSegs, outputDir):
    """
    This function uses the glue.datafind library to obtain the location of all
    the frame files that will be needed to cover the analysis of the data
    given in scienceSegs. This function will not check if the returned frames
    cover the whole time requested, such sanity checks are done in the
    pycbc.ahope.setup_datafind_workflow entry function. As opposed to 
    setup_datafind_runtime_generated this call will only run one call to
    datafind per ifo, spanning the whole time.

    Parameters
    -----------
    cp : ConfigParser.ConfigParser instance
        This contains a representation of the information stored within the
        ahope configuration files
    scienceSegs : Dictionary of ifo keyed glue.segment.segmentlist instances
        This contains the times that ahope is expected to analyse.
    outputDir : path
        All output files written by datafind processes will be written to this
        directory.

    Returns
    --------
    datafindOuts : AhopeOutGroupList
        List of all the datafind output files for use later in the pipeline.

    """
    # First job is to do setup for the datafind jobs
    # First get the server name
    logging.info("Setting up connection to datafind server.")
    connection = setup_datafind_server_connection(cp)

    # We want to ignore gaps as the detectors go up and down and calling this
    # way will give gaps. See the setup_datafind_runtime_generated function
    # for datafind calls that only query for data that will exist
    cp.set("datafind","on_gaps","ignore")

    # Now ready to loop over the input segments
    datafindouts = []
    datafindcaches = []
    ifos = scienceSegs.keys()
    jobTag = "DATAFIND"
    logging.info("Querying datafind server for all science segments.")
    for ifo, scienceSegsIfo in scienceSegs.items():
        observatory = ifo[0].upper()
        frameType = cp.get("ahope-datafind", "datafind-%s-frame-type"%(ifo))
        # This REQUIRES a coalesced segment list to work
        startTime = int(scienceSegsIfo[0][0])
        endTime = int(scienceSegsIfo[-1][1])
        cache, cache_file = run_datafind_instance(cp, outputDir, connection,
                                       observatory, frameType, startTime,
                                       endTime, ifo, jobTag)
        datafindouts.append(cache_file)
        datafindcaches.append(cache)
    return datafindcaches, datafindouts

def get_science_segs_from_datafind_outs(datafindcaches):
    """
    This function will calculate the science segments that are covered in
    the AhopeOutGroupList containing the frame files returned by various
    calls to the datafind server. This can then be used to check whether this
    list covers what it is expected to cover.

    Parameters
    ----------
    datafindOuts : AhopeOutGroupList
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
                # NOTE: This .coalesce probably isn't needed as the segments should
                # be disjoint. If speed becomes an issue maybe remove it?
                newScienceSegs[ifo].coalesce()
    return newScienceSegs

def get_missing_segs_from_frame_file_cache(datafindcaches):
    """
    This function will use os.path.isfile to determine if all the frame files
    returned by the local datafind server actually exist on the disk. This can
    then be used to update the science times if needed.
   
    Parameters
    -----------
    datafindOuts : AhopeOutGroupList
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



def setup_datafind_server_connection(cp):
    """
    This function is resposible for setting up the connection with the datafind
    server.
    """
    if cp.has_option("ahope-datafind", "datafind-ligo-datafind-server"):
        datafindServer = cp.get("ahope-datafind",
                                "datafind-ligo-datafind-server")
    else:
        # Get the server name from the environment
        if os.environ.has_key("LIGO_DATAFIND_SERVER"):
            datafindServer = os.environ["LIGO_DATAFIND_SERVER"]
        else:
            errMsg = "Trying to obtain the ligo datafind server url from "
            errMsg += "the environment, ${LIGO_DATAFIND_SERVER}, but that "
            errMsg += "variable is not populated."
            raise ValueError(errMsg)

    # verify authentication options
    if not datafindServer.endswith("80"):
        cert_file, key_file = datafind.find_credential()
    else:
        cert_file, key_file = None, None

    # Is a port specified in the server URL
    server, port = datafindServer.split(':',1)
    if port == "":
        port = None
    else:
        port = int(port)

    # Open connection to the datafind server
    if cert_file and key_file:
        #HTTPS connection
        connection =\
            datafind.GWDataFindHTTPSConnection(host=server, port=port, 
                                   cert_file=cert_file, key_file=key_file)
    else:
        # HTTP connection
        connection =\
            datafind.GWDataFindHTTPConnection(host=server, port=port)
    return connection

def get_segment_summary_times(scienceFile, segmentName):
    """
    ADD DOCUMENTATION
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
    xmldoc = utils.load_filename(scienceFile.path,
                             gz=scienceFile.path.endswith("gz"),
                             contenthandler=ligolw.DefaultLIGOLWContentHandler)

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
                          startTime, endTime, ifo, jobTag):
    """
    This function will query the datafind server once to find frames between
    the specified times for the specified frame type and observatory.

    Properties
    ----------
    cp : ConfigParser instance
        This is used to find any kwArgs that should be sent to the datafind
        module.
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
    jobTag : string
        Used in the naming of the output files.
    """
    seg = segments.segment([startTime, endTime])
    # Take the datafind KWargs from config (usually urltype=file is
    # given).
    dfKwargs = {}
    # By default ignore missing frames, this case is dealt with outside of here
    dfKwargs['on_gaps'] = 'ignore'
    for item, value in cp.items("datafind"):
        dfKwargs[item] = value
    # It is useful to print the corresponding command to the logs
    # directory to check if this was expected.
    log_datafind_command(observatory, frameType, startTime, endTime,
                         os.path.join(outputDir,'logs'), **dfKwargs)
    logging.debug("Asking datafind server for frames.")
    dfCache = connection.find_frame_urls(observatory, frameType, 
                                        startTime, endTime, **dfKwargs)
    logging.debug("Frames returned")
    # ahope format output file
    cache_file = AhopeFile(ifo, jobTag, seg, 
                           extension='lcf', directory=outputDir)
    dfCache.ifo = ifo
    # Dump output to file
    fP = open(cache_file.path, "w")
    dfCache.tofile(fP)
    fP.close()
    
    return dfCache, cache_file


    
def log_datafind_command(observatory, frameType, startTime, endTime,
                         outputDir, **dfKwargs):
    """
    This command will print an equivalent gw_data_find command to disk that
    can be used to debug why the internal datafind module is not working.
    """
    gw_command = ['gw_data_find', '--observatory', observatory,
                  '--type', frameType,
                  '--gps-start-time', str(startTime),
                  '--gps-end-time', str(endTime)]

    for name, value in dfKwargs.items():
        gw_command.append("--" + name)
        gw_command.append(str(value))
  
    fileName = "%s-%s-%d-%d.sh" \
               %(observatory, frameType, startTime, endTime-startTime)
    filePath = os.path.join(outputDir, fileName)
    fP = open(filePath, 'w')
    fP.write(' '.join(gw_command))
    fP.close()
