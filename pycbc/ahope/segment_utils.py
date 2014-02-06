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
This module is responsible for setting up the segment generation stage of ahope
workflows. For details about this module and its capabilities see here:
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/ahope/segments.html
"""


import os,sys,shutil
import subprocess
import logging
import pycurl
from glue import segments, pipeline
from glue.ligolw import utils, table, lsctables, ligolw
from pycbc.ahope.ahope_utils import *

def setup_segment_generation(workflow, ifos, start_time, end_time, out_dir):
    """
    This function is the gateway for setting up the segment generation steps in an
    ahope workflow. It is designed to be able to support multiple ways of obtaining
    these segments and to combine/edit such files as necessary for analysis.
    The current modules have the capability to generate files at runtime or to
    generate files that are not needed for workflow generation within the workflow.
    
    Parameters
    -----------
    Workflow : ahope.Workflow
        The ahope workflow instance that the coincidence jobs will be added to.
    ifos : list of ifo strings
        The ifos for which to attempt to obtain segments for this analysis.
    start_time : gps time (either int/LIGOTimeGPS)
        The time at which to begin searching for segments.
    end_time : gps time (either int/LIGOTimeGPS)
        The time at which to stop searching for segments.
    out_dir : path
        The directory in which output will be stored.    
    maxVetoCat : int, optional (default = 5)
        Generate veto files up to this category. If we move to a model where veto
        categories are not explicitly cumulative, this will be removed.
    minSegLength : int, optional (default = 0)
        This specifies the minimum length of science data to consider. If data
        segments are shorter than this length they will be discarded at this stage.
    generate_coincident_segs : boolean, optional (default = True)
        If given this module will generate a set of coincident, cumulative veto
        files that can be used with ligolw_thinca and pipedown.

    Returns
    -------
    segsToAnalyse : dictionay of ifo-keyed glue.segment.segmentlist instances
        This will contain the times that your code should analyse. By default this
        is science time - CAT_1 vetoes. (This default could be changed if desired)
    segFilesDict : dictionary of ahope.AhopeSegFile instances
        These are representations of the various segment files that were constructed
        at this stage of the workflow and may be needed at later stages of the
        analysis (e.g. for performing DQ vetoes). If the file was generated at
        run-time the segment lists contained within these files will be an attribute
        of the instance. (If it will be generated in the workflow it will not be
        because I am not psychic).
    """
    logging.info("Entering segment generation module")
    make_analysis_dir(out_dir)
    
    cp = workflow.cp

    # Parse for options in ini file
    segmentsMethod = cp.get("ahope-segments","segments-method")
    # These only needed if calling setup_segment_gen_mixed
    if segmentsMethod in ['AT_RUNTIME','CAT2_PLUS_DAG','CAT3_PLUS_DAG',
                          'CAT4_PLUS_DAG']:
        maxVetoCat = cp.get("ahope-segments", "segments-maximum-veto-category")
        maxVetoCat = int(maxVetoCat)
        veto_categories = range(1,maxVetoCat+1)
        if cp.has_option("ahope-segments",
                         "segments-generate-coincident-segments"):
            generate_coincident_segs = True
        else:
            generate_coincident_segs = False
        # Need to curl the veto-definer file
        vetoDefUrl = cp.get("ahope-segments", "segments-veto-definer-url")
        vetoDefBaseName = os.path.basename(vetoDefUrl)
        vetoDefNewPath = os.path.join(out_dir, vetoDefBaseName)
        fp = open(vetoDefNewPath, "wb")
        curl = pycurl.Curl()
        curl.setopt(pycurl.URL, vetoDefUrl)
        curl.setopt(pycurl.WRITEDATA, fp)
        curl.perform()
        curl.close()
        fp.close()
        # and update location
        cp.set("ahope-segments", "segments-veto-definer-file", vetoDefNewPath)

    
    if cp.has_option("ahope-segments", "segments-minimum-segment-length"):
        minSegLength = cp.get("ahope-segments", 
                              "segments-minimum-segment-length")
        minSegLength = int(minSegLength)
    else:
        minSegLength = 0

    if segmentsMethod == "AT_RUNTIME":
        logging.info("Generating segments with setup_segment_gen_mixed")
        segFilesDict = setup_segment_gen_mixed(workflow, ifos, veto_categories, 
                             start_time, end_time, out_dir, 1000,
                             generate_coincident_segs=generate_coincident_segs)
    elif segmentsMethod == "CAT2_PLUS_DAG":
        logging.info("Generating segments with setup_segment_gen_mixed")
        segFilesDict = setup_segment_gen_mixed(workflow, ifos, veto_categories, 
                             start_time, end_time, out_dir, 1,
                             generate_coincident_segs=generate_coincident_segs)
    elif segmentsMethod == "CAT3_PLUS_DAG":
        logging.info("Generating segments with setup_segment_gen_mixed")
        segFilesDict = setup_segment_gen_mixed(workflow, ifos, veto_categories, 
                             start_time, end_time, out_dir, 2,
                             generate_coincident_segs=generate_coincident_segs)
    elif segmentsMethod == "CAT4_PLUS_DAG":
        logging.info("Generating segments with setup_segment_gen_mixed")
        segFilesDict = setup_segment_gen_mixed(workflow, ifos, veto_categories, 
                             start_time, end_time, out_dir, 3,
                             generate_coincident_segs=generate_coincident_segs)
    else:
        msg = "Entry segments-method in [ahope-segments] does not have "
        msg += "expected value. Valid values are AT_RUNTIME, CAT4_PLUS_DAG, "
        msg += "CAT2_PLUS_DAG or CAT3_PLUS_DAG."
        raise ValueError(msg)
    logging.info("Segments obtained")

    # This creates the segsToAnalyse from the segFilesDict. Currently it uses
    # the 'ANALYSED' segFilesDict, which is science - CAT_1 in
    # setup_segment_gen_mixed.
    # This also applies the minimum science length
    segsToAnalyse = {}
    for ifo in ifos:
        if segFilesDict[ifo]['ANALYSED'].segmentList:
            if minSegLength:
                segFilesDict[ifo]['ANALYSED'].removeShortSciSegs(minSegLength)
                segFilesDict[ifo]['ANALYSED'].toSegmentXml()
            segsToAnalyse[ifo] = segFilesDict[ifo]['ANALYSED'].segmentList
        else:
            msg = "No science segments found for ifo %s. " %(ifo)
            msg += "If this is unexpected check the files that were dumped "
            msg += "in the %s directory. Also the " %(out_dir)
            msg += "commands that can be used to reproduce some of these "
            msg += "in %s/*.sh" %(os.path.join(out_dir,'logs'))
            logging.warn(msg)
            
    logging.info("Leaving segment generation module")
    return segsToAnalyse, segFilesDict

def setup_segment_gen_mixed(workflow, ifos, veto_categories, start_time,
                            end_time, out_dir, maxVetoAtRunTime,
                            generate_coincident_segs=True):
    """
    This function will generate veto files for each ifo and for each veto category.
    It can generate these vetoes at run-time or in the workflow (or do some at
    run-time and some in the workflow). However, the CAT_1 vetoes and science time
    must be generated at run time as they are needed to plan the workflow. CATs 2
    and higher *may* be needed for other workflow construction.
    It can also combine these files to create a set of cumulative, multi-detector
    veto files, which can be used in ligolw_thinca and in pipedown. Again these can
    be created at run time or within the workflow.

    Parameters
    -----------
    Workflow : hope.Workflow
        The ahope workflow instance that the coincidence jobs will be added to.
    ifos : list of ifo strings
        The ifos for which to attempt to obtain segments for this analysis.
    veto_categories : list of ints
        List of veto categories to generate segments for. If this stops being
        integers, this can be changed here.
    start_time : gps time (either int/LIGOTimeGPS)
        The time at which to begin searching for segments.
    end_time : gps time (either int/LIGOTimeGPS)
        The time at which to stop searching for segments.
    out_dir : path
        The directory in which output will be stored.    
    maxVetoAtRunTime : int
        Generate veto files at run time up to this category. Veto categories beyond
        this in veto_categories will be generated in the workflow.
        If we move to a model where veto
        categories are not explicitly cumulative, this will be rethought.
    generate_coincident_segs : boolean, optional (default = True)
        If given this module will generate a set of coincident, cumulative veto
        files that can be used with ligolw_thinca and pipedown.

    Returns
    -------
    segFilesDict : dictionary of ahope.AhopeSegFile instances
        These are representations of the various segment files that were constructed
        at this stage of the workflow and may be needed at later stages of the
        analysis (e.g. for performing DQ vetoes). If the file was generated at
        run-time the segment lists contained within these files will be an attribute
        of the instance. (If it will be generated in the workflow it will not be
        because I am not psychic).

    """
    # FIXME: With the tags structures added to the AhopeFile class, I see no reason
    # why segFilesDict cannot be a standard AhopeFileList
    cp = workflow.cp
    segFilesDict = {}
    segValidSeg = segments.segment([start_time,end_time])
    # Will I need to add some jobs to the workflow?
    if max(veto_categories) > maxVetoAtRunTime:
        vetoGenJob = create_segs_from_cats_job(cp, out_dir)
    for ifo in ifos:
        logging.info("Generating science segments for ifo %s" %(ifo))
        segFilesDict[ifo] = {}
        currSciSegs, currSciXmlFile = get_science_segments(ifo, cp, 
                                          start_time, end_time, out_dir)
        segFilesDict[ifo]['SCIENCE'] = currSciXmlFile

        # FIXME: Do I really need these?
        vetoSegs = {}
        vetoXmlFiles = {}
        for category in veto_categories:
            currTags=['VETO_CAT%d' %(category)]
            if category <= maxVetoAtRunTime:
                logging.info("Generating CAT_%d segments for ifo %s." \
                             %(category,ifo))
                vetoSegs[category],vetoXmlFiles[category] = \
                    get_veto_segs_at_runtime(ifo, category, cp, start_time, 
                                             end_time, out_dir, currTags)
            else:
                msg = "Adding creation of CAT_%d segments " %(category)
                msg += "for ifo %s to workflow." %(ifo)
                logging.info(msg)
                vetoXmlFiles[category] = get_veto_segs_in_workflow(ifo, 
                                  category, start_time, end_time, out_dir,
                                  workflow, vetoGenJob)        
                # Don't know what the segments are as they haven't been
                # calculated yet!
                vetoSegs[category] = None
            segFilesDict[ifo][currTags[0]] = vetoXmlFiles[category] 
                
        analysedSegs = currSciSegs - vetoSegs[1]
        analysedSegs.coalesce()
        analysedXmlFile = sciXmlFile = os.path.join(out_dir,
                             "%s-ANALYSED_SEGMENTS.xml" %(ifo.upper()) )
        currUrl = urlparse.urlunparse(['file', 'localhost', analysedXmlFile,
                          None, None, None])
        segFilesDict[ifo]['ANALYSED'] = AhopeOutSegFile(ifo, 'SEGMENTS',
                                 segValidSeg, currUrl, segList=analysedSegs,
                                 tags = ['ANALYSED'])
        segFilesDict[ifo]['ANALYSED'].toSegmentXml()


    if generate_coincident_segs:
        # Need to make some combined category veto files to use when vetoing
        # segments and triggers.
        ifoString = ''.join(ifos)
        segFilesDict[ifoString] = {}
        for category in veto_categories:
            # Set file name in ahope standard
            cumulativeVetoFile = os.path.join(out_dir,
                                   '%s-CUMULATIVE_CAT_%d_VETO_SEGMENTS.xml' \
                                   %(ifoString, category) )
            currUrl = urlparse.urlunparse(['file', 'localhost',
                                         cumulativeVetoFile, None, None, None])
            currSegFile = AhopeOutSegFile(ifoString, 'SEGMENTS',
                                   segValidSeg, currUrl, segList=analysedSegs,
                                   tags=['CUMULATIVE_CAT_%d' %(category)])
            # And actually make the file (or queue it in the workflow)
            if category <= maxVetoAtRunTime:
                logging.info("Generating combined, cumulative CAT_%d segments."\
                             %(category))
                get_cumulative_segs_at_runtime(ifos, currSegFile, category, cp,
                                               segFilesDict, out_dir)
            else:
                errMsg = "Generating segments in the workflow is temporarily "
                errMsg += "disabled as ligolw_segments_compat cannot be added "
                errMsg += "to the ahope workflow without breaking pegasus."
                raise NotImplementedError(errMsg)
            segFilesDict[ifoString]['CUMULATIVE_CAT_%d' %(category)]\
                                                                 = currSegFile

    return segFilesDict

#FIXME: Everything below here uses the S6 segment architecture. This is going
# to be replaced in aLIGO with a new architecture. When this is done all of
# the code that follows will need to be replaced with the new version.

def get_science_segments(ifo, cp, start_time, end_time, out_dir):
    """
    Obtain science segments for the selected ifo

    Properties
    -----------
    ifo : string
        The string describing the ifo to obtain science times for.
    start_time : gps time (either int/LIGOTimeGPS)
        The time at which to begin searching for segments.
    end_time : gps time (either int/LIGOTimeGPS)
        The time at which to stop searching for segments.
    out_dir : path
        The directory in which output will be stored.    

    Returns
    --------
    sciSegs : glue.segments.segmentlist
        The segmentlist generated by this call
    sciXmlFile : ahope.AhopeSegFile
        The ahope File object corresponding to this science segments file.

    """
    segValidSeg = segments.segment([start_time,end_time])
    sciSegName = cp.get("ahope-segments","segments-%s-science-name" \
                        %(ifo.lower()) ) 
    sciSegUrl = cp.get("ahope-segments","segments-database-url")
    sciXmlFilePath = os.path.join(out_dir, "%s-SCIENCE_SEGMENTS.xml" \
                                       %(ifo.upper()) )

    segFindCall = [ cp.get("executables","segment_query"),
        "--query-segments",
        "--segment-url", sciSegUrl,
        "--gps-start-time", str(start_time),
        "--gps-end-time", str(end_time),
        "--include-segments", sciSegName,
        "--output-file", sciXmlFilePath ]
   
    make_external_call(segFindCall, outDir=os.path.join(out_dir,'logs'),
                            outBaseName='%s-science-call' %(ifo.lower()) )

    # Yes its yucky to generate a file and then read it back in. This will be
    # fixed when the new API for segment generation is ready.
    sciXmlFP = open(sciXmlFilePath,'r')
    sciSegs = fromsegmentxml(sciXmlFP)
    sciXmlFP.close()
    currUrl = urlparse.urlunparse(['file', 'localhost', sciXmlFilePath,
                                   None, None, None])
    sciXmlFile = AhopeOutSegFile(ifo, 'SEGMENTS',
                                  segValidSeg, currUrl, segList=sciSegs,
                                  tags = ['SCIENCE'])

    return sciSegs, sciXmlFile

def get_veto_segs_at_runtime(ifo, category, cp, start_time, end_time, out_dir, tags):
    """
    Obtain veto segments for the selected ifo and veto category

    Properties
    -----------
    ifo : string
        The string describing the ifo to generate vetoes for.
    category : int
        The veto category to generate vetoes for.
    start_time : gps time (either int/LIGOTimeGPS)
        The time at which to begin searching for segments.
    end_time : gps time (either int/LIGOTimeGPS)
        The time at which to stop searching for segments.
    out_dir : path
        The directory in which output will be stored.    
    tags : list of strings
        The tags that describe this job

    Returns
    --------
    vetoSegs : glue.segments.segmentlist
        The segmentlist generated by this call
    vetoXmlFile : ahope.AhopeSegFile
        The ahope File object corresponding to this DQ veto file.
    """
    segValidSeg = segments.segment([start_time,end_time])
    segServerUrl = cp.get("ahope-segments", "segments-database-url")
    vetoDefFile = cp.get("ahope-segments", "segments-veto-definer-file")

    segFromCatsCall = [ cp.get("executables","segments_from_cats"),
#       FIXME: I want to use separate categories here, but to do so needs some
#       extra code added to create the cumulative files in the format thinca
#       and pipedown expect
#        "--separate-categories", 
        "--cumulative-categories",
        "--segment-url", segServerUrl,
        "--veto-file", vetoDefFile,
        "--output-dir", out_dir,
        # FIXME: And set this to just str(category) when not using cumulative
        "--veto-categories", ','.join([str(i+1) for i in range(category)]),
        "--ifo-list", ifo,
        "--gps-start-time", str(start_time),
        "--gps-end-time", str(end_time)]

    make_external_call(segFromCatsCall, outDir=os.path.join(out_dir,'logs'),
              outBaseName='%s-veto-cats-%d-call' %(ifo.lower(),category) )

    vetoXmlFileName = "%s-VETOTIME_CAT%d-%d-%d.xml" \
                         %(ifo, category, start_time, end_time-start_time)
    vetoXmlFilePath = os.path.join(out_dir, vetoXmlFileName)
    currUrl = urlparse.urlunparse(['file', 'localhost',
                                   vetoXmlFilePath, None, None, None])
    vetoXmlFile = AhopeOutSegFile(ifo, 'SEGMENTS', segValidSeg, currUrl, tags=tags)

    # Yes its yucky to generate a file and then read it back in. This will be
    # fixed when the new API for segment generation is ready.
    vetoXmlFP = open(vetoXmlFilePath, 'r')
    vetoSegs = fromsegmentxml(vetoXmlFP)
    vetoXmlFP.close()
    vetoXmlFile = AhopeOutSegFile(ifo, 'SEGMENTS', segValidSeg, currUrl, tags=tags,
                                  segList=vetoSegs)

    return vetoSegs, vetoXmlFile

def get_veto_segs_in_workflow(ifo, category, start_time, end_time, out_dir,
                              workflow, vetoGenJob, tags):
    """
    Obtain veto segments for the selected ifo and veto category and add the job
    to generate this to the workflow.

    Properties
    -----------
    ifo : string
        The string describing the ifo to generate vetoes for.
    category : int
        The veto category to generate vetoes for.
    start_time : gps time (either int/LIGOTimeGPS)
        The time at which to begin searching for segments.
    end_time : gps time (either int/LIGOTimeGPS)
        The time at which to stop searching for segments.
    out_dir : path
        The directory in which output will be stored.    
    workflow : ahope.Workflow
        The ahope workflow instance that the DQ generation Node will be added to.
    vetoGenJob : ahope.Job
        The veto generation Job class that will be used to create the Node.
    tags : list of strings
        The tags that describe this job

    Returns
    --------
    veto_def_file : ahope.AhopeSegFile
        The ahope File object corresponding to this DQ veto file.
    """
    segValidSeg = segments.segment([start_time,end_time])
    node = Node(vetoGenJob)
    node.add_var_opt('veto-categories', str(category))
    node.add_var_opt('ifo-list', ifo)
    node.add_var_opt('gps-start-time', str(start_time))
    node.add_var_opt('gps-end-time', str(end_time))
    vetoXmlFileName = "%s-VETOTIME_CAT%d-%d-%d.xml" \
                         %(ifo, category, start_time, end_time-start_time)
    vetoXmlFilePath = os.path.join(out_dir, vetoXmlFileName)
    currUrl = urlparse.urlunparse(['file', 'localhost',
                                   vetoXmlFilePath, None, None, None])
    vetoXmlFile = AhopeOutSegFile(ifo, 'SEGMENTS', segValidSeg, currUrl, tags=tags)
    node.add_output(vetoXmlFile)
    workflow.add_node(node)
    return vetoXmlFile

def create_segs_from_cats_job(cp, out_dir):
    """
    This function creates the CondorDAGJob that will be used to run 
    ligolw_segments_from_cats as part of the workflow

    Parameters
    -----------
    cp : ahope.AhopeConfigParser
        The in-memory representation of the configuration (.ini) files
    out_dir : path
        Directory in which to put output files

    Returns
    --------
    job : ahope.Job instance
        The Job instance that will run segments_from_cats jobs
    """
    # FIXME: Why is this not an ahope.Job class?!?
    segServerUrl = cp.get("ahope-segments", "segments-database-url")
    vetoDefFile = cp.get("ahope-segments", "segments-veto-definer-file")
    exeName = cp.get("executables","segments_from_cats")
    tag = "vetogen"
    logtag = '$(cluster)-$(process)'
    job = pipeline.CondorDAGJob("vanilla", exeName)
    job.set_sub_file('%s.sub' %(tag))
    job.set_stderr_file(os.path.join(out_dir,'logs','%s-%s.err' %(tag,logtag)))
    job.set_stdout_file(os.path.join(out_dir,'logs','%s-%s.out' %(tag,logtag)))
    job.add_condor_cmd('getenv','True')
    job.add_opt('separate-categories', '')
    job.add_opt('output-dir', out_dir)
    job.add_opt('segment-url', segServerUrl)
    job.add_opt('veto-file', vetoDefFile)
    job.exe_name = exeName 
    # FIXME: Explore using the x509 condor commands
    # Set up proxy to be accessible in a NFS location
    # If the user has logged in with gsissh then X509_USER_PROXY will be set
    # However, certain users log in with an ssh key and then ligo-proxy-init
    # This route does not set X509_USER_PROXY, so use the default file location
    if os.environ.has_key('X509_USER_PROXY'):
        proxy = os.getenv('X509_USER_PROXY')
    else:
        proxy = "/tmp/x509up_u%d" % os.getuid()
    proxyfile = os.path.join(out_dir, 'x509up.file')
    shutil.copyfile(proxy, proxyfile)
    job.add_condor_cmd('environment',
                       'USER=$ENV(USER);X509_USER_PROXY=%s' % proxyfile)
    # For now this falls back to local universe as some clusters do not
    # allow nodes to have access to the WWW.
    job.set_universe('local')

    return job

def get_cumulative_segs_at_runtime(ifos, currSegFile, category, cp,
                                   segFilesDict, out_dir):
    """
    Function to generate one of the cumulative, multi-detector segment files at
    runtime.
   
    Parameters
    -----------
    ifos : list
        List of ifos contained in the output file
    currSegFile : ahope.AhopeSegFile
        The AhopeSegFile corresponding to this file that will be created.
    category : int
        The veto category to cumulatively include up to in this file.
    cp : ahope.AhopeConfigParser
        The in-memory representation of the configuration (.ini) files
    segFilesDict : Dictionary of ahopeSegFiles
        The list of segment files to be used as input for combining.
    out_dir : path
        The directory to write output to.
    """
    ifoString = ''.join(ifos)
    # First need to determine the input files
    inputs = get_cumulative_segs_input_files(ifos, segFilesDict, category)

    # Construct the call to ligolw_add
    ligolwAddCall = [ cp.get("executables","llwadd"),
        "--output",
        currSegFile.path]
    ligolwAddCall.extend([inp.path for inp in inputs])

    make_external_call(ligolwAddCall, outDir=os.path.join(out_dir,'logs'),
              outBaseName='%s-gen-cum-cats-%d-call' %(ifoString, category) )

    # FIXME: I want this removed. It cannot go into the workflow as it modifies
    # files in place. It also shouldn't be needed!
    compatVetoCall = [ cp.get("executables","ligolw_segments_compat"),
        currSegFile.path]

    make_external_call(compatVetoCall, outDir=os.path.join(out_dir,'logs'),
              outBaseName='%s-compat-veto-cats-%d-call' %(ifoString, category) )

def get_cumulative_segs_input_files(ifos, segFilesDict, category):
    """
    This function is responsible for identifying which files should be used as input
    when generating a combined, multi-detector veto file.

    Parameters
    -----------
    ifos : list
        List of ifos contained in the output file
    segFilesDict : Dictionary of ahopeSegFiles
        The list of segment files to be used as input for combining.
    category : int
        The veto category to cumulatively include up to in this file.

    Returns
    --------
    fileList : ahope.AhopeFileList
        List of files to use an inputs.
    """
    ifoString = ''.join(ifos)
    fileList = AhopeFileList([])
    for ifo in ifos:
        currTag='VETO_CAT%d' %(category)
        fileList.append(segFilesDict[ifo][currTag])
    # FIXME: Add this back in when not using cumulative categories
    #if category > 1:
    #    currTag = 'CUMULATIVE_CAT_%d' %(category - 1)
    #    fileList.append(segFilesDict[ifoString][currTag])
    return fileList
    

def fromsegmentxml(file, dict=False, id=None):
    """
    Read a glue.segments.segmentlist from the file object file containing an
    xml segment table.

    Parameters
    -----------
    file : file object
        file object for segment xml file
    dict : boolean, optional (default = False)
        returns a glue.segments.segmentlistdict containing coalesced
        glue.segments.segmentlists keyed by seg_def.name for each entry in the
        contained segment_def_table. Default False
    id : int, optional (default = None)
        returns a glue.segments.segmentlist object containing only those
        segments matching the given segment_def_id integer

    Returns
    --------
    segs : glue.segments.segmentlist instance
        The segment list contained in the file.
    """

    # load xmldocument and SegmentDefTable and SegmentTables
    xmldoc, digest = utils.load_fileobj(file, gz=file.name.endswith(".gz"),
                             contenthandler=ligolw.DefaultLIGOLWContentHandler)

    seg_def_table  = table.get_table(xmldoc,
                                     lsctables.SegmentDefTable.tableName)
    seg_table      = table.get_table(xmldoc, lsctables.SegmentTable.tableName)

    for seg in seg_table:
        pass

    if dict:
        segs = segments.segmentlistdict()
    else:
        segs = segments.segmentlist()

    seg_id = {}
    for seg_def in seg_def_table:
        seg_id[int(seg_def.segment_def_id)] = str(seg_def.name)
        if dict:
            segs[str(seg_def.name)] = segments.segmentlist()

    for seg in seg_table:
        if dict:
            segs[seg_id[int(seg.segment_def_id)]]\
                .append(segments.segment(seg.start_time, seg.end_time))
            continue
        if id and int(seg.segment_def_id)==id:
            segs.append(segments.segment(seg.start_time, seg.end_time))
            continue
        segs.append(segments.segment(seg.start_time, seg.end_time))

    if dict:
        for seg_name in seg_id.values():
            segs[seg_name] = segs[seg_name].coalesce()
    else:
        segs = segs.coalesce()

    xmldoc.unlink()

    return segs
