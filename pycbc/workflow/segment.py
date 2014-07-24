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
This module is responsible for setting up the segment generation stage of
workflows. For details about this module and its capabilities see here:
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/ahope/segments.html
"""


import os,sys,shutil
import subprocess
import logging
import urllib
from glue import segments, pipeline
from glue.ligolw import utils, table, lsctables, ligolw
from pycbc.workflow.workflow import *
from pycbc.workflow.jobsetup import *

class ContentHandler(ligolw.LIGOLWContentHandler):
        pass

lsctables.use_in(ContentHandler)

def setup_segment_generation(workflow, out_dir, tag=None):
    """
    This function is the gateway for setting up the segment generation steps in a
    workflow. It is designed to be able to support multiple ways of obtaining
    these segments and to combine/edit such files as necessary for analysis.
    The current modules have the capability to generate files at runtime or to
    generate files that are not needed for workflow generation within the workflow.
    
    Parameters
    -----------
    Workflow : workflow.Workflow
        The workflow instance that the coincidence jobs will be added to.
        This instance also contains the ifos for which to attempt to obtain
        segments for this analysis and the start and end times to search for
        segments over.
    out_dir : path
        The directory in which output will be stored.    
    tag : string, optional (default=None)
        Use this to specify a tag. This can be used if this module is being
        called more than once to give call specific configuration (by setting
        options in [workflow-datafind-${TAG}] rather than [workflow-datafind]). This
        is also used to tag the WorkflowFiles returned by the class to uniqueify
        the WorkflowFiles and uniqueify the actual filename.
        FIXME: Filenames may not be unique with current codes!

    Returns
    -------
    segsToAnalyse : dictionay of ifo-keyed glue.segment.segmentlist instances
        This will contain the times that your code should analyse. By default this
        is science time - CAT_1 vetoes. (This default could be changed if desired)
    segFilesList : workflow.WorkflowFileList of workflow.WorkflowSegFile instances
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
    segmentsMethod = cp.get_opt_tag("workflow-segments", "segments-method", tag)
    # These only needed if calling setup_segment_gen_mixed
    if segmentsMethod in ['AT_RUNTIME','CAT2_PLUS_DAG','CAT3_PLUS_DAG',
                          'CAT4_PLUS_DAG']:
        maxVetoCat = cp.get_opt_tag("workflow-segments",
                                    "segments-maximum-veto-category", tag)
        maxVetoCat = int(maxVetoCat)
        veto_categories = range(1,maxVetoCat+1)
        if cp.has_option_tag("workflow-segments",
                         "segments-generate-coincident-segments", tag):
            generate_coincident_segs = True
        else:
            generate_coincident_segs = False
        # Need to curl the veto-definer file
        vetoDefUrl = cp.get_opt_tag("workflow-segments",
                                    "segments-veto-definer-url", tag)
        vetoDefBaseName = os.path.basename(vetoDefUrl)
        vetoDefNewPath = os.path.join(out_dir, vetoDefBaseName)
        urllib.urlretrieve (vetoDefUrl, vetoDefNewPath)
        # and update location
        cp.set("workflow-segments", "segments-veto-definer-file", vetoDefNewPath)

    
    if cp.has_option_tag("workflow-segments", "segments-minimum-segment-length",
                         tag):
        minSegLength = cp.get_opt_tag("workflow-segments", 
                              "segments-minimum-segment-length", tag)
        minSegLength = int(minSegLength)
    else:
        minSegLength = 0

    if segmentsMethod == "AT_RUNTIME":
        max_veto = 1000
    elif segmentsMethod == "CAT2_PLUS_DAG":
        max_veto = 1
    elif segmentsMethod == "CAT3_PLUS_DAG":
        max_veto = 2
    elif segmentsMethod == "CAT4_PLUS_DAG":
        max_veto = 3
    else:
        msg = "Entry segments-method in [workflow-segments] does not have "
        msg += "expected value. Valid values are AT_RUNTIME, CAT4_PLUS_DAG, "
        msg += "CAT2_PLUS_DAG or CAT3_PLUS_DAG."
        raise ValueError(msg)
        
    logging.info("Generating segments with setup_segment_gen_mixed")
    segFilesList = setup_segment_gen_mixed(workflow, veto_categories, 
                             out_dir, max_veto, tag=tag,
                             generate_coincident_segs=generate_coincident_segs)
    logging.info("Segments obtained")

    # This creates the segsToAnalyse from the segFilesList. Currently it uses
    # the 'SCIENCE_OK' segFilesList, which is science - CAT_1 in
    # setup_segment_gen_mixed.
    # This also applies the minimum science length
    segsToAnalyse = {}
    for ifo in workflow.ifos:
        analSegs = segFilesList.find_output_with_ifo(ifo)
        analSegs = analSegs.find_output_with_tag('SCIENCE_OK')
        assert len(analSegs) == 1
        analSegs = analSegs[0]
        if analSegs.segmentList:
            if minSegLength:
                analSegs.removeShortSciSegs(minSegLength)
                analSegs.toSegmentXml()
            segsToAnalyse[ifo] = analSegs.segmentList
        else:
            msg = "No science segments found for ifo %s. " %(ifo)
            msg += "If this is unexpected check the files that were dumped "
            msg += "in the %s directory. Also the " %(out_dir)
            msg += "commands that can be used to reproduce some of these "
            msg += "in %s/*.sh" %(os.path.join(out_dir,'logs'))
            logging.warn(msg)
            
    logging.info("Leaving segment generation module")
    return segsToAnalyse, segFilesList

def setup_segment_gen_mixed(workflow, veto_categories, out_dir, 
                            maxVetoAtRunTime, tag=None,
                            generate_coincident_segs=True):
    """
    This function will generate veto files for each ifo and for each veto
    category.
    It can generate these vetoes at run-time or in the workflow (or do some at
    run-time and some in the workflow). However, the CAT_1 vetoes and science
    time must be generated at run time as they are needed to plan the workflow.
    CATs 2 and higher *may* be needed for other workflow construction.
    It can also combine these files to create a set of cumulative,
    multi-detector veto files, which can be used in ligolw_thinca and in
    pipedown. Again these can be created at run time or within the workflow.

    Parameters
    -----------
    Workflow : workflow.Workflow
        The workflow instance that the coincidence jobs will be added to.
        This instance also contains the ifos for which to attempt to obtain
        segments for this analysis and the start and end times to search for
        segments over.
    veto_categories : list of ints
        List of veto categories to generate segments for. If this stops being
        integers, this can be changed here.
    out_dir : path
        The directory in which output will be stored.    
    maxVetoAtRunTime : int
        Generate veto files at run time up to this category. Veto categories
        beyond this in veto_categories will be generated in the workflow.
        If we move to a model where veto
        categories are not explicitly cumulative, this will be rethought.
    tag : string, optional (default=None)
        Use this to specify a tag. This can be used if this module is being
        called more than once to give call specific configuration (by setting
        options in [workflow-datafind-${TAG}] rather than [workflow-datafind]). This
        is also used to tag the WorkflowFiles returned by the class to uniqueify
        the WorkflowFiles and uniqueify the actual filename.
        FIXME: Filenames may not be unique with current codes!
    generate_coincident_segs : boolean, optional (default = True)
        If given this module will generate a set of coincident, cumulative veto
        files that can be used with ligolw_thinca and pipedown.

    Returns
    -------
    segFilesList : dictionary of workflow.WorkflowSegFile instances
        These are representations of the various segment files that were
        constructed
        at this stage of the workflow and may be needed at later stages of the
        analysis (e.g. for performing DQ vetoes). If the file was generated at
        run-time the segment lists contained within these files will be an
        attribute
        of the instance. (If it will be generated in the workflow it will 
        not be because I am not psychic).
    """
    cp = workflow.cp
    segFilesList = WorkflowFileList([])
    start_time = workflow.analysis_time[0]
    end_time = workflow.analysis_time[1]
    segValidSeg = workflow.analysis_time
    # Will I need to add some jobs to the workflow?
    vetoGenJob = create_segs_from_cats_job(cp, out_dir, workflow.ifo_string)
    
    for ifo in workflow.ifos:
        logging.info("Generating science segments for ifo %s" %(ifo))
        currSciSegs, currSciXmlFile = get_science_segments(ifo, cp, start_time,
                                                    end_time, out_dir, tag=tag)
        segFilesList.append(currSciXmlFile)

        for category in veto_categories:
            if category > maxVetoAtRunTime:
                msg = "Adding creation of CAT_%d segments " %(category)
                msg += "for ifo %s to workflow." %(ifo)
                logging.info(msg)
                execute_status = False
                                 
            if category <= maxVetoAtRunTime:
                logging.info("Generating CAT_%d segments for ifo %s." \
                             %(category,ifo))
                execute_status = True

            currVetoXmlFile = get_veto_segs(workflow, ifo, category, 
                                                start_time, end_time, out_dir,
                                                vetoGenJob, 
                                                execute_now=execute_status)  

            segFilesList.append(currVetoXmlFile) 
            # Store the CAT_1 veto segs for use below
            if category == 1:
                # Yes its yucky to generate a file and then read it back in. 
                #This will be
                # fixed when the new API for segment generation is ready.
                vetoXmlFP = open(currVetoXmlFile.storage_path, 'r')
                cat1Segs = fromsegmentxml(vetoXmlFP)
                vetoXmlFP.close()
                
        analysedSegs = currSciSegs - cat1Segs
        analysedSegs.coalesce()
        analysedXmlFile = os.path.join(out_dir,
                             "%s-SCIENCE_OK_SEGMENTS.xml" %(ifo.upper()) )
        currUrl = urlparse.urlunparse(['file', 'localhost', analysedXmlFile,
                          None, None, None])
        if tag:
            currTags = [tag, 'SCIENCE_OK']
        else:
            currTags = ['SCIENCE_OK']
        currFile = WorkflowOutSegFile(ifo, 'SEGMENTS',
                                   segValidSeg, currUrl, segment_list=analysedSegs,
                                   tags = currTags)
        segFilesList.append(currFile)
        currFile.toSegmentXml()


    if generate_coincident_segs:
        # Need to make some combined category veto files to use when vetoing
        # segments and triggers.
        ifo_string = workflow.ifo_string
        categories = []
        for category in veto_categories:
            categories.append(category)
            # Set file name in workflow standard
            if tag:
                currTags = [tag, 'CUMULATIVE_CAT_%d' %(category)]
            else:
                currTags = ['CUMULATIVE_CAT_%d' %(category)]

            cumulativeVetoFile = os.path.join(out_dir,
                                   '%s-CUMULATIVE_CAT_%d_VETO_SEGMENTS.xml' \
                                   %(ifo_string, category) )
            currUrl = urlparse.urlunparse(['file', 'localhost',
                                         cumulativeVetoFile, None, None, None])
            currSegFile = WorkflowOutSegFile(ifo_string, 'SEGMENTS',
                                   segValidSeg, currUrl, segment_list=analysedSegs,
                                   tags=currTags)
            # And actually make the file (or queue it in the workflow)
            logging.info("Generating combined, cumulative CAT_%d segments."\
                             %(category))
            if category <= maxVetoAtRunTime:
                execute_status = True
            else:
                execute_status = False
            get_cumulative_segs(workflow, currSegFile,  categories,
                                segFilesList, out_dir, 
                                execute_now=execute_status)

            segFilesList.append(currSegFile)

    return segFilesList

#FIXME: Everything below here uses the S6 segment architecture. This is going
# to be replaced in aLIGO with a new architecture. When this is done all of
# the code that follows will need to be replaced with the new version.

def get_science_segments(ifo, cp, start_time, end_time, out_dir, tag=None):
    """
    Obtain science segments for the selected ifo

    Parameters
    -----------
    ifo : string
        The string describing the ifo to obtain science times for.
    start_time : gps time (either int/LIGOTimeGPS)
        The time at which to begin searching for segments.
    end_time : gps time (either int/LIGOTimeGPS)
        The time at which to stop searching for segments.
    out_dir : path
        The directory in which output will be stored.    
    tag : string, optional (default=None)
        Use this to specify a tag. This can be used if this module is being
        called more than once to give call specific configuration (by setting
        options in [workflow-datafind-${TAG}] rather than [workflow-datafind]). This
        is also used to tag the WorkflowFiles returned by the class to uniqueify
        the WorkflowFiles and uniqueify the actual filename.

    Returns
    --------
    sciSegs : glue.segments.segmentlist
        The segmentlist generated by this call
    sciXmlFile : workflow.WorkflowSegFile
        The workflow File object corresponding to this science segments file.

    """
    segValidSeg = segments.segment([start_time,end_time])
    sciSegName = cp.get_opt_tag("workflow-segments","segments-%s-science-name" \
                        %(ifo.lower()), tag ) 
    sciSegUrl = cp.get_opt_tag("workflow-segments","segments-database-url", tag)
    if tag:
        sciXmlFilePath = os.path.join(out_dir, "%s-SCIENCE_SEGMENTS_%s.xml" \
                                       %(ifo.upper(), tag) )
        tagList=[tag, 'SCIENCE']
    else:
        sciXmlFilePath = os.path.join(out_dir, "%s-SCIENCE_SEGMENTS.xml" \
                                           %(ifo.upper()) )
        tagList = ['SCIENCE']

    segFindCall = [ cp.get("executables","segment_query"),
        "--query-segments",
        "--segment-url", sciSegUrl,
        "--gps-start-time", str(start_time),
        "--gps-end-time", str(end_time),
        "--include-segments", sciSegName,
        "--output-file", sciXmlFilePath ]
   
    make_external_call(segFindCall, out_dir=os.path.join(out_dir,'logs'),
                            out_basename='%s-science-call' %(ifo.lower()) )

    # Yes its yucky to generate a file and then read it back in. This will be
    # fixed when the new API for segment generation is ready.
    sciXmlFP = open(sciXmlFilePath,'r')
    sciSegs = fromsegmentxml(sciXmlFP)
    sciXmlFP.close()
    currUrl = urlparse.urlunparse(['file', 'localhost', sciXmlFilePath,
                                   None, None, None])
    sciXmlFile = WorkflowOutSegFile(ifo, 'SEGMENTS',
                                  segValidSeg, currUrl, segment_list=sciSegs,
                                  tags=tagList)

    return sciSegs, sciXmlFile

def get_veto_segs(workflow, ifo, category, start_time, end_time, out_dir, 
                  vetoGenJob, tag=None, execute_now=False):
    """
    Obtain veto segments for the selected ifo and veto category and add the job
    to generate this to the workflow.

    Parameters
    -----------
    workflow: Workflow
        An instance of the workflow.Workflow class that manages the workflow.
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
    workflow : workflow.Workflow
        The workflow instance that the DQ generation Node will be added
        to.
    vetoGenJob : workflow.Job
        The veto generation Job class that will be used to create the Node.
    tag : string, optional (default=None)
        Use this to specify a tag. This can be used if this module is being
        called more than once to give call specific configuration (by setting
        options in [workflow-datafind-${TAG}] rather than [workflow-datafind]). This
        is also used to tag the WorkflowFiles returned by the class to uniqueify
        the WorkflowFiles and uniqueify the actual filename.
        FIXME: Filenames may not be unique with current codes!
    execute_now : boolean, optional
        If true, jobs are executed immediately. If false, they are added to the
        workflow to be run later.

    Returns
    --------
    veto_def_file : workflow.WorkflowSegFile
        The workflow File object corresponding to this DQ veto file.
    """
    segValidSeg = segments.segment([start_time,end_time])
    node = WorkflowNode(vetoGenJob)
    node.add_opt('--veto-categories', str(category))
    node.add_opt('--ifo-list', ifo)
    node.add_opt('--gps-start-time', str(start_time))
    node.add_opt('--gps-end-time', str(end_time))
    vetoXmlFileName = "%s-VETOTIME_CAT%d-%d-%d.xml" \
                         %(ifo, category, start_time, end_time-start_time)
    vetoXmlFilePath = os.path.join(out_dir, vetoXmlFileName)
    currUrl = urlparse.urlunparse(['file', 'localhost',
                                   vetoXmlFilePath, None, None, None])
    if tag:
        currTags = [tag, 'VETO_CAT%d' %(category)]
    else:
        currTags = ['VETO_CAT%d' %(category)]

    vetoXmlFile = WorkflowOutSegFile(ifo, 'SEGMENTS', segValidSeg, currUrl,
                                  tags=currTags)
    node._add_output(vetoXmlFile)
    
    if execute_now:
        workflow.execute_node(node)
    else:
        workflow.add_node(node)
    return vetoXmlFile

def create_segs_from_cats_job(cp, out_dir, ifo_string, tag=None):
    """
    This function creates the CondorDAGJob that will be used to run 
    ligolw_segments_from_cats as part of the workflow

    Parameters
    -----------
    cp : workflow.WorkflowConfigParser
        The in-memory representation of the configuration (.ini) files
    out_dir : path
        Directory in which to put output files
    ifo_string : string
        String containing all active ifos, ie. "H1L1V1"
    tag : string, optional (default=None)
        Use this to specify a tag. This can be used if this module is being
        called more than once to give call specific configuration (by setting
        options in [workflow-datafind-${TAG}] rather than [workflow-datafind]). This
        is also used to tag the WorkflowFiles returned by the class to uniqueify
        the WorkflowFiles and uniqueify the actual filename.
        FIXME: Filenames may not be unique with current codes!

    Returns
    --------
    job : workflow.Job instance
        The Job instance that will run segments_from_cats jobs
    """
    segServerUrl = cp.get_opt_tag("workflow-segments", "segments-database-url",
                                  tag)
    vetoDefFile = cp.get_opt_tag("workflow-segments", "segments-veto-definer-file",
                                 tag)
    if tag:
        currTags = [tag]
    else:
        currTags = []
    job = WorkflowExecutable(cp, 'segments_from_cats', universe='local',
                               ifos=ifo_string, out_dir=out_dir, tags=currTags)
    job.add_opt('--separate-categories')
    job.add_opt('--segment-url', segServerUrl)
    job.add_opt('--veto-file', vetoDefFile)
    # FIXME: Would like the proxy in the Workflow instance
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
    try:
        shutil.copyfile(proxy, proxyfile)
    except IOError:
        raise RuntimeError('Cannot find certificate in %s. '
                           'Make sure that ligo-proxy-init ' 
                           'has been run.' % proxy)
                           
        
    job.add_profile('condor', 'environment',
                       'USER=$ENV(USER);X509_USER_PROXY=%s' % proxyfile)

    return job
    
def get_cumulative_segs(workflow, currSegFile, categories,
                                   segFilesList, out_dir, tags=[],
                                   execute_now=False):
    """
    Function to generate one of the cumulative, multi-detector segment files
    as part of the workflow.
   
    Parameters
    -----------
    workflow: workflow.Workflow
        An instance of the Workflow class that manages the workflow.
    currSegFile : workflow.WorkflowSegFile
        The WorkflowSegFile corresponding to this file that will be created.
    categories : int
        The veto categories to include in this cumulative veto.
    segFilesList : Listionary of workflowSegFiles
        The list of segment files to be used as input for combining.
    out_dir : path
        The directory to write output to.
    tags : list of strings, optional
        A list of strings that is used to identify this job
    execute_now : boolean, optional
        If true, jobs are executed immediately. If false, they are added to the
        workflow to be run later.
    """
    add_inputs = WorkflowFileList([])
    valid_segment = currSegFile.segment
    segment_name = segment_name = 'VETO_CAT%d_CUMULATIVE' % (categories[-1])
    cp = workflow.cp
    # calculate the cumulative veto files for a given ifo
    for ifo in workflow.ifos:
        cum_job = LigoLWCombineSegsExecutable(cp, 'ligolw_combine_segments', 
                       out_dir=out_dir, tags=tags + [segment_name], ifos=ifo)
        inputs = []
        files = segFilesList.find_output_with_ifo(ifo)
        for category in categories:
            fileList = files.find_output_with_tag('VETO_CAT%d' %(category))
            inputs+=fileList                                                      
        
        cum_node = cum_job.create_node(valid_segment, inputs, segment_name)
        if execute_now:
            workflow.execute_node(cum_node)
        else:
            workflow.add_node(cum_node)
        add_inputs += cum_node.output_files
            
    # add cumulative files for each ifo together
    add_job = LigolwAddExecutable(cp, 'llwadd', ifo=ifo, out_dir=out_dir, tags=tags)
    add_node = add_job.create_node(valid_segment, add_inputs,
                                   output=currSegFile)   
    if execute_now:
        workflow.execute_node(add_node)
    else:
        workflow.add_node(add_node)

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
                             contenthandler=ContentHandler)

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
        # Here we want to encode ifo, channel name and version
        full_channel_name = ':'.join([str(seg_def.ifos), str(seg_def.name), 
                                                         str(seg_def.version)])
        seg_id[int(seg_def.segment_def_id)] = full_channel_name
        if dict:
            segs[full_channel_name] = segments.segmentlist()

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
