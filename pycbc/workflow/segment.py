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
import logging
import urllib2, urlparse
import lal
from glue import segments, segmentsUtils
from glue.ligolw import utils, table, lsctables, ligolw
from pycbc.workflow.core import Executable, FileList, Node, OutSegFile, make_analysis_dir, make_external_call, File
from pycbc.workflow.jobsetup import LigolwAddExecutable, LigoLWCombineSegsExecutable

class ContentHandler(ligolw.LIGOLWContentHandler):
    pass

lsctables.use_in(ContentHandler)

# This variable is global to this module and is used to tell the code when to
# regenerate files. NOTE: This only applies to files generated during run-time
# within this module. For files generated within the workflow use the standard
# pegasus reuse facility if you need it. NOTE: It is recommended for most
# applications to use the default option and regenerate all segment files at
# runtime. Only use this facility if you need it
# Options are:
# generate_segment_files='always' : DEFAULT: All files will be generated
#                                     even if they already exist.
# generate_segment_files='if_not_present': Files will be generated if they do
#                                   not already exist. Pre-existing files will
#                                   be read in and used.
# generate_segment_files='error_on_duplicate': Files will be generated if they
#                                   do not already exist. Pre-existing files
#                                   will raise a failure.
# generate_segment_files='never': Pre-existing files will be read in and used.
#                                 If no file exists the code will fail.
global generate_segment_files
generate_segment_files='always'

def setup_segment_generation(workflow, out_dir, tag=None):
    """
    This function is the gateway for setting up the segment generation steps in a
    workflow. It is designed to be able to support multiple ways of obtaining
    these segments and to combine/edit such files as necessary for analysis.
    The current modules have the capability to generate files at runtime or to
    generate files that are not needed for workflow generation within the workflow.
    
    Parameters
    -----------
    workflow : pycbc.workflow.core.Workflow
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
        is also used to tag the Files returned by the class to uniqueify
        the Files and uniqueify the actual filename.
        FIXME: Filenames may not be unique with current codes!

    Returns
    -------
    segsToAnalyse : dictionay of ifo-keyed glue.segment.segmentlist instances
        This will contain the times that your code should analyse. By default this
        is science time - CAT_1 vetoes. (This default could be changed if desired)
    segFilesList : pycbc.workflow.core.FileList of SegFile instances
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
    segmentsMethod = cp.get_opt_tags("workflow-segments",
                                     "segments-method", [tag])
    # These only needed if calling setup_segment_gen_mixed
    if segmentsMethod in ['AT_RUNTIME','CAT2_PLUS_DAG','CAT3_PLUS_DAG',
                          'CAT4_PLUS_DAG']:
        veto_cats = cp.get_opt_tags("workflow-segments",
                                    "segments-veto-categories", [tag])
        max_veto_cat = max([int(c) for c in veto_cats.split(',')])
        veto_categories = range(1, max_veto_cat + 1)
        if cp.has_option_tags("workflow-segments",
                              "segments-generate-coincident-segments", [tag]):
            generate_coincident_segs = True
        else:
            generate_coincident_segs = False
        # Need to curl the veto-definer file
        vetoDefUrl = cp.get_opt_tags("workflow-segments",
                                     "segments-veto-definer-url", [tag])
        vetoDefBaseName = os.path.basename(vetoDefUrl)
        vetoDefNewPath = os.path.join(out_dir, vetoDefBaseName)
        response = urllib2.urlopen(vetoDefUrl)
        html = response.read()
        out_file = open(vetoDefNewPath, 'w')
        out_file.write(html)
        out_file.close()
        # and update location
        cp.set("workflow-segments", "segments-veto-definer-file",
                vetoDefNewPath)

    
    if cp.has_option_tags("workflow-segments",
                          "segments-minimum-segment-length", [tag]):
        minSegLength = int( cp.get_opt_tags("workflow-segments",
                           "segments-minimum-segment-length", [tag]) )
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

    global generate_segment_files
    if cp.has_option_tags("workflow-segments",
                          "segments-generate-segment-files", [tag]):
        value = cp.get_opt_tags("workflow-segments", 
                                   "segments-generate-segment-files", [tag])
        generate_segment_files = value

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
    workflow : pycbc.workflow.core.Workflow
        The Workflow instance that the coincidence jobs will be added to.
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
        is also used to tag the Files returned by the class to uniqueify
        the Files and uniqueify the actual filename.
        FIXME: Filenames may not be unique with current codes!
    generate_coincident_segs : boolean, optional (default = True)
        If given this module will generate a set of coincident, cumulative veto
        files that can be used with ligolw_thinca and pipedown.

    Returns
    -------
    segFilesList : dictionary of pycbc.workflow.core.SegFile instances
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
    segFilesList = FileList([])
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
        currFile = OutSegFile(ifo, 'SEGMENTS',
                              segValidSeg, currUrl, segment_list=analysedSegs,
                              tags = currTags)
        segFilesList.append(currFile)
        currFile.toSegmentXml()


    if generate_coincident_segs:
        # Need to make some combined category veto files to use when vetoing
        # segments and triggers.
        ifo_string = workflow.ifo_string
        categories = []
        cum_cat_files = []
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
            currSegFile = OutSegFile(ifo_string, 'SEGMENTS',
                                   segValidSeg, currUrl, tags=currTags)
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
            cum_cat_files.append(currSegFile)
        # Create a combined file
        # Set file tag in workflow standard
        if tag:
            currTags = [tag, 'COMBINED_CUMULATIVE_SEGMENTS']
        else:
            currTags = ['COMBINED_CUMULATIVE_SEGMENTS']

        combined_veto_file = os.path.join(out_dir,
                               '%s-CUMULATIVE_ALL_CATS_SEGMENTS.xml' \
                               %(ifo_string) )
        curr_url = urlparse.urlunparse(['file', 'localhost',
                                       combined_veto_file, None, None, None])
        curr_file = OutSegFile(ifo_string, 'SEGMENTS',
                               segValidSeg, curr_url, tags=currTags)

        for category in veto_categories:
            if category <= maxVetoAtRunTime:
                execute_status = True
                break
        else:
            execute_status = False
        add_cumulative_files(workflow, curr_file, cum_cat_files, out_dir,
                             execute_now=execute_status)
        segFilesList.append(curr_file)

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
        is also used to tag the Files returned by the class to uniqueify
        the Files and uniqueify the actual filename.

    Returns
    --------
    sciSegs : glue.segments.segmentlist
        The segmentlist generated by this call
    sciXmlFile : pycbc.workflow.core.OutSegFile
        The workflow File object corresponding to this science segments file.

    """
    segValidSeg = segments.segment([start_time,end_time])
    sciSegName = cp.get_opt_tags(
        "workflow-segments", "segments-%s-science-name" %(ifo.lower()), [tag])
    sciSegUrl = cp.get_opt_tags(
        "workflow-segments", "segments-database-url", [tag])
    if tag:
        sciXmlFilePath = os.path.join(
            out_dir, "%s-SCIENCE_SEGMENTS_%s.xml" %(ifo.upper(), tag) )
        tagList=[tag, 'SCIENCE']
    else:
        sciXmlFilePath = os.path.join(
            out_dir, "%s-SCIENCE_SEGMENTS.xml" %(ifo.upper()) )
        tagList = ['SCIENCE']

    if file_needs_generating(sciXmlFilePath):
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
    sciXmlFilePath = os.path.abspath(sciXmlFilePath)
    sciSegs = fromsegmentxml(sciXmlFP)
    sciXmlFP.close()
    currUrl = urlparse.urlunparse(['file', 'localhost', sciXmlFilePath,
                                   None, None, None])
    sciXmlFile = OutSegFile(ifo, 'SEGMENTS',
                                  segValidSeg, currUrl, segment_list=sciSegs,
                                  tags=tagList)
    sciXmlFile.PFN(sciXmlFilePath, site='local')
    return sciSegs, sciXmlFile

def get_veto_segs(workflow, ifo, category, start_time, end_time, out_dir, 
                  vetoGenJob, tag=None, execute_now=False):
    """
    Obtain veto segments for the selected ifo and veto category and add the job
    to generate this to the workflow.

    Parameters
    -----------
    workflow: pycbc.workflow.core.Workflow
        An instance of the Workflow class that manages the workflow.
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
    vetoGenJob : Job
        The veto generation Job class that will be used to create the Node.
    tag : string, optional (default=None)
        Use this to specify a tag. This can be used if this module is being
        called more than once to give call specific configuration (by setting
        options in [workflow-datafind-${TAG}] rather than [workflow-datafind]). This
        is also used to tag the Files returned by the class to uniqueify
        the Files and uniqueify the actual filename.
        FIXME: Filenames may not be unique with current codes!
    execute_now : boolean, optional
        If true, jobs are executed immediately. If false, they are added to the
        workflow to be run later.

    Returns
    --------
    veto_def_file : pycbc.workflow.core.OutSegFile
        The workflow File object corresponding to this DQ veto file.
    """
    segValidSeg = segments.segment([start_time,end_time])
    node = Node(vetoGenJob)
    node.add_opt('--veto-categories', str(category))
    node.add_opt('--ifo-list', ifo)
    node.add_opt('--gps-start-time', str(start_time))
    node.add_opt('--gps-end-time', str(end_time))
    vetoXmlFileName = "%s-VETOTIME_CAT%d-%d-%d.xml" \
                         %(ifo, category, start_time, end_time-start_time)
    vetoXmlFilePath = os.path.abspath(os.path.join(out_dir, vetoXmlFileName))
    currUrl = urlparse.urlunparse(['file', 'localhost',
                                   vetoXmlFilePath, None, None, None])
    if tag:
        currTags = [tag, 'VETO_CAT%d' %(category)]
    else:
        currTags = ['VETO_CAT%d' %(category)]

    vetoXmlFile = OutSegFile(ifo, 'SEGMENTS', segValidSeg, currUrl,
                                  tags=currTags)
    node._add_output(vetoXmlFile)
    
    if execute_now:
        if file_needs_generating(vetoXmlFile.cache_entry.path):
            workflow.execute_node(node, verbatim_exe = True)
        else:
            node.executed = True
            for fil in node._outputs:
                fil.node = None
                fil.PFN(fil.storage_path, site='local')
    else:
        workflow.add_node(node)
    return vetoXmlFile

def create_segs_from_cats_job(cp, out_dir, ifo_string, tag=None):
    """
    This function creates the CondorDAGJob that will be used to run 
    ligolw_segments_from_cats as part of the workflow

    Parameters
    -----------
    cp : pycbc.workflow.configuration.WorkflowConfigParser
        The in-memory representation of the configuration (.ini) files
    out_dir : path
        Directory in which to put output files
    ifo_string : string
        String containing all active ifos, ie. "H1L1V1"
    tag : string, optional (default=None)
        Use this to specify a tag. This can be used if this module is being
        called more than once to give call specific configuration (by setting
        options in [workflow-datafind-${TAG}] rather than [workflow-datafind]). This
        is also used to tag the Files returned by the class to uniqueify
        the Files and uniqueify the actual filename.
        FIXME: Filenames may not be unique with current codes!

    Returns
    --------
    job : Job instance
        The Job instance that will run segments_from_cats jobs
    """
    segServerUrl = cp.get_opt_tags("workflow-segments", 
                                   "segments-database-url", [tag])
    vetoDefFile = cp.get_opt_tags("workflow-segments", 
                                  "segments-veto-definer-file", [tag])

    if tag:
        currTags = [tag]
    else:
        currTags = []
    job = Executable(cp, 'segments_from_cats', universe='local',
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
                                   execute_now=False, segment_name=None):
    """
    Function to generate one of the cumulative, multi-detector segment files
    as part of the workflow.
   
    Parameters
    -----------
    workflow: pycbc.workflow.core.Workflow
        An instance of the Workflow class that manages the workflow.
    currSegFile : pycbc.workflow.core.SegFile
        The SegFile corresponding to this file that will be created.
    categories : int
        The veto categories to include in this cumulative veto.
    segFilesList : Listionary of SegFiles
        The list of segment files to be used as input for combining.
    out_dir : path
        The directory to write output to.
    tags : list of strings, optional
        A list of strings that is used to identify this job
    execute_now : boolean, optional
        If true, jobs are executed immediately. If false, they are added to the
        workflow to be run later.
    """
    add_inputs = FileList([])
    valid_segment = currSegFile.segment
    if segment_name is None:
        segment_name = 'VETO_CAT%d_CUMULATIVE' % (categories[-1])
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
            if file_needs_generating(cum_node.output_files[0].cache_entry.path):
                workflow.execute_node(cum_node, verbatim_exe = True)
            else:
                cum_node.executed = True
                for fil in cum_node._outputs:
                    fil.node = None
                    fil.PFN(fil.storage_path, site='local')
        else:
            workflow.add_node(cum_node)
        add_inputs += cum_node.output_files
            
    # add cumulative files for each ifo together
    add_job = LigolwAddExecutable(cp, 'llwadd', ifo=ifo, out_dir=out_dir, tags=tags)
    add_node = add_job.create_node(valid_segment, add_inputs,
                                   output=currSegFile)   
    if execute_now:
        if file_needs_generating(add_node.output_files[0].cache_entry.path):
            workflow.execute_node(add_node, verbatim_exe = True)
        else:
            add_node.executed = True
            for fil in add_node._outputs:
                fil.node = None
                fil.PFN(fil.storage_path, site='local')
    else:
        workflow.add_node(add_node)
    return add_node.output_files[0]

def add_cumulative_files(workflow, output_file, input_files, out_dir,
                         execute_now=False, tags=[]):
    """
    Function to combine a set of segment files into a single one. This function
    will not merge the segment lists but keep each separate.

    Parameters
    -----------
    workflow: pycbc.workflow.core.Workflow
        An instance of the Workflow class that manages the workflow.
    output_file: pycbc.workflow.core.File
        The output file object
    input_files: pycbc.workflow.core.FileList
        This list of input segment files
    out_dir : path
        The directory to write output to.
    execute_now : boolean, optional
        If true, jobs are executed immediately. If false, they are added to the
        workflow to be run later.
    tags : list of strings, optional
        A list of strings that is used to identify this job
    """
    llwadd_job = LigolwAddExecutable(workflow.cp, 'llwadd', 
                       ifo=output_file.ifo_list, out_dir=out_dir, tags=tags)
    add_node = llwadd_job.create_node(output_file.segment, input_files,
                                   output=output_file)
    if execute_now:
        if file_needs_generating(add_node.output_files[0].cache_entry.path):
            workflow.execute_node(add_node, verbatim_exe = True)
        else:
            add_node.executed = True
            for fil in add_node._outputs:
                fil.node = None
                fil.PFN(fil.storage_path, site='local')
    else:
        workflow.add_node(add_node)
    return add_node.output_files[0]

def fromsegmentxml(xml_file, return_dict=False, select_seg_def_id=None):
    """
    Read a glue.segments.segmentlist from the file object file containing an
    xml segment table.

    Parameters
    -----------
    xml_file : file object
        file object for segment xml file
    return_dict : boolean, optional (default = False)
        returns a glue.segments.segmentlistdict containing coalesced
        glue.segments.segmentlists keyed by seg_def.name for each entry in the
        contained segment_def_table.
    select_seg_def_id : int, optional (default = None)
        returns a glue.segments.segmentlist object containing only those
        segments matching the given segment_def_id integer

    Returns
    --------
    segs : glue.segments.segmentlist instance
        The segment list contained in the file.
    """

    # load xmldocument and SegmentDefTable and SegmentTables
    xmldoc, digest = utils.load_fileobj(xml_file,
                                        gz=xml_file.name.endswith(".gz"),
                                        contenthandler=ContentHandler)
    seg_def_table = table.get_table(xmldoc,
                                    lsctables.SegmentDefTable.tableName)
    seg_table = table.get_table(xmldoc, lsctables.SegmentTable.tableName)

    if return_dict:
        segs = segments.segmentlistdict()
    else:
        segs = segments.segmentlist()

    seg_id = {}
    for seg_def in seg_def_table:
        # Here we want to encode ifo, channel name and version
        full_channel_name = ':'.join([str(seg_def.ifos),
                                      str(seg_def.name),
                                      str(seg_def.version)])
        seg_id[int(seg_def.segment_def_id)] = full_channel_name
        if return_dict:
            segs[full_channel_name] = segments.segmentlist()

    for seg in seg_table:
        seg_obj = segments.segment(
                lal.LIGOTimeGPS(seg.start_time, seg.start_time_ns),
                lal.LIGOTimeGPS(seg.end_time, seg.end_time_ns))
        if return_dict:
            segs[seg_id[int(seg.segment_def_id)]].append(seg_obj)
        elif select_seg_def_id is not None:
            if int(seg.segment_def_id) == select_seg_def_id:
                segs.append(seg_obj)
        else:
            segs.append(seg_obj)

    if return_dict:
        for seg_name in seg_id.values():
            segs[seg_name] = segs[seg_name].coalesce()
    else:
        segs = segs.coalesce()

    xmldoc.unlink()

    return segs

def find_playground_segments(segs):
    '''Finds playground time in a list of segments.

      Playground segments include the first 600s of every 6370s stride starting
      at GPS time 729273613.

      Parameters
      ----------
      segs : segmentfilelist
          A segmentfilelist to find playground segments.
        
      Returns
      -------
      outlist : segmentfilelist
          A segmentfilelist with all playground segments during the input
          segmentfilelist (ie. segs).
    '''

    # initializations
    start_s2 = 729273613
    playground_stride = 6370
    playground_length = 600
    outlist = segments.segmentlist()

    # loop over segments
    for seg in segs:
        start = seg[0]
        end = seg[1]

        # the first playground segment whose end is after the start of seg
        playground_start = start_s2 + playground_stride * ( 1 + \
                     int(start-start_s2-playground_length) / playground_stride)

        while playground_start < end:
            # find start of playground segment
            if playground_start > start:
                ostart = playground_start
            else:
                ostart = start

            playground_end = playground_start + playground_length

            # find end of playground segment
            if playground_end < end:
                oend = playground_end
            else:
                oend = end

            # append segment
            x = segments.segment(ostart, oend)
            outlist.append(x)

            # increment
            playground_start = playground_start + playground_stride

    return outlist

def get_triggered_coherent_segment(workflow, out_dir, sciencesegs, tag=None):
    """
    Construct the coherent network on and off source segments.

    Parameters
    -----------
    workflow : pycbc.workflow.core.Workflow
        The workflow instance that the coincidence jobs will be added to.
        This instance also contains the ifos for which to attempt to obtain
        segments for this analysis and the start and end times to search for
        segments over.
    out_dir : path
        The directory in which output will be stored.
    sciencesegs : dictionary
        Dictionary of science segments produced by
        ahope.setup_segment_generation()
    tag : string, optional (default=None)
        Use this to specify a tag.

    Returns
    --------
    onsource : glue.segments.segmentlistdict
        A dictionary containing the on source segments for network IFOs

    offsource : glue.segments.segmentlistdict
        A dictionary containing the off source segments for network IFOs
    """
    logging.info("Calculating optimal coherent segment.")

    # Load parsed workflow config options
    cp = workflow.cp
    ra = float(os.path.basename(cp.get('workflow', 'ra')))
    dec = float(os.path.basename(cp.get('workflow', 'dec')))
    triggertime = int(os.path.basename(cp.get('workflow', 'trigger-time')))
    
    minbefore = int(os.path.basename(cp.get('workflow-exttrig_segments',
                                            'min-before')))
    minafter = int(os.path.basename(cp.get('workflow-exttrig_segments',
                                           'min-after')))
    minduration = int(os.path.basename(cp.get('workflow-exttrig_segments',
                                              'min-duration')))
    maxduration = int(os.path.basename(cp.get('workflow-exttrig_segments',
                                              'max-duration')))
    onbefore = int(os.path.basename(cp.get('workflow-exttrig_segments',
                                           'on-before')))
    onafter = int(os.path.basename(cp.get('workflow-exttrig_segments',
                                          'on-after')))
    padding = int(os.path.basename(cp.get('workflow-exttrig_segments',
                                          'pad-data')))
    quanta = int(os.path.basename(cp.get('workflow-exttrig_segments',
                                         'quanta')))
    bufferleft = int(cp.get('workflow-exttrig_segments', 'num-buffer-before'))
    bufferright = int(cp.get('workflow-exttrig_segments', 'num-buffer-after'))

    # Check available data segments meet criteria specified in arguments
    sciencesegs = segments.segmentlistdict(sciencesegs)
    sciencesegs = sciencesegs.extract_common(sciencesegs.keys())
    if triggertime not in sciencesegs[sciencesegs.keys()[0]]:
        logging.error("Trigger is not contained within any available segment."
                      " Exiting.")
        sys.exit()

    offsrclist = sciencesegs[sciencesegs.keys()[0]]
    if len(offsrclist) > 1:
        logging.info("Removing network segments that do not contain trigger "
                     "time")
        for seg in offsrclist:
            if triggertime in seg:
                offsrc = seg
    else:
        offsrc = offsrclist[0]

    if (triggertime - minbefore - padding not in offsrc) or (
            triggertime + minafter + padding not in offsrc):
        logging.error("Not enough data either side of trigger time. Exiting.")
        sys.exit()

    if abs(offsrc) < minduration + 2 * padding:
        logging.error("Available network segment shorter than minimum allowed "
                      "duration. Exiting.")
        sys.exit()

    # Will segment duration be the maximum desired length or not?
    if abs(offsrc) >= maxduration + 2 * padding:
        logging.info("Available network science segment duration (%ds) is "
                     "greater than the maximum allowed segment length (%ds). "
                     "Truncating..." % (abs(offsrc), maxduration))
    else:
        logging.info("Available network science segment duration (%ds) is "
                     "less than the maximum allowed segment length (%ds)."
                     % (abs(offsrc), maxduration))

    logging.info("%ds of padding applied at beginning and end of segment."
                 % padding)


    # Construct on-source
    onstart = triggertime - onbefore
    onend = triggertime + onafter
    oncentre = onstart + ((onbefore + onafter) / 2)
    onsrc = segments.segment(onstart, onend)
    logging.info("Constructed ON-SOURCE: duration %ds (%ds before to %ds after"
                 " trigger)."
                 % (abs(onsrc), triggertime - onsrc[0],
                    onsrc[1] - triggertime))
    onsrc = segments.segmentlist([onsrc])

    # Maximal, centred coherent network segment
    idealsegment = segments.segment(int(oncentre - padding -
                                    0.5 * maxduration),
                                    int(oncentre + padding +
                                    0.5 * maxduration))

    # Construct off-source
    if (idealsegment in offsrc):
        offsrc = idealsegment

    elif idealsegment[1] not in offsrc:
        offsrc &= segments.segment(offsrc[1] - maxduration - 2 * padding,
                                   offsrc[1])

    elif idealsegment[0] not in offsrc:
        offsrc &= segments.segment(offsrc[0],
                                   offsrc[0] + maxduration + 2 * padding)

    # Trimming off-source
    excess = (abs(offsrc) - 2 * padding) % quanta
    if excess != 0:
        logging.info("Trimming %ds excess time to make OFF-SOURCE duration a "
                     "multiple of %ds" % (excess, quanta))
        offset = (offsrc[0] + abs(offsrc) / 2.) - oncentre
        if 2 * abs(offset) > excess:
            if offset < 0:
                offsrc &= segments.segment(offsrc[0] + excess,
                                           offsrc[1])
            elif offset > 0:
                offsrc &= segments.segment(offsrc[0],
                                           offsrc[1] - excess)
            assert abs(offsrc) % quanta == 2 * padding
        else:
            logging.info("This will make OFF-SOURCE symmetrical about trigger "
                         "time.")
            start = int(offsrc[0] - offset + excess / 2)
            end = int(offsrc[1] - offset - round(float(excess) / 2))
            offsrc = segments.segment(start, end)
            assert abs(offsrc) % quanta == 2 * padding

    logging.info("Constructed OFF-SOURCE: duration %ds (%ds before to %ds "
                 "after trigger)."
                 % (abs(offsrc) - 2 * padding,
                    triggertime - offsrc[0] - padding,
                    offsrc[1] - triggertime - padding))
    offsrc = segments.segmentlist([offsrc])

    # Put segments into segmentlistdicts
    onsource = segments.segmentlistdict()
    offsource = segments.segmentlistdict()
    ifos = ''
    for iifo in sciencesegs.keys():
        ifos += str(iifo)
        onsource[iifo] = onsrc
        offsource[iifo] = offsrc

    # Write off-source to xml file
    XmlFile = os.path.join(out_dir,
                           "%s-COH_OFFSOURCE_SEGMENT.xml" % ifos.upper())
    currUrl = urlparse.urlunparse(['file', 'localhost', XmlFile, None, None,
                                   None])
    currFile = OutSegFile(ifos, 'COH-OFFSOURCE', offsource[iifo], currUrl,
                          offsource[iifo])
    currFile.toSegmentXml()
    logging.info("Optimal coherent segment calculated.")

    #FIXME: For legacy coh_PTF_post_processing we must write out the segments
    #       to segwizard files (with hardcoded names!)
    offsourceSegfile = os.path.join(out_dir, "offSourceSeg.txt")
    segmentsUtils.tosegwizard(open(offsourceSegfile, "w"), offsrc)

    onsourceSegfile = os.path.join(out_dir, "onSourceSeg.txt")
    segmentsUtils.tosegwizard(file(onsourceSegfile, "w"), onsrc)

    onlen = abs(onsrc[0])
    bufferSegment = segments.segment(onstart - bufferleft * onlen,
                                     onend + bufferright * onlen)
    bufferSegfile = os.path.join(out_dir, "bufferSeg.txt")
    segmentsUtils.tosegwizard(file(bufferSegfile, "w"),
                              segments.segmentlist([bufferSegment]))

    return onsource, offsource

def save_veto_definer(cp, out_dir, tags=[]):
    """ Retrieve the veto definer file and save it locally
    
    Parameters
    -----------
    cp : ConfigParser instance
    out_dir : path
    tags : list of strings
        Used to retrieve subsections of the ini file for
        configuration options.
    """
    make_analysis_dir(out_dir)
    vetoDefUrl = cp.get_opt_tags("workflow-segments",
                                 "segments-veto-definer-url", tags)
    vetoDefBaseName = os.path.basename(vetoDefUrl)
    vetoDefNewPath = os.path.abspath(os.path.join(out_dir, vetoDefBaseName))
    response = urllib2.urlopen(vetoDefUrl)
    html = response.read()
    out_file = open(vetoDefNewPath, 'w')
    out_file.write(html)
    out_file.close()

    # and update location
    cp.set("workflow-segments", "segments-veto-definer-file", vetoDefNewPath)
    
def parse_cat_ini_opt(cat_str):
    """ Parse a cat str from the ini file into a list of sets """
    if cat_str == "":
        return []    
        
    cat_groups = cat_str.split(',')
    cat_sets = []
    for group in cat_groups:
        group = group.strip()
        cat_sets += [set(c for c in group)]
    return cat_sets  
    
def cat_to_pipedown_cat(val):
    """ Convert a category character to the pipedown equivelant
    
    Parameters
    ----------
    str : single character string
        The input category character
        
    Returns
    -------
    pipedown_str : str
        The pipedown equivelant notation that can be passed to programs 
    that expect this definition.
    """
    if val == '1':
        return 1
    if val == '2':
        return 2
    if val == '3':
        return 4
    if val == 'H':
        return 3 
    else:
        raise ValueError('Invalid Category Choice')

    return cat_sets

def get_analyzable_segments(workflow, out_dir, tags=[]):
    """
    Get the analyzable segments after applying ini specified vetoes.

    Parameters
    -----------
    workflow : Workflow object
        Instance of the workflow object
    out_dir : path
        Location to store output files
    tags : list of strings
        Used to retrieve subsections of the ini file for
        configuration options.

    Returns
    --------
    segs : glue.segments.segmentlist instance
        The segment list specifying the times to analyze
    data_segs : glue.segments.segmentlist
        The segment list specifying the time where data exists
    seg_files : workflow.core.FileList instance
        The cumulative segment files from each ifo that determined the
        analyzable time.
    """
    from pycbc.events import segments_to_file
    logging.info('Entering generation of science segments')
    segments_method = workflow.cp.get_opt_tags("workflow-segments", 
                                      "segments-method", tags)
    
    make_analysis_dir(out_dir)
    start_time = workflow.analysis_time[0]
    end_time = workflow.analysis_time[1]
    save_veto_definer(workflow.cp, out_dir, tags)
    
    cat_sets = parse_cat_ini_opt(workflow.cp.get_opt_tags('workflow-segments',
                                                'segments-science-veto', tags))
    if len(cat_sets) > 1: 
        raise ValueError('Provide only 1 category group to determine'
                         ' analyzable segments')
    cat_set = cat_sets[0]
    
    veto_gen_job = create_segs_from_cats_job(workflow.cp, out_dir, 
                                             workflow.ifo_string) 
    sci_segs, data_segs = {}, {}
    seg_files = FileList()
    for ifo in workflow.ifos:
        sci_segs[ifo], sci_xml = get_science_segments(ifo, workflow.cp, 
                                                 start_time, end_time, out_dir) 
        seg_files += [sci_xml]    
        data_segs[ifo] = sci_segs[ifo]  
        for category in cat_set:
            curr_veto_file = get_veto_segs(workflow, ifo, 
                                        cat_to_pipedown_cat(category), 
                                        start_time, end_time, out_dir,
                                        veto_gen_job, execute_now=True)  
            f = open(curr_veto_file.storage_path, 'r')
            cat_segs = fromsegmentxml(f)
            f.close()    
            sci_segs[ifo] -= cat_segs
            
        sci_segs[ifo].coalesce()
        seg_ok_path = os.path.abspath(os.path.join(out_dir,
                                                   '%s-SCIENCE-OK.xml' % ifo))
        curr_url = urlparse.urlunparse(['file', 'localhost', seg_ok_path,
                          None, None, None])
        curr_file = OutSegFile(ifo, 'SCIENCE_OK', workflow.analysis_time,
                               curr_url, segment_list=sci_segs[ifo],
                               tags=tags)
        curr_file.PFN(seg_ok_path, 'local')
        curr_file.toSegmentXml()
        seg_files += [curr_file]

    if segments_method == 'ALL_SINGLE_IFO_TIME':
        pass
    elif segments_method == 'COINC_TIME':
        cum_segs = None
        for ifo in sci_segs:
            if cum_segs is not None:
                cum_segs = (cum_segs & sci_segs[ifo]).coalesce() 
            else:
                cum_segs = sci_segs[ifo]
                
        for ifo in sci_segs:
            sci_segs[ifo] = cum_segs 
    else:
        raise ValueError("Invalid segments-method, %s. Options are "
                         "ALL_SINGLE_IFO_TIME and COINC_TIME" % segments_method)

    logging.info('Leaving generation of science segments')
    return sci_segs, data_segs, seg_files
    
def get_cumulative_veto_group_files(workflow, option, out_dir, tags=[]):
    """
    Get the cumulative veto files that define the different backgrounds 
    we want to analyze, defined by groups of vetos.

    Parameters
    -----------
    workflow : Workflow object
        Instance of the workflow object
    option : str
        ini file option to use to get the veto groups
    out_dir : path
        Location to store output files
    tags : list of strings
        Used to retrieve subsections of the ini file for
        configuration options.

    Returns
    --------
    seg_files : workflow.core.FileList instance
        The cumulative segment files for each veto group.   
    cat_files : workflow.core.FileList instance
        The list of individual category veto files
    """
    make_analysis_dir(out_dir)
    start_time = workflow.analysis_time[0]
    end_time = workflow.analysis_time[1]

    cat_sets = parse_cat_ini_opt(workflow.cp.get_opt_tags('workflow-segments',
                                            option, tags))
    veto_gen_job = create_segs_from_cats_job(workflow.cp, out_dir,
                                             workflow.ifo_string) 
    cats = set()
    for cset in cat_sets:
        cats = cats.union(cset)
    
    cat_files = FileList()
    for ifo in workflow.ifos:
        for category in cats:
            cat_files.append(get_veto_segs(workflow, ifo,
                                        cat_to_pipedown_cat(category), 
                                        start_time, end_time, out_dir,
                                        veto_gen_job, execute_now=True))

    cum_seg_files = FileList()     
    names = []   
    for cat_set in cat_sets:
        segment_name = "CUMULATIVE_CAT_%s" % (''.join(sorted(cat_set)))
        logging.info('getting information for %s' % segment_name)
        categories = [cat_to_pipedown_cat(c) for c in cat_set]
        path = os.path.join(out_dir, '%s-%s_VETO_SEGMENTS.xml' \
                            % (workflow.ifo_string, segment_name))
        path = os.path.abspath(path)
        url = urlparse.urlunparse(['file', 'localhost', path, None, None, None])
        seg_file = File(workflow.ifos, 'CUM_VETOSEGS', workflow.analysis_time,
                        file_url=url, tags=[segment_name])
                        
        cum_seg_files += [get_cumulative_segs(workflow, seg_file,  categories,
              cat_files, out_dir, execute_now=True, segment_name=segment_name)]
        names.append(segment_name)
              
    return cum_seg_files, names, cat_files

def file_needs_generating(file_path):
    """
    This job tests the file location and determines if the file should be
    generated now or if an error should be raised. This uses the 
    generate_segment_files variable, global to this module, which is described
    above and in the documentation.
   
    Parameters
    -----------
    file_path : path
        Location of file to check

    Returns
    --------
    int
        1 = Generate the file. 0 = File already exists, use it. Other cases
        will raise an error.
    """
    global generate_segment_files
    # Does the file exist
    if os.path.isfile(file_path):
        if generate_segment_files in ['if_not_present', 'never']:
            return 0
        elif generate_segment_files == 'always':
            err_msg = "File %s already exists. " %(file_path,)
            err_msg += "Regenerating and overwriting."
            logging.warn(err_msg)
            return 1
        elif generate_segment_files == 'error_on_duplicate':
            err_msg = "File %s already exists. " %(file_path,)
            err_msg += "Refusing to overwrite file and exiting."
            raise ValueError(err_msg)
        else:
            err_msg = 'Global variable generate_segment_files must be one of '
            err_msg += '"always", "if_not_present", "error_on_duplicate", '
            err_msg += '"never". Got %s.' %(generate_segment_files,)
            raise ValueError(err_msg)
    else:
        if generate_segment_files in ['always', 'if_not_present', 
                                      'error_on_duplicate']:
            return 1
        elif generate_segment_files == 'never':
            err_msg = 'File %s does not exist. ' %(file_path,)
            raise ValueError(err_msg)
        else:
            err_msg = 'Global variable generate_segment_files must be one of '
            err_msg += '"always", "if_not_present", "error_on_duplicate", '
            err_msg += '"never". Got %s.' %(generate_segment_files,)
            raise ValueError(err_msg)
