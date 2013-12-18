import os,sys,shutil
import subprocess
import logging
from glue import segments, pipeline
from glue.ligolw import utils, table, lsctables, ligolw
from pycbc.ahope.ahope_utils import *

def setup_segment_generation(workflow, ifos, start_time, end_time, out_dir, 
                             maxVetoCat = 5, minSegLength=0):
    """
    Setup the segment generation needed in an ahope workflow
    FIXME: Add more DOCUMENTATION
    """
    logging.info("Entering segment generation module")
    veto_categories = range(1,maxVetoCat)
    
    cp = workflow.cp

    if cp.get("ahope-segments","segments-method") == "AT_RUNTIME":
        logging.info("Generating segments with setup_segment_gen_runtime")
        segFilesDict = setup_segment_gen_runtime(cp, ifos, veto_categories, 
                                 start_time, end_time, out_dir)
    elif cp.get("ahope-segments","segments-method") == "CAT2_PLUS_DAG":
        logging.info("Generating segments with setup_segment_gen_mixed")
        segFilesDict = setup_segment_gen_mixed(cp, ifos, veto_categories, 
                                 start_time, end_time, out_dir, workflow, 1)
    elif cp.get("ahope-segments","segments-method") == "CAT3_PLUS_DAG":
        logging.info("Generating segments with setup_segment_gen_mixed")
        segFilesDict = setup_segment_gen_mixed(cp, ifos, veto_categories, 
                                 start_time, end_time, out_dir, workflow, 2)
    elif cp.get("ahope-segments","segments-method") == "CAT4_PLUS_DAG":
        logging.info("Generating segments with setup_segment_gen_mixed")
        segFilesDict = setup_segment_gen_mixed(cp, ifos, veto_categories, 
                                 start_time, end_time, out_dir, workflow, 3)
    else:
        msg = "Entry segments-method in [ahope-segments] does not have "
        msg += "expected value. Valid values are AT_RUNTIME, CAT4_PLUS_DAG."
        raise ValueError(msg)
    logging.info("Segments obtained")

    segsToAnalyse = {}
    for ifo in ifos:
        if segFilesDict[ifo]['ANALYSED'].segmentList:
            if minSegLength:
                segFilesDict[ifo]['ANALYSED'].removeShortSciSegs(minSegLength)
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


def setup_segment_gen_runtime(cp, ifos, veto_categories, start_time,
                              end_time, out_dir):
    """
    ADD DOCUMENTATION
    """
    segFilesDict = {}
    segValidSeg = segments.segment([start_time,end_time])
    baseTag='SEGS'
    for ifo in ifos:
        logging.info("Generating science segments for ifo %s" %(ifo))
        ifoTag=baseTag + "_%s" %(ifo.upper())
        segFilesDict[ifo] = {}
        currSciSegs, currSciXmlFile = get_science_segments(ifo, cp, 
                                          start_time, end_time, out_dir)
        currTag = 'SCIENCE'
        currUrl = urlparse.urlunparse(['file', 'localhost', currSciXmlFile,
                          None, None, None])
        segFilesDict[ifo][currTag] = AhopeOutSegFile(ifo, 
                                 '%s_%s' %(ifoTag, currTag), 
                                 segValidSeg, currUrl, segList=currSciSegs)

        vetoSegs = {}
        vetoXmlFiles = {} 
        for category in veto_categories:
            logging.info("Generating CAT_%d segments for ifo %s." \
                         %(category,ifo))
            vetoSegs[category],vetoXmlFiles[category] = \
                get_veto_segs_at_runtime(ifo, category, cp, start_time, 
                                         end_time, out_dir)
            currTag='VETO_CAT%d' %(category)
            currUrl = urlparse.urlunparse(['file', 'localhost',
                          vetoXmlFiles[category], None, None, None])
            segFilesDict[ifo][currTag] = AhopeOutSegFile(ifo, 
                                 '%s_%s' %(ifoTag, currTag), 
                                 segValidSeg, currUrl, \
                                 segList=vetoSegs[category])
        analysedSegs = currSciSegs - vetoSegs[1]
        analysedSegs.coalesce()
        analysedXmlFile = sciXmlFile = os.path.join(out_dir,
                             "%s-ANALYSED_SEGMENTS.xml" %(ifo.upper()) ) 
        currUrl = urlparse.urlunparse(['file', 'localhost', analysedXmlFile,
                          None, None, None])
        currTag='ANALYSED'
        segFilesDict[ifo][currTag] = AhopeOutSegFile(ifo, 
                                 '%s_%s' %(ifoTag, currTag), 
                                 segValidSeg, currUrl, segList=analysedSegs)
        segFilesDict[ifo][currTag].toSegmentXml()
    return segFilesDict

def setup_segment_gen_mixed(cp, ifos, veto_categories, start_time,
                            end_time, out_dir, workflow, maxVetoAtRunTime):
    """
    ADD DOCUMENTATION
    """
    segFilesDict = {}
    segValidSeg = segments.segment([start_time,end_time])
    baseTag='SEGS'
    vetoGenJob = create_segs_from_cats_job(cp, out_dir)
    for ifo in ifos:
        logging.info("Generating science segments for ifo %s" %(ifo))
        ifoTag=baseTag + "_%s" %(ifo.upper())
        segFilesDict[ifo] = {}
        currSciSegs, currSciXmlFile = get_science_segments(ifo, cp, 
                                          start_time, end_time, out_dir)
        currTag = 'SCIENCE'
        currUrl = urlparse.urlunparse(['file', 'localhost', currSciXmlFile,
                          None, None, None])
        segFilesDict[ifo][currTag] = AhopeOutSegFile(ifo, 
                                 '%s_%s' %(ifoTag, currTag), 
                                 segValidSeg, currUrl, segList=currSciSegs)

        vetoSegs = {}
        vetoXmlFiles = {}
        for category in veto_categories:
            currTag='VETO_CAT%d' %(category)
            if category <= maxVetoAtRunTime:
                logging.info("Generating CAT_%d segments for ifo %s." \
                             %(category,ifo))
                vetoSegs[category],vetoXmlFiles[category] = \
                    get_veto_segs_at_runtime(ifo, category, cp, start_time, 
                                             end_time, out_dir)
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
            currUrl = urlparse.urlunparse(['file', 'localhost',
                              vetoXmlFiles[category], None, None, None])
            segFilesDict[ifo][currTag] = AhopeOutSegFile(ifo, 
                                 '%s_%s' %(ifoTag, currTag), 
                                 segValidSeg, currUrl,
                                 segList=vetoSegs[category])
                
        analysedSegs = currSciSegs - vetoSegs[1]
        analysedSegs.coalesce()
        analysedXmlFile = sciXmlFile = os.path.join(out_dir,
                             "%s-ANALYSED_SEGMENTS.xml" %(ifo.upper()) )
        currUrl = urlparse.urlunparse(['file', 'localhost', analysedXmlFile,
                          None, None, None])
        currTag='ANALYSED'
        segFilesDict[ifo][currTag] = AhopeOutSegFile(ifo, 
                                 '%s_%s' %(ifoTag, currTag), 
                                 segValidSeg, currUrl, segList=analysedSegs)
        segFilesDict[ifo][currTag].toSegmentXml()
    return segFilesDict

#FIXME: Everything below here uses the S6 segment architecture. This is going
# to be replaced in aLIGO with a new architecture. When this is done all of
# the code that follows will need to be replaced with the new version.

def get_science_segments(ifo, cp, start_time, end_time, out_dir):
    """
    Obtain science segments for the selected ifo
    """
    sciSegName = cp.get("ahope-segments","segments-%s-science-name" \
                        %(ifo.lower()) ) 
    sciSegUrl = cp.get("ahope-segments","segments-database-url")
    sciXmlFile = os.path.join(out_dir, "%s-SCIENCE_SEGMENTS.xml" \
                                       %(ifo.upper()) )

    segFindCall = [ cp.get("executables","segment_query"),
        "--query-segments",
        "--segment-url", sciSegUrl,
        "--gps-start-time", str(start_time),
        "--gps-end-time", str(end_time),
        "--include-segments", sciSegName,
        "--output-file", sciXmlFile ]
   
    make_external_call(segFindCall, outDir=os.path.join(out_dir,'logs'),
                            outBaseName='%s-science-call' %(ifo.lower()) )

    # Yes its yucky to generate a file and then read it back in. This will be
    # fixed when the new API for segment generation is ready.
    sciXmlFP = open(sciXmlFile,'r')
    sciSegs = fromsegmentxml(sciXmlFP)
    sciXmlFP.close()

    return sciSegs, sciXmlFile

def get_veto_segs_at_runtime(ifo, category, cp, start_time, end_time, out_dir):
    """
    Obtain veto segments for the selected ifo and veto category
    """
    segServerUrl = cp.get("ahope-segments", "segments-database-url")
    vetoDefFile = cp.get("ahope-segments", "segments-veto-definer-file")

    segFromCatsCall = [ cp.get("executables","segments_from_cats"),
        "--separate-categories", 
        "--segment-url", segServerUrl,
        "--veto-file", vetoDefFile,
        "--output-dir", out_dir,
        "--veto-categories", str(category),
        "--ifo-list", ifo,
        "--gps-start-time", str(start_time),
        "--gps-end-time", str(end_time)]

    make_external_call(segFromCatsCall, outDir=os.path.join(out_dir,'logs'),
              outBaseName='%s-veto-cats-%d-call' %(ifo.lower(),category) )

    vetoDefXmlFileName = "%s-VETOTIME_CAT%d-%d-%d.xml" \
                         %(ifo, category, start_time, end_time-start_time)

    vetoDefXmlFile = os.path.join(out_dir, vetoDefXmlFileName)

    # Yes its yucky to generate a file and then read it back in. This will be
    # fixed when the new API for segment generation is ready.
    vetoDefXmlFP = open(vetoDefXmlFile,'r')
    vetoSegs = fromsegmentxml(vetoDefXmlFP)
    vetoDefXmlFP.close()

    return vetoSegs, vetoDefXmlFile

def get_veto_segs_in_workflow(ifo, category, start_time, end_time, out_dir,
                              workflow, vetoGenJob):
    """
    Obtain veto segments for the selected ifo and veto category and add the job
    to generate this to the workflow.
    """
    node = Node(vetoGenJob)
    node.add_var_opt('veto-categories', str(category))
    node.add_var_opt('ifo-list', ifo)
    node.add_var_opt('gps-start-time', str(start_time))
    node.add_var_opt('gps-end-time', str(end_time))
    veto_def_file = AhopeFile(ifo, 'VETOTIME_CAT%d' % category, 
                              segments.segment(start_time, end_time),
                              extension='.xml', directory=out_dir)
    node.add_output(veto_def_file)
    workflow.add_node(node)
    return veto_def_file.path

def create_segs_from_cats_job(cp, out_dir):
    """
    This function creates the CondorDAGJob that will be used to run 
    ligolw_segments_from_cats as part of the workflow
    """
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
    # set up proxy to be accessible in a NFS location
    proxy = os.getenv('X509_USER_PROXY')
#    if os.path.exists(proxy):
    # FIXME: For now this falls back to local universe as some clusters do not
    # allow nodes to have access to the WWW.
    if 0:
        proxyfile = os.path.join(out_dir, 'x509up.file')
        shutil.copyfile(proxy, proxyfile)
        job.add_condor_cmd('environment',\
                              'USER=$ENV(USER);X509_USER_PROXY=%s' % proxyfile)
    else:
        # Fall back is to run in local universe
        job.set_universe('local')

    return job


# Function to load segments from an xml file taken from pylal/dq
# FIXME: Use the pylal/pylal/dq function.
# Need to understand the load_fileobj warning first
def fromsegmentxml(file, dict=False, id=None):

    """
    Read a glue.segments.segmentlist from the file object file containing an
    xml segment table.

    Arguments:

      file : file object
        file object for segment xml file

    Keyword Arguments:

      dict : [ True | False ]
        returns a glue.segments.segmentlistdict containing coalesced
        glue.segments.segmentlists keyed by seg_def.name for each entry in the
        contained segment_def_table. Default False
      id : int
        returns a glue.segments.segmentlist object containing only those
        segments matching the given segment_def_id integer
        
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
