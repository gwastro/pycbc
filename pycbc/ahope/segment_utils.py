import os
import subprocess
from glue import segments
from glue.ligolw import utils,table,lsctables



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

    segFindCall = [ "ligolw_segment_query",
        "--query-segments",
        "--segment-url", config.get("segfind", "segment-url"),
        "--gps-start-time", start,
        "--gps-end-time", end,
        "--include-segments", sciSegName,
        "--output-file", segFindXML ]
   
    errCode = make_external_call(segFindCall, outDir=os.path.join(out_dir, \
                       'logs'), outBaseName='%s-science-call' %(ifo.lower()) )
    if errCode:
        raise subprocess.CalledProcessError(errCode, ' '.join(segFindCall))

    # Yes its yucky to generate a file and then read it back in. This will be
    # fixed when the new API for segment generation is ready.
    sciSegs = fromsegmentxml(sciXmlFile)

    return sciSegs, sciXmlFile

# Function to load segments from an xml file taken from pylal/dq

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
    xmldoc, digest = utils.load_fileobj(file, gz=file.name.endswith(".gz"))
    seg_def_table  = table.get_table(xmldoc, \
                                     lsctables.SegmentDefTable.tableName)
    seg_table      = table.get_table(xmldoc, lsctables.SegmentTable.tableName)

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

