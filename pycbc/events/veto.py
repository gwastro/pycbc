""" This module contains utilities to manipulate trigger lists based on 
segment.
"""
import numpy
from pycbc_glue.ligolw import table, lsctables, utils as ligolw_utils
from pycbc_glue.segments import segment, segmentlist


def start_end_to_segments(start, end):
    return segmentlist([segment(s, e) for s, e in zip(start, end)])

def segments_to_start_end(segs):
    segs.coalesce()
    return (numpy.array([s[0] for s in segs]), 
            numpy.array([s[1] for s in segs]))

def start_end_from_segments(segment_file):
    """
    Return the start and end time arrays from a segment file.
    
    Parameters
    ----------
    segment_file: xml segment file
    
    Returns
    -------
    start: numpy.ndarray
    end: numpy.ndarray
    """
    from pycbc_glue.ligolw.ligolw import LIGOLWContentHandler as h; lsctables.use_in(h)
    indoc = ligolw_utils.load_filename(segment_file, False, contenthandler=h)
    segment_table  = table.get_table(indoc, lsctables.SegmentTable.tableName)
    start = numpy.array(segment_table.getColumnByName('start_time'))
    start_ns = numpy.array(segment_table.getColumnByName('start_time_ns'))
    end = numpy.array(segment_table.getColumnByName('end_time'))
    end_ns = numpy.array(segment_table.getColumnByName('end_time_ns'))
    return start + start_ns * 1e-9, end + end_ns * 1e-9


def indices_within_times(times, start, end):
    """
    Return an index array into times that lie within the durations defined by start end arrays
    
    Parameters
    ----------
    times: numpy.ndarray
        Array of times
    start: numpy.ndarray
        Array of duration start times
    end: numpy.ndarray 
        Array of duration end times
    
    Returns
    -------
    indices: numpy.ndarray
        Array of indices into times
    """
    # coalesce the start/end segments
    start, end = segments_to_start_end(start_end_to_segments(start, end).coalesce())

    tsort = times.argsort()
    times_sorted = times[tsort]
    left = numpy.searchsorted(times_sorted, start)
    right = numpy.searchsorted(times_sorted, end)

    if len(left) == 0:
        return numpy.array([], dtype=numpy.uint32)

    return tsort[numpy.hstack(numpy.r_[s:e] for s, e in zip(left, right))]

def indices_outside_times(times, start, end):
    """
    Return an index array into times that like outside the durations defined by start end arrays
    
    Parameters
    ----------
    times: numpy.ndarray
        Array of times
    start: numpy.ndarray
        Array of duration start times
    end: numpy.ndarray 
        Array of duration end times
    
    Returns
    -------
    indices: numpy.ndarray
        Array of indices into times
    """
    exclude = indices_within_times(times, start, end)
    indices = numpy.arange(0, len(times))
    return numpy.delete(indices, exclude)

def select_segments_by_definer(segment_file, segment_name=None, ifo=None):
    """ Return the list of segments that match the segment name
    
    Parameters
    ----------
    segment_file: str
        path to segment xml file
    
    segment_name: str
        Name of segment
    ifo: str, optional
    
    Returns
    -------
    seg: list of segments
    """
    from pycbc_glue.ligolw.ligolw import LIGOLWContentHandler as h; lsctables.use_in(h)
    indoc = ligolw_utils.load_filename(segment_file, False, contenthandler=h)
    segment_table  = table.get_table(indoc, 'segment')

    seg_def_table = table.get_table(indoc, 'segment_definer')
    def_ifos = seg_def_table.getColumnByName('ifos')
    def_names = seg_def_table.getColumnByName('name')
    def_ids = seg_def_table.getColumnByName('segment_def_id')
    
    valid_id = []
    for def_ifo, def_name, def_id in zip(def_ifos, def_names, def_ids):
        if ifo and ifo != def_ifo:
            continue        
        if segment_name and segment_name != def_name:
            continue
        valid_id += [def_id]

    start = numpy.array(segment_table.getColumnByName('start_time'))
    start_ns = numpy.array(segment_table.getColumnByName('start_time_ns'))
    end = numpy.array(segment_table.getColumnByName('end_time'))
    end_ns = numpy.array(segment_table.getColumnByName('end_time_ns'))
    start, end = start + 1e-9 * start_ns, end + 1e-9 * end_ns
    did = segment_table.getColumnByName('segment_def_id')
    
    keep = numpy.array([d in valid_id for d in did])
    if sum(keep) > 0:
        return start_end_to_segments(start[keep], end[keep])
    else:
        return segmentlist([])

def indices_within_segments(times, segment_files, ifo=None, segment_name=None):
    """ Return the list of indices that should be vetoed by the segments in the
    list of veto_files.
    
    Parameters
    ----------
    times: numpy.ndarray of integer type
        Array of gps start times
    segment_files: string or list of strings
        A string or list of strings that contain the path to xml files that
        contain a segment table
    ifo: string, optional
        The ifo to retrieve segments for from the segment files
    segment_name: str, optional
        name of segment       
    Returns
    -------
    indices: numpy.ndarray
        The array of index values within the segments
    segmentlist: 
        The segment list corresponding to the selected time.
    """
    veto_segs = segmentlist([])
    indices = numpy.array([], dtype=numpy.uint32)   
    for veto_file in segment_files:
        veto_segs += select_segments_by_definer(veto_file, segment_name, ifo)
    veto_segs.coalesce()  
    
    start, end = segments_to_start_end(veto_segs)
    if len(start) > 0:
        idx = indices_within_times(times, start, end)
        indices = numpy.union1d(indices, idx)

    return indices, veto_segs.coalesce()
 
def indices_outside_segments(times, segment_files, ifo=None, segment_name=None):
    """ Return the list of indices that are outside the segments in the
    list of segment files.
    
    Parameters
    ----------
    times: numpy.ndarray of integer type
        Array of gps start times
    segment_files: string or list of strings
        A string or list of strings that contain the path to xml files that
        contain a segment table
    ifo: string, optional
        The ifo to retrieve segments for from the segment files
    segment_name: str, optional
        name of segment               
    Returns
    --------
    indices: numpy.ndarray
        The array of index values outside the segments
    segmentlist: 
        The segment list corresponding to the selected time.
    """
    exclude, segs = indices_within_segments(times, segment_files,
                                         ifo=ifo, segment_name=segment_name)
    indices = numpy.arange(0, len(times))
    return numpy.delete(indices, exclude), segs

def get_segment_definer_comments(xml_file, include_version=True):
    """Returns a dict with the comment column as the value for each segment"""

    from pycbc_glue.ligolw.ligolw import LIGOLWContentHandler as h
    lsctables.use_in(h)

    # read segment definer table
    xmldoc, _ = ligolw_utils.load_fileobj(xml_file,
                                        gz=xml_file.name.endswith(".gz"),
                                        contenthandler=h)
    seg_def_table = table.get_table(xmldoc,
                                    lsctables.SegmentDefTable.tableName)

    # put comment column into a dict
    comment_dict = {}
    for seg_def in seg_def_table:
        if include_version:
            full_channel_name = ':'.join([str(seg_def.ifos),
                                          str(seg_def.name),
                                          str(seg_def.version)])
        else:
            full_channel_name = ':'.join([str(seg_def.ifos),
                                          str(seg_def.name)])

        comment_dict[full_channel_name] = seg_def.comment

    return comment_dict
