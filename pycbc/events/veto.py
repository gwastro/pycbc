""" This module contains utilities to manipulate trigger lists based on 
segment.
"""
import numpy, urlparse, os.path
from sys import argv
from glue.ligolw import ligolw, table, lsctables, utils as ligolw_utils
from glue import segments
from glue.segments import segment, segmentlist
from glue.ligolw.lsctables import LIGOTimeGPS
from glue.ligolw.utils import segments as ligolw_segments

def start_end_to_segments(start, end):
    return segmentlist([segment(s, e) for s, e in zip(start, end)])

def segments_to_start_end(segs):
    segs.coalesce()
    return (numpy.array([s[0] for s in segs]), 
            numpy.array([s[1] for s in segs]))

def segments_to_file(segs, filename, name, ifo=""):
    """ Save segments to an xml file
    
    Parameters
    ----------
    segs : glue.segments.segmentlist
        List of segments to write to disk
    filename : str
        name of the output file
    name : 
        name of the segmentlist
        
    Returns
    -------
    File : Return a pycbc.core.File reference to the file
    """
    from pycbc.workflow.core import File

    # create XML doc and add process table
    outdoc = ligolw.Document()
    outdoc.appendChild(ligolw.LIGO_LW())
    process = ligolw_utils.process.register_to_xmldoc(outdoc, argv[0], {})

    # cast segment values into LIGOTimeGPS for glue library utils
    if type(segs[0][0]) != LIGOTimeGPS:
        fsegs = [( LIGOTimeGPS(segs[i][0]), LIGOTimeGPS(segs[i][1]) ) for i in range(len(segs))]
    else:
        fsegs = segs

    # add segments, segments summary, and segment definer tables using glue library
    with ligolw_segments.LigolwSegments(outdoc, process) as xmlsegs:
        xmlsegs.insert_from_segmentlistdict({ifo : fsegs}, name)

    # write file
    ligolw_utils.write_filename(outdoc, filename)

    # return a File instance
    url = urlparse.urlunparse(['file', 'localhost', filename, None, None, None])
    f = File(ifo, name, segs, file_url=url, tags=[name])
    f.PFN(os.path.abspath(filename), site='local')
    return f
    

def start_end_from_segments(segment_file):
    """ Return the start and end time arrays from a segment file.
    
    Parameters
    ----------
    segment_file: xml segment file
    
    Return
    ------
    start: numpy.ndarray
    end: numpy.ndarray
    """
    # dummy class needed for loading LIGOLW files
    class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
        pass
    lsctables.use_in(LIGOLWContentHandler)
    indoc = ligolw_utils.load_filename(segment_file, False, 
                                       contenthandler=LIGOLWContentHandler)
    segment_table  = table.get_table(indoc, lsctables.SegmentTable.tableName)
    start = numpy.array(segment_table.getColumnByName('start_time'))
    start_ns = numpy.array(segment_table.getColumnByName('start_time_ns'))
    end = numpy.array(segment_table.getColumnByName('end_time'))
    end_ns = numpy.array(segment_table.getColumnByName('end_time_ns'))
    return start + start_ns * 1e-9, end + end_ns * 1e-9


def indices_within_times(times, start, end):
    """ Return the an index array into times that give the values within the 
    durations defined by the start and end arrays
    
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
    tsort = times.argsort()
    times_sorted = times[tsort]
    left = numpy.searchsorted(times_sorted, start)
    right = numpy.searchsorted(times_sorted, end)
    return tsort[numpy.hstack(numpy.r_[s:e] for s, e in zip(left, right))]

def indices_outside_times(times, start, end):
    """ Return the an index array into times that give the values outside the 
    durations defined by the start and end arrays
    
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
    

def indices_within_segments(times, ifo, segment_files):
    """ Return the list of indices that should be vetoed by the segments in the
    lsit of veto_files.
    
    Parameters
    ----------
    times: numpy.ndarray of integer type
        This contains the arry of gps start times
    ifo: string
        The ifo to retrieve segments for from the segment files
    segment_files: string or list of strings
        A string or list of strings that contain the path to xml files that
        contain a segment table
        
    Returns
    --------
    indices: numpy.ndarray
        The array of index values within the segments
    segmentlist: 
        The segment list corresponding to the selected time.
    """
    # dummy class needed for loading LIGOLW files
    class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
        pass
    lsctables.use_in(LIGOLWContentHandler)

    time_sorting = numpy.argsort(times)
    times = times[time_sorting]
    indices = numpy.array([], dtype=numpy.uint32)
   
    for veto_file in segment_files:
        indoc = ligolw_utils.load_filename(veto_file, False, 
                                           contenthandler=LIGOLWContentHandler)
        segment_table  = table.get_table(indoc, 
                                           lsctables.SegmentTable.tableName)
        
        seg_def_table = table.get_table(indoc, 
                                           lsctables.SegmentDefTable.tableName)
        def_ifos = seg_def_table.getColumnByName('ifos')
        def_ids = seg_def_table.getColumnByName('segment_def_id')
        ifo_map =  {}
        for def_ifo, def_id in zip(def_ifos, def_ids):
            ifo_map[def_id] = def_ifo
        
        start = numpy.array(segment_table.getColumnByName('start_time')) + \
             numpy.array(segment_table.getColumnByName('start_time_ns')) * 1e-9
        end = numpy.array(segment_table.getColumnByName('end_time')) + \
             numpy.array(segment_table.getColumnByName('end_time_ns')) * 1e-9
        ifos = [ifo_map[v] for v in segment_table.getColumnByName('segment_def_id')]
        
        veto_segs = segments.segmentlist()
        for s, e, ifo_row in zip(start, end, ifos):
            if ifo != ifo_row:
                continue            
            veto_segs += [segments.segment(s, e)]

        veto_segs.coalesce()        

        left = numpy.searchsorted(times, start, side='left')
        right = numpy.searchsorted(times, end, side='right')
        for li, ri, ifo_row in zip(left, right, ifos):
            if ifo != ifo_row:
                continue
                
            seg_indices = numpy.arange(li, ri, 1).astype(numpy.uint32)
            indices=numpy.union1d(seg_indices, indices)  
    return time_sorting[indices], veto_segs
 
def indices_outside_segments(times, ifo, segment_files):
    """ Return the list of indices that are outside the segments in the
    list of segment files.
    
    Parameters
    ----------
    times: numpy.ndarray of integer type
        This contains the arry of gps start times
    ifo: string
        The ifo to retrieve segments for from the segment files
    segment_files: string or list of strings
        A string or list of strings that contain the path to xml files that
        contain a segment table
        
    Returns
    --------
    indices: numpy.ndarray
        The array of index values outside the segments
    segmentlist: 
        The segment list corresponding to the selected time.
    """
    exclude, segs = indices_within_segments(times, ifo, segment_files)
    indices = numpy.arange(0, len(times))
    return numpy.delete(indices, exclude), segs
    
