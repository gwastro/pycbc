""" This module contains utilities to manipulate trigger lists based on 
segment.
"""
import numpy
from glue.ligolw import ligolw, table, lsctables, utils as ligolw_utils
from glue import segments

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
    return numpy.delete(indices, exclude)
