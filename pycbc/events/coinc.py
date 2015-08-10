# Copyright (C) 2015 Alex Nitz
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
""" This modules contains functions for calculating and manipulating
coincident triggers.
"""
import numpy, logging, h5py
from itertools import izip
from scipy.interpolate import interp1d  

def load_coincs(coinc_files):
    import pycbc.io
    class StatmapData(pycbc.io.DictArray):
        def __init__(self, data=None, segments=None, attrs=None,
                           files=None, groups=None):
            super().__init__(self, data=data, files=files, groups=groups)
            
            if data:
                self.segments=segments
                self.attrs=attrs
            elif files:
                f = h5py.File(files[0], "r")
                self.segments = f['segments']
                self.attrs = f.attrs
    
        def _return(self, data)
            return self.__class__(data=data, 
                                  attrs=self.attrs, 
                                  segments=self.segments)
    
        def cluster(self, window):
            """ Cluster the dict array, assuming it has the relevant Coinc colums,
            time1, time2, stat, and timeslide_id
            """
            interval = self.attrs['timeslide_interval']
            cid = cluster_coincs(self.stat, self.time1, self.time2,
                                     self.timeslide_id, interval, window)
            return self.select(cid) 

    columns = ['stat', 'time1', 'time2', 'trigger_id1', 'trigger_id2', 
               'template_id', 'decimation_factor', 'timeslide_id']
    return StatmapData(files=coinc_files, groups=columns)
   
def calculate_fan_map(combined_stat, dec):
    """ Return a function to map between false alarm number (FAN) and the
    combined ranking statistic.
    """
    stat_sorting = combined_stat.argsort()    
    combined_stat = combined_stat[stat_sorting]
    fan = dec[stat_sorting][::-1].cumsum()[::-1]    
    return interp1d(combined_stat, fan, fill_value=1, bounds_error=False) 

def timeslide_durations(start1, start2, end1, end2, timeslide_offsets):
    """ Find the coincident time for each timeslide.
    
    Find the coincident time for each timeslide, where the first time vector
    is slid to the right by the offset in the given timeslide_offsets 
    vector.
    
    Parameters
    ----------
    start1: numpy.ndarray
        Array of the start of valid analyzed times for detector 1
    start2: numpy.ndarray
        Array of the start of valid analyzed times for detector 2
    end1: numpy.ndarray
        Array of the end of valid analyzed times for detector 1
    end2: numpy.ndarray
        Array of the end of valid analyzed times for detector 2
    timseslide_offset: numpy.ndarray
        Array of offsets (in seconds) for each timeslide
        
    Returns
    --------
    durations: numpy.ndarray
        Array of coincident time for each timeslide in the offset array
    """
    from . import veto
    durations = []
    seg2 = veto.start_end_to_segments(start2, end2)
    for offset in timeslide_offsets:
        seg1 = veto.start_end_to_segments(start1 + offset, end1 + offset)
        durations.append(abs((seg1 & seg2).coalesce()))
    return numpy.array(durations)   
    
def time_coincidence(t1, t2, window, slide_step=0):
    """ Find coincidences by time window
    
    Parameters
    ----------
    t1 : numpy.ndarray
        Array of trigger times from the first detector
    t2 : numpy.ndarray
        Array of trigger times from the second detector
    window : float
        The coincidence window in seconds
    slide_step : optional, {None, float}
        If calculating background coincidences, the interval between background
        slides in seconds.
        
    Returns
    -------
    idx1 : numpy.ndarray
        Array of indices into the t1 array.
    idx2 : numpy.ndarray 
        Array of indices into the t2 array.
    slide : numpy.ndarray
        Array of slide ids 
    """
    if slide_step:
        fold1 = t1 % slide_step
        fold2 = t2 % slide_step
    else:
        fold1 = t1
        fold2 = t2
        
    sort1 = fold1.argsort()
    sort2 = fold2.argsort()    
    fold1 = fold1[sort1]
    fold2 = fold2[sort2]
    
    if slide_step:
        fold2 = numpy.concatenate([fold2 - slide_step, fold2, fold2 + slide_step])
        sort2 = numpy.concatenate([sort2, sort2, sort2])

    left = numpy.searchsorted(fold2, fold1 - window)
    right = numpy.searchsorted(fold2, fold1 + window)

    idx1 = numpy.repeat(sort1, right-left)
    idx2 = [sort2[l:r] for l,r in zip(left, right)]

    if len(idx2) > 0:
        idx2 = numpy.concatenate(idx2)
    else:
        idx2 = numpy.array([])
    
    if slide_step:
        diff = ((t1 / slide_step)[idx1] - (t2 / slide_step)[idx2])
        slide = numpy.rint(diff)
    else:
        slide = numpy.zeros(len(idx1))
        
    return idx1, idx2, slide


def cluster_coincs(stat, time1, time2, timeslide_id, slide, window):
    """Cluster coincident events for each timeslide separately, across 
    templates, based on the ranking statistic 

    Parameters
    ----------
    stat: numpy.ndarray
        vector of ranking values to maximize
    time1: numpy.ndarray
        first time vector
    time2: numpy.ndarray
        second time vector
    timeslide_id: numpy.ndarray
        vector that determines the timeslide offset
    slide: float
        length of the timeslides offset interval
    window: float
        length to cluster over

    Returns
    -------
    cindex: numpy.ndarray 
        The set of indices corresponding to the surviving coincidences.
    """
    
    logging.info('clustering coinc triggers over %ss window' % window)
    
    indices = []
    if numpy.isfinite(slide):
        time = (time2 + (time1 + timeslide_id * slide)) / 2
    else:
        time = 0.5 * (time2 + time1)
        
    tslide = timeslide_id.astype(numpy.float128)
    time = time.astype(numpy.float128)
    span = (time.max() - time.min()) + window * 10
    time = time + span * tslide
    
    time_sorting = time.argsort()
    stat = stat[time_sorting]
    time = time[time_sorting]
    tslide = tslide[time_sorting]
    
    logging.info('sorting...')
    left = numpy.searchsorted(time, time - window)
    right = numpy.searchsorted(time, time + window)
    logging.info('done sorting')
    indices = numpy.zeros(len(left), dtype=numpy.uint32)
    
    # i is the index we are inspecting, j is the next one to save
    i = 0
    j = 0
    while i < len(left):
        l = left[i]
        r = right[i]
        
        # If there are no other points to compare it is obviosly the max
        if (r - l) == 1:
            indices[j] = i
            j += 1
            i += 1            
            continue            
        
        # Find the location of the maximum within the time interval around i
        max_loc = stat[l:r].argmax() + l 
        
        # If this point is the max, we can skip to the right boundary
        if max_loc == i:
            indices[j] = i
            i = r
            j += 1
        
        # If the max is later than i, we can skip to it
        elif max_loc > i:
            i = max_loc 
            
        elif max_loc < i:
            i += 1

    indices = indices[:j]
            
    logging.info('done clustering coinc triggers: %s triggers remaining' % len(indices))
    return time_sorting[indices]

