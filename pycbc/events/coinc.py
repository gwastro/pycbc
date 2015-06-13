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
import numpy, logging
from itertools import izip

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
    indices = []
    for i, (l, r) in enumerate(izip(left, right)):
        if stat[l:r].argmax() + l == i:
            indices += [i]
    logging.info('done clustering coinc triggers: %s triggers remaining' % len(indices))
    return time_sorting[numpy.array(indices)]

