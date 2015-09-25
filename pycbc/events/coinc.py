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
import numpy, logging, h5py, pycbc.pnutils
from itertools import izip
from scipy.interpolate import interp1d  

def background_bin_from_string(background_bins, data):
    """ Return template ids for each bin as defined by the format string
    
    Parameters
    ----------
    bins: list of strings
        List of strings which define how a background bin is taken from the
        list of templates.
    data: dict of numpy.ndarrays
        Dict with parameter key values and numpy.ndarray values which define
        the parameters of the template bank to bin up.
    
    Returns
    -------
    bins: dict
        Dictionary of location indices indexed by a bin name
    
    """
    used = numpy.array([], dtype=numpy.uint32)
    bins = {}
    for mbin in background_bins:
        name, bin_type, boundary = tuple(mbin.split(':'))
        if bin_type == 'component':
            locs = numpy.maximum(data['mass1'], data['mass2']) < float(boundary)
        elif bin_type == 'total':
            locs = data['mass1'] + data['mass2'] < float(boundary)
        elif bin_type == 'chirp':
            locs = pycbc.pnutils.mass1_mass2_to_mchirp_eta(data['mass1'], data['mass2'])[0] < float(boundary)
        elif bin_type == 'SEOBNRv2Peak':
            locs = pycbc.pnutils.get_freq('fSEOBNRv2Peak', data['mass1'],
                data['mass2'], data['spin1z'], data['spin2z']) < float(boundary)
        else:
            raise ValueError('Invalid bin type %s' % bin_type)    
        
        # make sure we don't reuse anythign from an earlier bin
        locs = numpy.where(locs)[0]
        locs = numpy.delete(locs, numpy.where(numpy.in1d(locs, used))[0])
        used = numpy.concatenate([used, locs])   
        bins[name] = locs
    return bins
           

def calculate_n_louder(bstat, fstat, dec):
    """ Calculate for each foreground event the number of background events
    that are louder than it.

    Parameters
    ----------
    bstat: numpy.ndarray
        Array of the background statistic values
    fstat: numpy.ndarray
        Array of the foreground statitsic values
    dec: numpy.ndarray
        Array of the decimation factors for the background statistics

    
    Returns
    ------- 
    cum_back_num: numpy.ndarray
        The cumulative array of background triggers 
    fore_n_louder: numpy.ndarray
        The number of background triggers above each foreground trigger
    """
    sort = bstat.argsort()
    unsort = sort.argsort()
    bstat = bstat[sort]
    dec = dec[sort]
    
    # calculate cumulative number of triggers louder than the trigger in 
    # a given index. We need to subtract the decimation factor, as the cumsum
    # includes itself in the first sum (it is inclusive of the first value)
    n_louder = dec[::-1].cumsum()[::-1] - dec
    
    # Determine how many values are louder than the foreground ones
    # We need to subtract one from the index, to be consistent with the definition
    # of n_louder, as here we do want to include the background value at the
    # found index
    fore_n_louder = n_louder[numpy.searchsorted(bstat, fstat, side='left') - 1]
    back_cum_num = n_louder[unsort]
    return back_cum_num, fore_n_louder

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
        
    return idx1.astype(numpy.uint32), idx2.astype(numpy.uint32), slide.astype(numpy.int32)


def cluster_coincs(stat, time1, time2, timeslide_id, slide, window, argmax=numpy.argmax):
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

    if len(time1) == 0 or len(time2) == 0:
        logging.info('No coinc triggers in one, or both, ifos.')
        return numpy.array([])
    
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
        max_loc = argmax(stat[l:r]) + l 
        
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

