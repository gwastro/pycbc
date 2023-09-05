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
""" This module contains functions for calculating and manipulating
coincident triggers.
"""

import numpy, logging, pycbc.pnutils, pycbc.conversions, copy, lal
from pycbc.detector import Detector, ppdets
from .eventmgr_cython import coincbuffer_expireelements
from .eventmgr_cython import coincbuffer_numgreater
from .eventmgr_cython import timecoincidence_constructidxs
from .eventmgr_cython import timecoincidence_constructfold
from .eventmgr_cython import timecoincidence_getslideint
from .eventmgr_cython import timecoincidence_findidxlen
from .eventmgr_cython import timecluster_cython

# Mapping used in background_bin_from_string to select approximant for
# duration function, if duration-based binning is used.
_APPROXIMANT_DURATION_MAP = {
    'SEOBNRv2duration': 'SEOBNRv2',
    'SEOBNRv4duration': 'SEOBNRv4',
    'SEOBNRv5duration': 'SEOBNRv5_ROM'
}


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
    # Some duration/peak frequency functions are expensive.
    # Do not want to recompute many times, if using lots of bins.
    cached_values = {}
    for mbin in background_bins:
        locs = None
        name, bin_type_list, boundary_list = tuple(mbin.split(':'))

        bin_type_list = bin_type_list.split(',')
        boundary_list = boundary_list.split(',')

        for bin_type, boundary in zip(bin_type_list, boundary_list):
            if boundary[0:2] == 'lt':
                member_func = lambda vals, bd=boundary : vals < float(bd[2:])
            elif boundary[0:2] == 'gt':
                member_func = lambda vals, bd=boundary : vals > float(bd[2:])
            else:
                raise RuntimeError("Can't parse boundary condition! Must begin "
                                   "with 'lt' or 'gt'")

            if bin_type in cached_values:
                vals = cached_values[bin_type]
            elif bin_type == 'component' and boundary[0:2] == 'lt':
                # maximum component mass is less than boundary value
                vals = numpy.maximum(data['mass1'], data['mass2'])
            elif bin_type == 'component' and boundary[0:2] == 'gt':
                # minimum component mass is greater than bdary
                vals = numpy.minimum(data['mass1'], data['mass2'])
            elif bin_type == 'total':
                vals = data['mass1'] + data['mass2']
            elif bin_type == 'chirp':
                vals = pycbc.pnutils.mass1_mass2_to_mchirp_eta(
                                                   data['mass1'], data['mass2'])[0]
            elif bin_type == 'ratio':
                vals = pycbc.conversions.q_from_mass1_mass2(
                                                   data['mass1'], data['mass2'])
            elif bin_type == 'eta':
                vals = pycbc.pnutils.mass1_mass2_to_mchirp_eta(
                                                   data['mass1'], data['mass2'])[1]
            elif bin_type == 'chi_eff':
                vals = pycbc.conversions.chi_eff(data['mass1'], data['mass2'],
                                                 data['spin1z'], data['spin2z'])
            elif bin_type in ['SEOBNRv2Peak', 'SEOBNRv4Peak']:
                vals = pycbc.pnutils.get_freq(
                    'f' + bin_type,
                    data['mass1'],
                    data['mass2'],
                    data['spin1z'],
                    data['spin2z']
                )
                cached_values[bin_type] = vals
            elif bin_type in _APPROXIMANT_DURATION_MAP:
                vals = pycbc.pnutils.get_imr_duration(
                    data['mass1'],
                    data['mass2'],
                    data['spin1z'],
                    data['spin2z'],
                    data['f_lower'],
                    approximant=_APPROXIMANT_DURATION_MAP[bin_type]
                )
                cached_values[bin_type] = vals
            else:
                raise ValueError('Invalid bin type %s' % bin_type)

            sub_locs = member_func(vals)
            sub_locs = numpy.where(sub_locs)[0]
            if locs is not None:
                # find intersection of boundary conditions
                locs = numpy.intersect1d(locs, sub_locs)
            else:
                locs = sub_locs

        # make sure we don't reuse anything from an earlier bin
        locs = numpy.delete(locs, numpy.where(numpy.in1d(locs, used))[0])
        used = numpy.concatenate([used, locs])
        bins[name] = locs

    return bins


def timeslide_durations(start1, start2, end1, end2, timeslide_offsets):
    """ Find the coincident time for each timeslide.

    Find the coincident time for each timeslide, where the first time vector
    is slid to the right by the offset in the given timeslide_offsets vector.

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
        Coincidence window maximum time difference, arbitrary units (usually s)
    slide_step : float (default 0)
        If calculating background coincidences, the interval between background
        slides, arbitrary units (usually s)

    Returns
    -------
    idx1 : numpy.ndarray
        Array of indices into the t1 array for coincident triggers
    idx2 : numpy.ndarray
        Array of indices into the t2 array
    slide : numpy.ndarray
        Array of slide ids
    """
    if slide_step:
        length1 = len(t1)
        length2 = len(t2)
        fold1 = numpy.zeros(length1, dtype=numpy.float64)
        fold2 = numpy.zeros(length2, dtype=numpy.float64)
        timecoincidence_constructfold(fold1, fold2, t1, t2, slide_step,
                                      length1, length2)
    else:
        fold1 = t1
        fold2 = t2

    sort1 = fold1.argsort()
    sort2 = fold2.argsort()
    fold1 = fold1[sort1]
    fold2 = fold2[sort2]

    if slide_step:
        # FIXME explain this
        fold2 = numpy.concatenate([fold2 - slide_step, fold2,
                                   fold2 + slide_step])

    left = fold2.searchsorted(fold1 - window)
    right = fold2.searchsorted(fold1 + window)

    lenidx = timecoincidence_findidxlen(left, right, len(left))
    idx1 = numpy.zeros(lenidx, dtype=numpy.uint32)
    idx2 = numpy.zeros(lenidx, dtype=numpy.uint32)
    timecoincidence_constructidxs(idx1, idx2, sort1, sort2, left, right,
                                  len(left), len(sort2))

    slide = numpy.zeros(lenidx, dtype=numpy.int32)
    if slide_step:
        timecoincidence_getslideint(slide, t1, t2, idx1, idx2, slide_step)
    else:
        slide = numpy.zeros(len(idx1))

    return idx1, idx2, slide


def time_multi_coincidence(times, slide_step=0, slop=.003,
                           pivot='H1', fixed='L1'):
    """ Find multi detector coincidences.

    Parameters
    ----------
    times: dict of numpy.ndarrays
        Dictionary keyed by ifo of single ifo trigger times
    slide_step: float
        Interval between time slides
    slop: float
        The amount of time to add to the TOF between detectors for coincidence
    pivot: str
        The ifo to which time shifts are applied in first stage coincidence
    fixed: str
        The other ifo used in first stage coincidence, subsequently used as a
        time reference for additional ifos. All other ifos are not time shifted
        relative to this ifo

    Returns
    -------
    ids: dict of arrays of int
        Dictionary keyed by ifo with ids of trigger times forming coincidences.
        Coincidence is tested for every pair of ifos that can be formed from
        the input dict: only those tuples of times passing all tests are
        recorded
    slide: array of int
        Slide ids of coincident triggers in pivot ifo
    """
    def win(ifo1, ifo2):
        d1 = Detector(ifo1)
        d2 = Detector(ifo2)
        return d1.light_travel_time_to_detector(d2) + slop

    # Find coincs between the 'pivot' and 'fixed' detectors as in 2-ifo case
    pivot_id, fix_id, slide = time_coincidence(times[pivot], times[fixed],
                                               win(pivot, fixed),
                                               slide_step=slide_step)

    # Additional detectors do not slide independently of the 'fixed' one
    # Each trigger in an additional detector must be concident with both
    # triggers in an existing coincidence

    # Slide 'pivot' trigger times to be coincident with trigger times in
    # 'fixed' detector
    fixed_time = times[fixed][fix_id]
    pivot_time = times[pivot][pivot_id] - slide_step * slide
    ctimes = {fixed: fixed_time, pivot: pivot_time}
    ids = {fixed: fix_id, pivot: pivot_id}

    dep_ifos = [ifo for ifo in times.keys() if ifo != fixed and ifo != pivot]
    for ifo1 in dep_ifos:
        # FIXME - make this loop into a function?

        # otime is extra ifo time in original trigger order
        otime = times[ifo1]
        # tsort gives ordering from original order to time sorted order
        tsort = otime.argsort()
        time1 = otime[tsort]

        # Find coincidences between dependent ifo triggers and existing coincs
        # - Cycle over fixed and pivot
        # - At the 1st iteration, the fixed and pivot triggers are reduced to
        #  those for which the first out of fixed/pivot forms a coinc with ifo1
        # - At the 2nd iteration, we are left with triggers for which both
        #  fixed and pivot are coincident with ifo1
        # - If there is more than 1 dependent ifo, ones that were previously
        #  tested against fixed and pivot are now present for testing with new
        #  dependent ifos
        for ifo2 in ids:
            logging.info('added ifo %s, testing against %s' % (ifo1, ifo2))
            w = win(ifo1, ifo2)
            left = time1.searchsorted(ctimes[ifo2] - w)
            right = time1.searchsorted(ctimes[ifo2] + w)
            # Any times within time1 coincident with the time in ifo2 have
            # indices between 'left' and 'right'
            # 'nz' indexes into times in ifo2 which have coincidences with ifo1
            # times
            nz = (right - left).nonzero()
            if len(right - left):
                rlmax = (right - left).max()
            if len(nz[0]) and rlmax > 1:
                # We expect at most one coincident time in ifo1, assuming
                #  trigger spacing in ifo1 > time window.
                # However there are rare corner cases at starts/ends of inspiral
                #  jobs. For these, arbitrarily keep the first trigger and
                #  discard the second (and any subsequent ones).
                logging.warning('Triggers in %s are closer than coincidence '
                                'window, 1 or more coincs will be discarded. '
                                'This is a warning, not an error.' % ifo1)
            # identify indices of times in ifo1 that form coincs with ifo2
            dep_ids = left[nz]
            # slide is array of slide ids attached to pivot ifo
            slide = slide[nz]

            for ifo in ctimes:
                # cycle over fixed and pivot & any previous additional ifos
                # reduce times and IDs to just those forming a coinc with ifo1
                ctimes[ifo] = ctimes[ifo][nz]
                ids[ifo] = ids[ifo][nz]

        # undo time sorting on indices of ifo1 triggers, add ifo1 ids and times
        # to dicts for testing against any additional detectrs
        ids[ifo1] = tsort[dep_ids]
        ctimes[ifo1] = otime[ids[ifo1]]

    return ids, slide


def cluster_coincs(stat, time1, time2, timeslide_id, slide, window, **kwargs):
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
    if len(time1) == 0 or len(time2) == 0:
        logging.info('No coinc triggers in one, or both, ifos.')
        return numpy.array([])

    if numpy.isfinite(slide):
        # for a time shifted coinc, time1 is greater than time2 by approximately timeslide_id*slide
        # adding this quantity gives a mean coinc time located around time1
        time = (time1 + time2 + timeslide_id * slide) / 2
    else:
        time = 0.5 * (time2 + time1)

    tslide = timeslide_id.astype(numpy.longdouble)
    time = time.astype(numpy.longdouble)

    span = (time.max() - time.min()) + window * 10
    time = time + span * tslide
    logging.info('Clustering events over %s s window', window)
    cidx = cluster_over_time(stat, time, window, **kwargs)
    logging.info('%d triggers remaining', len(cidx))
    return cidx


def cluster_coincs_multiifo(stat, time_coincs, timeslide_id, slide, window,
                            **kwargs):
    """Cluster coincident events for each timeslide separately, across
    templates, based on the ranking statistic

    Parameters
    ----------
    stat: numpy.ndarray
        vector of ranking values to maximize
    time_coincs: tuple of numpy.ndarrays
        trigger times for each ifo, or -1 if an ifo does not participate in a coinc
    timeslide_id: numpy.ndarray
        vector that determines the timeslide offset
    slide: float
        length of the timeslides offset interval
    window: float
        duration of clustering window in seconds

    Returns
    -------
    cindex: numpy.ndarray
        The set of indices corresponding to the surviving coincidences
    """
    time_coinc_zip = list(zip(*time_coincs))
    if len(time_coinc_zip) == 0:
        logging.info('No coincident triggers.')
        return numpy.array([])

    time_avg_num = []
    #find number of ifos and mean time over participating ifos for each coinc
    for tc in time_coinc_zip:
        time_avg_num.append(mean_if_greater_than_zero(tc))

    time_avg, num_ifos = zip(*time_avg_num)

    time_avg = numpy.array(time_avg)
    num_ifos = numpy.array(num_ifos)

    # shift all but the pivot ifo by (num_ifos-1) * timeslide_id * slide
    # this leads to a mean coinc time located around pivot time
    if numpy.isfinite(slide):
        nifos_minusone = (num_ifos - numpy.ones_like(num_ifos))
        time_avg = time_avg + (nifos_minusone * timeslide_id * slide)/num_ifos

    tslide = timeslide_id.astype(numpy.longdouble)
    time_avg = time_avg.astype(numpy.longdouble)

    span = (time_avg.max() - time_avg.min()) + window * 10
    time_avg = time_avg + span * tslide
    logging.info('Clustering events over %s s window', window)
    cidx = cluster_over_time(stat, time_avg, window, **kwargs)
    logging.info('%d triggers remaining', len(cidx))

    return cidx


def mean_if_greater_than_zero(vals):
    """ Calculate mean over numerical values, ignoring values less than zero.
    E.g. used for mean time over coincident triggers when timestamps are set
    to -1 for ifos not included in the coincidence.

    Parameters
    ----------
    vals: iterator of numerical values
        values to be mean averaged

    Returns
    -------
    mean: float
        The mean of the values in the original vector which are
        greater than zero
    num_above_zero: int
        The number of entries in the vector which are above zero
    """
    vals = numpy.array(vals)
    above_zero = vals > 0
    return vals[above_zero].mean(), above_zero.sum()


def cluster_over_time(stat, time, window, method='python',
                      argmax=numpy.argmax):
    """Cluster generalized transient events over time via maximum stat over a
    symmetric sliding window

    Parameters
    ----------
    stat: numpy.ndarray
        vector of ranking values to maximize
    time: numpy.ndarray
        time to use for clustering
    window: float
        length to cluster over
    method: string
        Either "cython" to use the cython implementation, or "python" to use
        the pure python version.
    argmax: function
        the function used to calculate the maximum value

    Returns
    -------
    cindex: numpy.ndarray
        The set of indices corresponding to the surviving coincidences.
    """

    indices = []
    time_sorting = time.argsort()
    stat = stat[time_sorting]
    time = time[time_sorting]

    left = time.searchsorted(time - window)
    right = time.searchsorted(time + window)
    indices = numpy.zeros(len(left), dtype=numpy.uint32)

    logging.debug('%d triggers before clustering', len(time))

    if method == 'cython':
        j = timecluster_cython(indices, left, right, stat, len(left))
    elif method == 'python':
        # i is the index we are inspecting, j is the next one to save
        i = 0
        j = 0
        while i < len(left):
            l = left[i]
            r = right[i]

            # If there are no other points to compare it is obviously the max
            if (r - l) == 1:
                indices[j] = i
                j += 1
                i += 1
                continue

            # Find the location of the maximum within the time interval
            # around i
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
    else:
        raise ValueError(f'Do not recognize method {method}')

    indices = indices[:j]

    logging.debug('%d triggers remaining', len(indices))
    return time_sorting[indices]


class MultiRingBuffer(object):
    """Dynamic size n-dimensional ring buffer that can expire elements."""

    def __init__(self, num_rings, max_time, dtype, min_buffer_size=16,
                 buffer_increment=8, resize_invalid_fraction=0.4):
        """
        Parameters
        ----------
        num_rings: int
            The number of ring buffers to create. They all will have the same
            intrinsic size and will expire at the same time.
        max_time: int
            The maximum "time" an element can exist in each ring.
        dtype: numpy.dtype
            The type of each element in the ring buffer.
        min_buffer_size: int (optional: default=16)
            All ring buffers will be initialized to this length. If a buffer is
            made larger it will no smaller than this value. Buffers may become
            smaller than this length at any given time as triggers are expired.
        buffer_increment: int (optional: default=8)
            When increasing ring buffers, add this many points. Be careful if
            changing this and min_buffer_size from default values, it is
            possible to get stuck in a mode where the buffers are always being
            resized.
        resize_invalid_fraction: float (optional:default=0.4)
            If this fraction of any buffer contains unused data points then
            resize it to contain only valid points. As with the previous two
            options, be careful changing default values, it is
            possible to get stuck in a mode where the buffers are always being
            resized.
        """
        self.max_time = max_time
        self.buffer = []
        self.buffer_expire = []
        self.valid_ends = []
        self.valid_starts = []
        self.min_buffer_size = min_buffer_size
        self.buffer_increment = buffer_increment
        self.resize_invalid_fraction = resize_invalid_fraction
        for _ in range(num_rings):
            self.buffer.append(numpy.zeros(self.min_buffer_size, dtype=dtype))
            self.buffer_expire.append(numpy.zeros(self.min_buffer_size,
                                                  dtype=int))
            self.valid_ends.append(0)
            self.valid_starts.append(0)
        self.time = 0

    @property
    def filled_time(self):
        return min(self.time, self.max_time)

    def num_elements(self):
        count = 0
        for idx, a in enumerate(self.buffer):
            vals = self.valid_starts[idx]
            vale = self.valid_ends[idx]
            count += len(a[vals:vale])
        return count

    @property
    def nbytes(self):
        return sum([a.nbytes for a in self.buffer])

    def discard_last(self, indices):
        """Discard the triggers added in the latest update"""
        for i in indices:
            self.valid_ends[i] -= 1

    def advance_time(self):
        """Advance the internal time increment by 1, expiring any triggers
        that are now too old.
        """
        self.time += 1

    def add(self, indices, values):
        """Add triggers in 'values' to the buffers indicated by the indices
        """
        for i, v in zip(indices, values):
            # Expand ring buffer size if needed
            if self.valid_ends[i] == len(self.buffer[i]):
                # First clear out any old triggers before resizing
                self.update_valid_start(i)
                self.check_expired_triggers(i)

                # Then increase arrays by buffer_increment
                self.buffer[i] = numpy.resize(
                    self.buffer[i],
                    max(
                        len(self.buffer[i]) + self.buffer_increment,
                        self.min_buffer_size
                    )
                )
                self.buffer_expire[i] = numpy.resize(
                    self.buffer_expire[i],
                    max(
                        len(self.buffer[i]) + self.buffer_increment,
                        self.min_buffer_size
                    )
                )
            curr_pos = self.valid_ends[i]
            self.buffer[i][curr_pos] = v
            self.buffer_expire[i][curr_pos] = self.time
            self.valid_ends[i] = self.valid_ends[i] + 1
        self.advance_time()

    def valid_slice(self, buffer_index):
        """Return the valid slice for this buffer index"""
        ret_slice = slice(
            self.valid_starts[buffer_index],
            self.valid_ends[buffer_index]
        )
        return ret_slice

    def expire_vector(self, buffer_index):
        """Return the expiration vector of a given ring buffer """
        return self.buffer_expire[buffer_index][self.valid_slice(buffer_index)]

    def update_valid_start(self, buffer_index):
        """Update the valid_start for the given buffer index"""
        expired = self.time - self.max_time
        exp = self.buffer_expire[buffer_index]
        j = self.valid_starts[buffer_index]
        while j < self.valid_ends[buffer_index]:
            # Everything before this j must be expired
            if exp[j] >= expired:
                break
            j += 1
        self.valid_starts[buffer_index] = j

    def check_expired_triggers(self, buffer_index):
        """Check if we should free memory for this buffer index.

        Check what fraction of triggers are expired in the specified buffer
        and if it is more than the allowed fraction (set by
        self.resize_invalid_fraction) resize the array to remove them.
        """
        val_start = self.valid_starts[buffer_index]
        val_end = self.valid_ends[buffer_index]
        buf_len = len(self.buffer[buffer_index])
        invalid_limit = self.resize_invalid_fraction * buf_len
        if (buf_len - val_end) + val_start > invalid_limit:
            # If self.resize_invalid_fraction of stored triggers are expired
            # or are not set, free up memory
            self.buffer_expire[buffer_index] = self.buffer_expire[buffer_index][val_start:val_end].copy()
            self.buffer[buffer_index] = self.buffer[buffer_index][val_start:val_end].copy()
            self.valid_ends[buffer_index] -= val_start
            self.valid_starts[buffer_index] = 0

    def data(self, buffer_index):
        """Return the data vector for a given ring buffer"""
        self.update_valid_start(buffer_index)
        self.check_expired_triggers(buffer_index)

        return self.buffer[buffer_index][self.valid_slice(buffer_index)]


class CoincExpireBuffer(object):
    """Unordered dynamic sized buffer that handles
    multiple expiration vectors.
    """

    def __init__(self, expiration, ifos,
                       initial_size=2**20, dtype=numpy.float32):
        """
        Parameters
        ----------
        expiration: int
            The 'time' in arbitrary integer units to allow to pass before
            removing an element.
        ifos: list of strs
            List of strings to identify the multiple data expiration times.
        initial_size: int, optional
            The initial size of the buffer.
        dtype: numpy.dtype
            The dtype of each element of the buffer.
        """

        self.expiration = expiration
        self.buffer = numpy.zeros(initial_size, dtype=dtype)
        self.index = 0
        self.ifos = ifos

        self.time = {}
        self.timer = {}
        for ifo in self.ifos:
            self.time[ifo] = 0
            self.timer[ifo] = numpy.zeros(initial_size, dtype=numpy.int32)

    def __len__(self):
        return self.index

    @property
    def nbytes(self):
        """Returns the approximate memory usage of self.
        """
        nbs = [self.timer[ifo].nbytes for ifo in self.ifos]
        nbs.append(self.buffer.nbytes)
        return sum(nbs)

    def increment(self, ifos):
        """Increment without adding triggers"""
        self.add([], [], ifos)

    def remove(self, num):
        """Remove the the last 'num' elements from the buffer"""
        self.index -= num

    def add(self, values, times, ifos):
        """Add values to the internal buffer

        Parameters
        ----------
        values: numpy.ndarray
            Array of elements to add to the internal buffer.
        times: dict of arrays
            The current time to use for each element being added.
        ifos: list of strs
            The set of timers to be incremented.
        """

        for ifo in ifos:
            self.time[ifo] += 1

        # Resize the internal buffer if we need more space
        if self.index + len(values) >= len(self.buffer):
            newlen = len(self.buffer) * 2
            for ifo in self.ifos:
                self.timer[ifo].resize(newlen)
            self.buffer.resize(newlen, refcheck=False)

        self.buffer[self.index:self.index+len(values)] = values
        if len(values) > 0:
            for ifo in self.ifos:
                self.timer[ifo][self.index:self.index+len(values)] = times[ifo]

            self.index += len(values)

        # Remove the expired old elements
        if len(ifos) == 2:
            # Cython version for two ifo case
            self.index = coincbuffer_expireelements(
                self.buffer,
                self.timer[ifos[0]],
                self.timer[ifos[1]],
                self.time[ifos[0]],
                self.time[ifos[1]],
                self.expiration,
                self.index
            )
        else:
            # Numpy version for >2 ifo case
            keep = None
            for ifo in ifos:
                kt = self.timer[ifo][:self.index] >= self.time[ifo] - self.expiration
                keep = numpy.logical_and(keep, kt) if keep is not None else kt

            self.buffer[:keep.sum()] = self.buffer[:self.index][keep]
            for ifo in self.ifos:
                self.timer[ifo][:keep.sum()] = self.timer[ifo][:self.index][keep]
            self.index = keep.sum()

    def num_greater(self, value):
        """Return the number of elements larger than 'value'"""
        return coincbuffer_numgreater(self.buffer, self.index, value)

    @property
    def data(self):
        """Return the array of elements"""
        return self.buffer[:self.index]


class LiveCoincTimeslideBackgroundEstimator(object):
    """Rolling buffer background estimation."""

    def __init__(self, num_templates, analysis_block, background_statistic,
                 sngl_ranking, stat_files, ifos,
                 ifar_limit=100,
                 timeslide_interval=.035,
                 coinc_threshold=.002,
                 return_background=False,
                 **kwargs):
        """
        Parameters
        ----------
        num_templates: int
            The size of the template bank
        analysis_block: int
            The number of seconds in each analysis segment
        background_statistic: str
            The name of the statistic to rank coincident events.
        sngl_ranking: str
            The single detector ranking to use with the background statistic
        stat_files: list of strs
            List of filenames that contain information used to construct
            various coincident statistics.
        ifos: list of strs
            List of ifo names that are being analyzed. At the moment this must
            be two items such as ['H1', 'L1'].
        ifar_limit: float
            The largest inverse false alarm rate in years that we would like to
            calculate.
        timeslide_interval: float
            The time in seconds between consecutive timeslide offsets.
        coinc_threshold: float
            Amount of time allowed to form a coincidence in addition to the
            time of flight in seconds.
        return_background: boolean
            If true, background triggers will also be included in the file
            output.
        kwargs: dict
            Additional options for the statistic to use. See stat.py
            for more details on statistic options.
        """
        from . import stat
        self.num_templates = num_templates
        self.analysis_block = analysis_block

        stat_class = stat.get_statistic(background_statistic)
        self.stat_calculator = stat_class(
            sngl_ranking,
            stat_files,
            ifos=ifos,
            **kwargs
        )

        self.timeslide_interval = timeslide_interval
        self.return_background = return_background

        self.ifos = ifos
        if len(self.ifos) != 2:
            raise ValueError("Only a two ifo analysis is supported at this time")

        self.lookback_time = (ifar_limit * lal.YRJUL_SI * timeslide_interval) ** 0.5
        self.buffer_size = int(numpy.ceil(self.lookback_time / analysis_block))

        det0, det1 = Detector(ifos[0]), Detector(ifos[1])
        self.time_window = det0.light_travel_time_to_detector(det1) + coinc_threshold
        self.coincs = CoincExpireBuffer(self.buffer_size, self.ifos)

        self.singles = {}

        # temporary array used in `_find_coincs()` to turn `trig_stat`
        # into an array much faster than using `numpy.resize()`
        self.trig_stat_memory = None

    @classmethod
    def pick_best_coinc(cls, coinc_results):
        """Choose the best two-ifo coinc by ifar first, then statistic if needed.

        This function picks which of the available double-ifo coincs to use.
        It chooses the best (highest) ifar. The ranking statistic is used as
        a tie-breaker.
        A trials factor is applied if multiple types of coincs are possible
        at this time given the active ifos.

        Parameters
        ----------
        coinc_results: list of coinc result dicts
            Dictionary by detector pair of coinc result dicts.

        Returns
        -------
        best: coinc results dict
            If there is a coinc, this will contain the 'best' one. Otherwise
            it will return the provided dict.
        """
        mstat = 0
        mifar = 0
        mresult = None

        # record the trials factor from the possible coincs we could
        # maximize over
        trials = 0
        for result in coinc_results:
            # Check that a coinc was possible. See the 'add_singles' method
            # to see where this flag was added into the results dict
            if 'coinc_possible' in result:
                trials += 1

                # Check that a coinc exists
                if 'foreground/ifar' in result:
                    ifar = result['foreground/ifar']
                    stat = result['foreground/stat']
                    if ifar > mifar or (ifar == mifar and stat > mstat):
                        mifar = ifar
                        mstat = stat
                        mresult = result

        # apply trials factor for the best coinc
        if mresult:
            mresult['foreground/ifar'] = mifar / float(trials)
            logging.info('Found %s coinc with ifar %s',
                         mresult['foreground/type'],
                         mresult['foreground/ifar'])
            return mresult
        # If no coinc, just return one of the results dictionaries. They will
        # all contain the same results (i.e. single triggers) in this case.
        else:
            return coinc_results[0]

    @classmethod
    def from_cli(cls, args, num_templates, analysis_chunk, ifos):
        from . import stat

        # Allow None inputs
        stat_files = args.statistic_files or []
        stat_keywords = args.statistic_keywords or []

        # flatten the list of lists of filenames to a single list (may be empty)
        stat_files = sum(stat_files, [])

        kwargs = stat.parse_statistic_keywords_opt(stat_keywords)

        return cls(num_templates, analysis_chunk,
                   args.ranking_statistic,
                   args.sngl_ranking,
                   stat_files,
                   return_background=args.store_background,
                   ifar_limit=args.background_ifar_limit,
                   timeslide_interval=args.timeslide_interval,
                   ifos=ifos,
                   **kwargs)

    @staticmethod
    def insert_args(parser):
        from . import stat

        stat.insert_statistic_option_group(parser)

        group = parser.add_argument_group('Coincident Background Estimation')
        group.add_argument('--store-background', action='store_true',
            help="Return background triggers with zerolag coincidencs")
        group.add_argument('--background-ifar-limit', type=float,
            help="The limit on inverse false alarm rate to calculate "
                 "background in years", default=100.0)
        group.add_argument('--timeslide-interval', type=float,
            help="The interval between timeslides in seconds", default=0.1)
        group.add_argument('--ifar-remove-threshold', type=float,
            help="NOT YET IMPLEMENTED", default=100.0)

    @property
    def background_time(self):
        """Return the amount of background time that the buffers contain"""
        time = 1.0 / self.timeslide_interval
        for ifo in self.singles:
            time *= self.singles[ifo].filled_time * self.analysis_block
        return time

    def save_state(self, filename):
        """Save the current state of the background buffers"""
        import pickle
        pickle.dump(self, filename)

    @staticmethod
    def restore_state(filename):
        """Restore state of the background buffers from a file"""
        import pickle
        return pickle.load(filename)

    def ifar(self, coinc_stat):
        """Map a given value of the coincident ranking statistic to an inverse
        false-alarm rate (IFAR) using the interally stored background sample.

        Parameters
        ----------
        coinc_stat: float
            Value of the coincident ranking statistic to be converted.

        Returns
        -------
        ifar: float
            Inverse false-alarm rate in unit of years.
        ifar_saturated: bool
            True if `coinc_stat` is larger than all the available background,
            in which case `ifar` is to be considered an upper limit.
        """
        n = self.coincs.num_greater(coinc_stat)
        ifar = self.background_time / lal.YRJUL_SI / (n + 1)
        return ifar, n == 0

    def set_singles_buffer(self, results):
        """Create the singles buffer

        This creates the singles buffer for each ifo. The dtype is determined
        by a representative sample of the single triggers in the results.

        Parameters
        ----------
        results: dict of dict
            Dict indexed by ifo and then trigger column.
        """
        # Determine the dtype from a sample of the data.
        self.singles_dtype = []
        data = False
        for ifo in self.ifos:
            if ifo in results and results[ifo] is not False \
                    and len(results[ifo]['snr']):
                data = results[ifo]
                break

        if data is False:
            return

        for key in data:
            self.singles_dtype.append((key, data[key].dtype))

        if 'stat' not in data:
            self.singles_dtype.append(('stat', self.stat_calculator.single_dtype))

        # Create a ring buffer for each template ifo combination
        for ifo in self.ifos:
            self.singles[ifo] = MultiRingBuffer(self.num_templates,
                                            self.buffer_size,
                                            self.singles_dtype)

    def _add_singles_to_buffer(self, results, ifos):
        """Add single detector triggers to the internal buffer

        Parameters
        ----------
        results: dict
            Dictionary of dictionaries indexed by ifo and keys such as 'snr',
            'chisq', etc. The specific format is determined by the
            LiveBatchMatchedFilter class.

        Returns
        -------
        updated_singles: dict of numpy.ndarrays
            Array of indices that have been just updated in the internal
            buffers of single detector triggers.
        """
        if len(self.singles.keys()) == 0:
            self.set_singles_buffer(results)
        # If this *still* didn't work, no triggers in first set, try next time
        if len(self.singles.keys()) == 0:
            return {}

        # convert to single detector trigger values
        # FIXME Currently configured to use pycbc live output
        # where chisq is the reduced chisq and chisq_dof is the actual DOF
        logging.info("adding singles to the background estimate...")
        updated_indices = {}
        for ifo in ifos:
            trigs = results[ifo]

            if len(trigs['snr'] > 0):
                trigsc = copy.copy(trigs)
                trigsc['chisq'] = trigs['chisq'] * trigs['chisq_dof']
                trigsc['chisq_dof'] = (trigs['chisq_dof'] + 2) / 2
                single_stat = self.stat_calculator.single(trigsc)
            else:
                single_stat = numpy.array([], ndmin=1,
                              dtype=self.stat_calculator.single_dtype)
            trigs['stat'] = single_stat

            # add each single detector trigger to the and advance the buffer
            data = numpy.zeros(len(single_stat), dtype=self.singles_dtype)
            for key, value in trigs.items():
                data[key] = value

            self.singles[ifo].add(trigs['template_id'], data)
            updated_indices[ifo] = trigs['template_id']
        return updated_indices

    def _find_coincs(self, results, valid_ifos):
        """Look for coincs within the set of single triggers

        Parameters
        ----------
        results: dict
            Dictionary of dictionaries indexed by ifo and keys such as 'snr',
            'chisq', etc. The specific format is determined by the
            LiveBatchMatchedFilter class.
        valid_ifos: list of strs
            List of ifos for which new triggers might exist. This must be a
            subset of self.ifos. If an ifo is in self.ifos but not in this list
            either the ifo is down, or its data has been flagged as "bad".

        Returns
        -------
        num_background: int
            Number of time shifted coincidences found.
        coinc_results: dict of arrays
            A dictionary of arrays containing the coincident results.
        """
        # For each new single detector trigger find the allowed coincidences
        # Record the template and the index of the single trigger that forms
        # each coincidence

        # Initialize
        cstat = [[]]
        offsets = []
        ctimes = {self.ifos[0]:[], self.ifos[1]:[]}
        single_expire = {self.ifos[0]:[], self.ifos[1]:[]}
        template_ids = [[]]
        trigger_ids = {self.ifos[0]:[[]], self.ifos[1]:[[]]}

        # Calculate all the permutations of coincident triggers for each
        # new single detector trigger collected
        # Currently only two detectors are supported.
        # For each ifo, check its newly added triggers for (zerolag and time
        # shift) coincs with all currently stored triggers in the other ifo.
        # Do this by keeping the ifo with new triggers fixed and time shifting
        # the other ifo. The list 'shift_vec' must be in the same order as
        # self.ifos and contain -1 for the shift_ifo / 0 for the fixed_ifo.
        for fixed_ifo, shift_ifo, shift_vec in zip(
            [self.ifos[0], self.ifos[1]],
            [self.ifos[1], self.ifos[0]],
            [[0, -1], [-1, 0]]
        ):
            if fixed_ifo not in valid_ifos:
                # This ifo is not online now, so no new triggers or coincs
                continue
            # Find newly added triggers in fixed_ifo
            trigs = results[fixed_ifo]
            # Loop over them one trigger at a time
            for i in range(len(trigs['end_time'])):
                trig_stat = trigs['stat'][i]
                trig_time = trigs['end_time'][i]
                template = trigs['template_id'][i]

                # Get current shift_ifo triggers in the same template
                times = self.singles[shift_ifo].data(template)['end_time']
                stats = self.singles[shift_ifo].data(template)['stat']

                # Perform coincidence. i1 is the list of trigger indices in the
                # shift_ifo which make coincs, slide is the corresponding slide
                # index.
                # (The second output would just be a list of zeroes as we only
                # have one trigger in the fixed_ifo.)
                i1, _, slide = time_coincidence(times,
                                 numpy.array(trig_time, ndmin=1,
                                 dtype=numpy.float64),
                                 self.time_window,
                                 self.timeslide_interval)

                # Make a copy of the fixed ifo trig_stat for each coinc.
                # NB for some statistics the "stat" entry holds more than just
                # a ranking number. E.g. for the phase time consistency test,
                # it must also contain the phase, time and sensitivity.
                if self.trig_stat_memory is None:
                    self.trig_stat_memory = numpy.zeros(
                        1,
                        dtype=trig_stat.dtype
                    )
                while len(self.trig_stat_memory) < len(i1):
                    self.trig_stat_memory = numpy.resize(
                        self.trig_stat_memory,
                        len(self.trig_stat_memory)*2
                    )
                self.trig_stat_memory[:len(i1)] = trig_stat

                # Force data into form needed by stat.py and then compute the
                # ranking statistic values.
                sngls_list = [[fixed_ifo, self.trig_stat_memory[:len(i1)]],
                              [shift_ifo, stats[i1]]]
                c = self.stat_calculator.rank_stat_coinc(
                    sngls_list,
                    slide,
                    self.timeslide_interval,
                    shift_vec
                )

                # Store data about new triggers: slide index, stat value and
                # times.
                offsets.append(slide)
                cstat.append(c)
                ctimes[shift_ifo].append(times[i1])
                ctimes[fixed_ifo].append(numpy.zeros(len(c),
                                         dtype=numpy.float64))
                ctimes[fixed_ifo][-1].fill(trig_time)

                # As background triggers are removed after a certain time, we
                # need to log when this will be for new background triggers.
                single_expire[shift_ifo].append(
                    self.singles[shift_ifo].expire_vector(template)[i1]
                )
                single_expire[fixed_ifo].append(numpy.zeros(len(c),
                                                dtype=numpy.int32))
                single_expire[fixed_ifo][-1].fill(
                    self.singles[fixed_ifo].time - 1
                )

                # Save the template and trigger ids to keep association
                # to singles. The trigger was just added so it must be in
                # the last position: we mark this with -1 so the
                # slicing picks the right point
                template_ids.append(numpy.zeros(len(c)) + template)
                trigger_ids[shift_ifo].append(i1)
                trigger_ids[fixed_ifo].append(numpy.zeros(len(c)) - 1)

        cstat = numpy.concatenate(cstat)
        template_ids = numpy.concatenate(template_ids).astype(numpy.int32)
        for ifo in valid_ifos:
            trigger_ids[ifo] = numpy.concatenate(trigger_ids[ifo]).astype(numpy.int32)

        logging.info(
            "%s: %s background and zerolag coincs",
            ppdets(self.ifos, "-"), len(cstat)
        )

        # Cluster the triggers we've found
        # (both zerolag and shifted are handled together)
        num_zerolag = 0
        num_background = 0
        if len(cstat) > 0:
            offsets = numpy.concatenate(offsets)
            ctime0 = numpy.concatenate(ctimes[self.ifos[0]]).astype(numpy.float64)
            ctime1 = numpy.concatenate(ctimes[self.ifos[1]]).astype(numpy.float64)
            logging.info("Clustering %s coincs", ppdets(self.ifos, "-"))
            cidx = cluster_coincs(cstat, ctime0, ctime1, offsets,
                                  self.timeslide_interval,
                                  self.analysis_block + 2*self.time_window,
                                  method='cython')
            offsets = offsets[cidx]
            zerolag_idx = (offsets == 0)
            bkg_idx = (offsets != 0)

            for ifo in self.ifos:
                single_expire[ifo] = numpy.concatenate(single_expire[ifo])
                single_expire[ifo] = single_expire[ifo][cidx][bkg_idx]

            self.coincs.add(cstat[cidx][bkg_idx], single_expire, valid_ifos)
            num_zerolag = zerolag_idx.sum()
            num_background = bkg_idx.sum()
        elif len(valid_ifos) > 0:
            self.coincs.increment(valid_ifos)

        # Collect coinc results for saving
        coinc_results = {}
        # Save information about zerolag triggers
        if num_zerolag > 0:
            idx = cidx[zerolag_idx][0]
            zerolag_cstat = cstat[cidx][zerolag_idx]
            ifar, ifar_sat = self.ifar(zerolag_cstat)
            zerolag_results = {
                'foreground/ifar': ifar,
                'foreground/ifar_saturated': ifar_sat,
                'foreground/stat': zerolag_cstat,
                'foreground/type': '-'.join(self.ifos)
            }
            template = template_ids[idx]
            for ifo in self.ifos:
                trig_id = trigger_ids[ifo][idx]
                single_data = self.singles[ifo].data(template)[trig_id]
                for key in single_data.dtype.names:
                    path = f'foreground/{ifo}/{key}'
                    zerolag_results[path] = single_data[key]
            coinc_results.update(zerolag_results)

        # Save some summary statistics about the background
        coinc_results['background/time'] = numpy.array([self.background_time])
        coinc_results['background/count'] = len(self.coincs.data)

        # Save all the background triggers
        if self.return_background:
            coinc_results['background/stat'] = self.coincs.data

        return num_background, coinc_results

    def backout_last(self, updated_singles, num_coincs):
        """Remove the recently added singles and coincs

        Parameters
        ----------
        updated_singles: dict of numpy.ndarrays
            Array of indices that have been just updated in the internal
            buffers of single detector triggers.
        num_coincs: int
            The number of coincs that were just added to the internal buffer
            of coincident triggers
        """
        for ifo in updated_singles:
            self.singles[ifo].discard_last(updated_singles[ifo])
        self.coincs.remove(num_coincs)

    def add_singles(self, results):
        """Add singles to the background estimate and find candidates

        Parameters
        ----------
        results: dict
            Dictionary of dictionaries indexed by ifo and keys such as 'snr',
            'chisq', etc. The specific format is determined by the
            LiveBatchMatchedFilter class.

        Returns
        -------
        coinc_results: dict of arrays
            A dictionary of arrays containing the coincident results.
        """
        # Let's see how large everything is
        logging.info(
            "%s: %s coincs, %s bytes",
            ppdets(self.ifos, "-"), len(self.coincs), self.coincs.nbytes
        )

        # If there are no results just return
        valid_ifos = [k for k in results.keys() if results[k] and k in self.ifos]
        if len(valid_ifos) == 0: return {}

        # Add single triggers to the internal buffer
        self._add_singles_to_buffer(results, ifos=valid_ifos)

        # Calculate zerolag and background coincidences
        _, coinc_results = self._find_coincs(results, valid_ifos=valid_ifos)

        # record if a coinc is possible in this chunk
        if len(valid_ifos) == 2:
            coinc_results['coinc_possible'] = True

        return coinc_results


__all__ = [
    "background_bin_from_string",
    "timeslide_durations",
    "time_coincidence",
    "time_multi_coincidence",
    "cluster_coincs",
    "cluster_coincs_multiifo",
    "mean_if_greater_than_zero",
    "cluster_over_time",
    "MultiRingBuffer",
    "CoincExpireBuffer",
    "LiveCoincTimeslideBackgroundEstimator"
]
