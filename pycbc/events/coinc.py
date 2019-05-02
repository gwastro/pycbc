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

import h5py, numpy, logging, pycbc.pnutils, copy, lal
from pycbc.detector import Detector


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

        if boundary[0:2] == 'lt':
            member_func = lambda vals, bd=boundary : vals < float(bd[2:])
        elif boundary[0:2] == 'gt':
            member_func = lambda vals, bd=boundary : vals > float(bd[2:])
        else:
            raise RuntimeError("Can't parse boundary condition! Must begin "
                               "with 'lt' or 'gt'")

        if bin_type == 'component' and boundary[0:2] == 'lt':
            # maximum component mass is less than boundary value
            vals = numpy.maximum(data['mass1'], data['mass2'])
        if bin_type == 'component' and boundary[0:2] == 'gt':
            # minimum component mass is greater than bdary
            vals = numpy.minimum(data['mass1'], data['mass2'])
        elif bin_type == 'total':
            vals = data['mass1'] + data['mass2']
        elif bin_type == 'chirp':
            vals = pycbc.pnutils.mass1_mass2_to_mchirp_eta(
                                               data['mass1'], data['mass2'])[0]
        elif bin_type == 'SEOBNRv2Peak':
            vals = pycbc.pnutils.get_freq('fSEOBNRv2Peak',
                  data['mass1'], data['mass2'], data['spin1z'], data['spin2z'])
        elif bin_type == 'SEOBNRv4Peak':
            vals = pycbc.pnutils.get_freq('fSEOBNRv4Peak', data['mass1'],
                                          data['mass2'], data['spin1z'],
                                          data['spin2z'])
        elif bin_type == 'SEOBNRv2duration':
            vals = pycbc.pnutils.get_imr_duration(data['mass1'], data['mass2'],
                               data['spin1z'], data['spin2z'], data['f_lower'],
                                                        approximant='SEOBNRv2')
        else:
            raise ValueError('Invalid bin type %s' % bin_type)

        locs = member_func(vals)
        del vals

        # make sure we don't reuse anything from an earlier bin
        locs = numpy.where(locs)[0]
        locs = numpy.delete(locs, numpy.where(numpy.in1d(locs, used))[0])
        used = numpy.concatenate([used, locs])
        bins[name] = locs

    return bins

def calculate_n_louder(bstat, fstat, dec, skip_background=False):
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
    skip_background: optional, {boolean, False}
        Skip calculating cumulative numbers for background triggers

    Returns
    -------
    cum_back_num: numpy.ndarray
        The cumulative array of background triggers. Does not return this
        argument if skip_background == True
    fore_n_louder: numpy.ndarray
        The number of background triggers above each foreground trigger
    """
    sort = bstat.argsort()
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
    idx = numpy.searchsorted(bstat, fstat, side='left') - 1

    # If the foreground are *quieter* than the background or at the same value
    # then the search sorted alorithm will choose position -1, which does not exist
    # We force it back to zero.
    if isinstance(idx, numpy.ndarray): # Handle the case where our input is an array
        idx[idx < 0] = 0
    else: # Handle the case where we are simply given a scalar value
        if idx < 0:
            idx = 0

    fore_n_louder = n_louder[idx]

    if not skip_background:
        unsort = sort.argsort()
        back_cum_num = n_louder[unsort]
        return back_cum_num, fore_n_louder
    else:
        return fore_n_louder

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
        idx2 = numpy.array([], dtype=numpy.int64)

    if slide_step:
        diff = ((t1 / slide_step)[idx1] - (t2 / slide_step)[idx2])
        slide = numpy.rint(diff)
    else:
        slide = numpy.zeros(len(idx1))

    return idx1.astype(numpy.uint32), idx2.astype(numpy.uint32), slide.astype(numpy.int32)


def time_multi_coincidence(times, slide_step=0, slop=.003,
                           pivot='H1', fixed='L1'):
    """ Find multi detector concidences.

    Parameters
    ----------
    times: dict of numpy.ndarrays
        Dictionary keyed by ifo of the times of each single detector trigger.
    slide_step: float
        The interval between time slides
    slop: float
        The amount of time to add to the TOF between detectors for coincidence
    pivot: str
        ifo used to test coincidence against in first stage
    fixed: str
        the other ifo used in the first stage coincidence which we'll use
        as a fixed time reference for coincident triggers. All other detectors
        are time slid by being fixed to this detector.
    """
    # pivots are used to determine standard coincidence triggers, we then
    # pair off additional detectors to those.
    def win(ifo1, ifo2):
        d1 = Detector(ifo1)
        d2 = Detector(ifo2)
        return d1.light_travel_time_to_detector(d2) + slop

    # Find coincs first between the two fully time-slid detectors
    pivot_id, fix_id, slide = time_coincidence(times[pivot], times[fixed],
                                               win(pivot, fixed),
                                               slide_step=slide_step)

    # additional detectors do not slide independently of the fixed one
    # Each trigger in an additional detector must be concident with an
    # existing coincident one. All times moved to 'fixed' relative time
    fixed_time = times[fixed][fix_id]
    pivot_time = times[pivot][pivot_id] - slide_step * slide

    ctimes = {fixed: fixed_time, pivot:pivot_time}
    ids = {fixed:fix_id, pivot:pivot_id}

    dep_ifos = [ifo for ifo in times.keys() if ifo != fixed and ifo != pivot]
    for ifo1 in dep_ifos:
        otime = times[ifo1]
        sort = times[ifo1].argsort()
        time = otime[sort]

        # Find coincidences between dependent ifo triggers and existing coinc.
        for ifo2 in ids.keys():
            # Currently assumes that additional detectors do not slide
            # independently of the 'fixed one'
            #
            # To modify that assumption, the code here would be modified
            # by adding a function that remaps the coinc time frame and unmaps
            # it and the end of this loop.
            # This remapping must ensure
            #    * function of the standard slide number
            #    * ensure all times remain within coincident segment
            #    * unbiased distribution of triggers after mapping.

            w = win(ifo1, ifo2)
            left = numpy.searchsorted(time, ctimes[ifo2] - w)
            right = numpy.searchsorted(time, ctimes[ifo2] + w)

            # remove elements that will not form a coinc
            # There is only at most one trigger for an existing coinc
            # (assumes triggers spaced > slide step)
            nz = (right - left).nonzero()
            dep_ids = left[nz]

            # The property that only one trigger can be within the window is ensured
            # by the peak finding algorithm we use for each template.
            # If that is modifed, this function may need to be
            # extended.
            if len(left) > 0 and (right - left).max() > 1:
                raise ValueError('Somehow triggers are closer than time-delay window')

            slide = slide[nz]
            for ifo in ctimes:
                ctimes[ifo] = ctimes[ifo][nz]
                ids[ifo] = ids[ifo][nz]

        # Add this detector now to the cumulative set and proceed to the next
        # ifo coincidence test
        ids[ifo1] = sort[dep_ids]
        ctimes[ifo1] = otime[ids[ifo1]]

    return ids, slide

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

    if numpy.isfinite(slide):
        # for a time shifted coinc, time1 is greater than time2 by approximately timeslide_id*slide
        # adding this quantity gives a mean coinc time located around time1
        time = (time1 + time2 + timeslide_id * slide) / 2
    else:
        time = 0.5 * (time2 + time1)

    tslide = timeslide_id.astype(numpy.float128)
    time = time.astype(numpy.float128)

    span = (time.max() - time.min()) + window * 10
    time = time + span * tslide
    cidx = cluster_over_time(stat, time, window, argmax)
    return cidx

def cluster_coincs_multiifo(stat, time_coincs, timeslide_id, slide, window, argmax=numpy.argmax):
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
    time_coinc_zip = zip(*time_coincs)
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

    tslide = timeslide_id.astype(numpy.float128)
    time_avg = time_avg.astype(numpy.float128)

    span = (time_avg.max() - time_avg.min()) + window * 10
    time_avg = time_avg + span * tslide
    cidx = cluster_over_time(stat, time_avg, window, argmax)

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

def cluster_over_time(stat, time, window, argmax=numpy.argmax):
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
    argmax: function
        the function used to calculate the maximum value

    Returns
    -------
    cindex: numpy.ndarray
        The set of indices corresponding to the surviving coincidences.
    """
    logging.info('Clustering events over %s s window', window)

    indices = []
    time_sorting = time.argsort()
    stat = stat[time_sorting]
    time = time[time_sorting]

    left = numpy.searchsorted(time, time - window)
    right = numpy.searchsorted(time, time + window)
    indices = numpy.zeros(len(left), dtype=numpy.uint32)

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

    logging.info('%d triggers remaining', len(indices))
    return time_sorting[indices]


class MultiRingBuffer(object):
    """Dynamic size n-dimensional ring buffer that can expire elements."""

    def __init__(self, num_rings, max_time, dtype):
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
        """
        self.max_time = max_time
        self.buffer = []
        self.buffer_expire = []
        for _ in range(num_rings):
            self.buffer.append(numpy.zeros(0, dtype=dtype))
            self.buffer_expire.append(numpy.zeros(0, dtype=int))
        self.time = 0

    @property
    def filled_time(self):
        return min(self.time, self.max_time)

    def num_elements(self):
        return sum([len(a) for a in self.buffer])

    @property
    def nbytes(self):
        return sum([a.nbytes for a in self.buffer])

    def discard_last(self, indices):
        """Discard the triggers added in the latest update"""
        for i in indices:
            self.buffer_expire[i] = self.buffer_expire[i][:-1]
            self.buffer[i] = self.buffer[i][:-1]

    def advance_time(self):
        """Advance the internal time increment by 1, expiring any triggers that
        are now too old.
        """
        self.time += 1

    def add(self, indices, values):
        """Add triggers in 'values' to the buffers indicated by the indices
        """
        for i, v in zip(indices, values):
            self.buffer[i] = numpy.append(self.buffer[i], v)
            self.buffer_expire[i] = numpy.append(self.buffer_expire[i], self.time)
        self.advance_time()

    def expire_vector(self, buffer_index):
        """Return the expiration vector of a given ring buffer """
        return self.buffer_expire[buffer_index]

    def data(self, buffer_index):
        """Return the data vector for a given ring buffer"""
        # Check for expired elements and discard if they exist
        expired = self.time - self.max_time     
        exp = self.buffer_expire[buffer_index]
        j = 0
        while j < len(exp):
            # Everything before this j must be expired
            if exp[j] >= expired:
                self.buffer_expire[buffer_index] = exp[j:].copy()
                self.buffer[buffer_index] = self.buffer[buffer_index][j:].copy()
                break
            j += 1
       
        return self.buffer[buffer_index]

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
        return self.buffer.nbytes

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
            self.buffer.resize(newlen)

        self.buffer[self.index:self.index+len(values)] = values
        if len(values) > 0:
            for ifo in self.ifos:
                self.timer[ifo][self.index:self.index+len(values)] = times[ifo]

            self.index += len(values)

        # Remove the expired old elements
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
        return (self.buffer[:self.index] > value).sum()

    @property
    def data(self):
        """Return the array of elements"""
        return self.buffer[:self.index]

class LiveCoincTimeslideBackgroundEstimator(object):
    """Rolling buffer background estimation."""

    def __init__(self, num_templates, analysis_block, background_statistic,
                 stat_files, ifos,
                 ifar_limit=100,
                 timeslide_interval=.035,
                 coinc_threshold=.002,
                 return_background=False):
        """
        Parameters
        ----------
        num_templates: int
            The size of the template bank
        analysis_block: int
            The number of seconds in each analysis segment
        background_statistic: str
            The name of the statistic to rank coincident events.
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
        """
        from . import stat
        self.num_templates = num_templates
        self.analysis_block = analysis_block

        # Only pass a valid stat file for this ifo pair
        for fname in stat_files:
            f = h5py.File(fname, 'r')
            ifos_set = set([f.attrs['ifo0'], f.attrs['ifo1']])
            f.close()
            if ifos_set == set(ifos):
                stat_files = [fname]
                logging.info('Setup ifos %s-%s with file %s and stat %s',
                             ifos[0], ifos[1], fname, background_statistic)

        self.stat_calculator = stat.get_statistic(background_statistic)(stat_files)

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
        return cls(num_templates, analysis_chunk,
                   args.background_statistic,
                   args.background_statistic_files,
                   return_background=args.store_background,
                   ifar_limit=args.background_ifar_limit,
                   timeslide_interval=args.timeslide_interval,
                   ifos=ifos)

    @staticmethod
    def insert_args(parser):
        from . import stat

        group = parser.add_argument_group('Coincident Background Estimation')
        group.add_argument('--background-statistic', default='newsnr',
            choices=sorted(stat.statistic_dict.keys()),
            help="Ranking statistic to use for candidate coincident events")
        group.add_argument('--background-statistic-files', nargs='+',
            help="Files containing precalculate values to calculate ranking"
                 " statistic values", default=[])
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
        from six.moves import cPickle
        cPickle.dump(self, filename)

    @staticmethod
    def restore_state(filename):
        """Restore state of the background buffers from a file"""
        from six.moves import cPickle
        return cPickle.load(filename)

    def ifar(self, coinc_stat):
        """Return the far that would be associated with the coincident given.
        """
        n = self.coincs.num_greater(coinc_stat)
        return self.background_time / lal.YRJUL_SI / (n + 1)

    def set_singles_buffer(self, results):
        """Create the singles buffer

        This creates the singles buffer for each ifo. The dtype is determined
        by a representative sample of the single triggers in the results.

        Parameters
        ----------
        restuls: dict of dict
            Dict indexed by ifo and then trigger column.
        """
        # Determine the dtype from a sample of the data.
        self.singles_dtype = []
        data = False
        for ifo in self.ifos:
            if ifo in results and results[ifo] is not False:
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
        results: dict of arrays
            Dictionary of dictionaries indexed by ifo and keys such as 'snr',
            'chisq', etc. The specific format it determined by the
            LiveBatchMatchedFilter class.

        Returns
        -------
        updated_singles: dict of numpy.ndarrays
            Array of indices that have been just updated in the internal
            buffers of single detector triggers.
        """
        if len(self.singles.keys()) == 0:
            self.set_singles_buffer(results)

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

    def _find_coincs(self, results, ifos):
        """Look for coincs within the set of single triggers

        Parameters
        ----------
        results: dict of arrays
            Dictionary of dictionaries indexed by ifo and keys such as 'snr',
            'chisq', etc. The specific format it determined by the
            LiveBatchMatchedFilter class.

        Returns
        -------
        coinc_results: dict of arrays
            A dictionary of arrays containing the coincident results.
        """
        # for each single detector trigger find the allowed coincidences
        # Record which template and the index of the single trigger
        # that forms each coincident trigger
        cstat = [[]]
        offsets = []
        ctimes = {self.ifos[0]:[], self.ifos[1]:[]}
        single_expire = {self.ifos[0]:[], self.ifos[1]:[]}
        template_ids = [[]]
        trigger_ids = {self.ifos[0]:[[]], self.ifos[1]:[[]]}

        # Calculate all the permutations of coincident triggers for each
        # new single detector trigger collected
        for ifo in ifos:
            trigs = results[ifo]

            oifo = self.ifos[1] if self.ifos[0] == ifo else self.ifos[0]

            for i in range(len(trigs['end_time'])):
                trig_stat = trigs['stat'][i]
                trig_time = trigs['end_time'][i]
                template = trigs['template_id'][i]

                times = self.singles[oifo].data(template)['end_time']
                stats = self.singles[oifo].data(template)['stat']

                i1, _, slide = time_coincidence(times,
                                 numpy.array(trig_time, ndmin=1,
                                 dtype=numpy.float64),
                                 self.time_window,
                                 self.timeslide_interval)
                trig_stat = numpy.resize(trig_stat, len(i1))
                c = self.stat_calculator.coinc(stats[i1], trig_stat,
                                               slide, self.timeslide_interval)
                offsets.append(slide)
                cstat.append(c)
                ctimes[oifo].append(times[i1])
                ctimes[ifo].append(numpy.zeros(len(c), dtype=numpy.float64))
                ctimes[ifo][-1].fill(trig_time)

                single_expire[oifo].append(self.singles[oifo].expire_vector(template)[i1])
                single_expire[ifo].append(numpy.zeros(len(c),
                                          dtype=numpy.int32))
                single_expire[ifo][-1].fill(self.singles[ifo].time - 1)

                # save the template and trigger ids to keep association
                # to singles. The trigger was just added so it must be in
                # the last position we mark this with -1 so the
                # slicing picks the right point
                template_ids.append(numpy.zeros(len(c)) + template)
                trigger_ids[oifo].append(i1)
                trigger_ids[ifo].append(numpy.zeros(len(c)) - 1)

        cstat = numpy.concatenate(cstat)
        template_ids = numpy.concatenate(template_ids).astype(numpy.int32)
        for ifo in ifos:
            trigger_ids[ifo] = numpy.concatenate(trigger_ids[ifo]).astype(numpy.int32)

        # cluster the triggers we've found
        # (both zerolag and non handled together)
        num_zerolag = 0
        num_background = 0

        logging.info('%s background and zerolag coincs', len(cstat))
        if len(cstat) > 0:
            offsets = numpy.concatenate(offsets)
            ctime0 = numpy.concatenate(ctimes[self.ifos[0]]).astype(numpy.float64)
            ctime1 = numpy.concatenate(ctimes[self.ifos[1]]).astype(numpy.float64)

            cidx = cluster_coincs(cstat, ctime0, ctime1, offsets,
                                      self.timeslide_interval,
                                      self.analysis_block)
            offsets = offsets[cidx]
            zerolag_idx = (offsets == 0)
            bkg_idx = (offsets != 0)

            for ifo in self.ifos:
                single_expire[ifo] = numpy.concatenate(single_expire[ifo])
                single_expire[ifo] = single_expire[ifo][cidx][bkg_idx]

            self.coincs.add(cstat[cidx][bkg_idx], single_expire, ifos)
            num_zerolag = zerolag_idx.sum()
            num_background = bkg_idx.sum()
        elif len(ifos) > 0:
            self.coincs.increment(ifos)

        ####################################Collect coinc results for saving
        coinc_results = {}
        # Save information about zerolag triggers
        if num_zerolag > 0:
            zerolag_results = {}
            idx = cidx[zerolag_idx][0]
            zerolag_cstat = cstat[cidx][zerolag_idx]
            zerolag_results['foreground/ifar'] = self.ifar(zerolag_cstat)
            zerolag_results['foreground/stat'] = zerolag_cstat
            template = template_ids[idx]
            for ifo in self.ifos:
                trig_id = trigger_ids[ifo][idx]
                single_data = self.singles[ifo].data(template)[trig_id]
                for key in single_data.dtype.names:
                    path = 'foreground/%s/%s' % (ifo, key)
                    zerolag_results[path] = single_data[key]

            zerolag_results['foreground/type'] = '-'.join(self.ifos)

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
        """Add singles to the bacckground estimate and find candidates

        Parameters
        ----------
        results: dict of arrays
            Dictionary of dictionaries indexed by ifo and keys such as 'snr',
            'chisq', etc. The specific format it determined by the
            LiveBatchMatchedFilter class.

        Returns
        -------
        coinc_results: dict of arrays
            A dictionary of arrays containing the coincident results.
        """
        # Let's see how large everything is
        logging.info('BKG Coincs %s stored %s bytes',
                     len(self.coincs), self.coincs.nbytes)

        # If there are no results just return
        valid_ifos = [k for k in results.keys() if results[k] and k in self.ifos]
        if len(valid_ifos) == 0: return {}

        # Add single triggers to the internal buffer
        self._add_singles_to_buffer(results, ifos=valid_ifos)

        # Calculate zerolag and background coincidences
        _, coinc_results = self._find_coincs(results, ifos=valid_ifos)

        # record if a coinc is possible in this chunk
        if len(valid_ifos) == 2:
            coinc_results['coinc_possible'] = True

        return coinc_results
