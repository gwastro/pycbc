# Copyright (C) 2017  Collin Capano
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
"""
This modules provides classes and functions for windowing data.
"""

import numpy
from scipy import signal
import lalsimulation as sim
from pycbc.types import Array, TimeSeries, FrequencySeries, float32, float64

# values for informing data whitening level
UNWHITENED = 0
WHITENED = 1
OVERWHITENED = 2

#
# =============================================================================
#
#                             LALSimulation window
#
# =============================================================================
#

# map between tapering string in sim_inspiral table or inspiral
# code option and lalsimulation constants
laltaper_map = {
    'TAPER_NONE'    : None,
    'TAPER_START'   : sim.SIM_INSPIRAL_TAPER_START,
    'start'         : sim.SIM_INSPIRAL_TAPER_START,
    'TAPER_END'     : sim.SIM_INSPIRAL_TAPER_END,
    'end'           : sim.SIM_INSPIRAL_TAPER_END,
    'TAPER_STARTEND': sim.SIM_INSPIRAL_TAPER_STARTEND,
    'startend'      : sim.SIM_INSPIRAL_TAPER_STARTEND
}

laltaper_func_map = {
    numpy.dtype(float32): sim.SimInspiralREAL4WaveTaper,
    numpy.dtype(float64): sim.SimInspiralREAL8WaveTaper
}

def laltaper_timeseries(tsdata, tapermethod=None, return_lal=False):
    """
    Taper either or both ends of a time series using wrapped
    LALSimulation functions

    Parameters
    ----------
    tsdata : TimeSeries
        Series to be tapered, dtype must be either float32 or float64
    tapermethod : string
        Should be one of ('TAPER_NONE', 'TAPER_START', 'TAPER_END',
        'TAPER_STARTEND', 'start', 'end', 'startend') - NB 'TAPER_NONE' will
        not change the series!
    return_lal : Boolean
        If True, return a wrapped LAL time series object, else return a
        PyCBC time series.
    """
    if tapermethod is None:
        raise ValueError("Must specify a tapering method (function was called"
                         "with tapermethod=None)")
    if tapermethod not in laltaper_map.keys():
        raise ValueError("Unknown tapering method %s, valid methods are %s" % \
                         (tapermethod, ", ".join(laltaper_map.keys())))
    if tsdata.dtype not in (float32, float64):
        raise TypeError("Strain dtype must be float32 or float64, not "
                    + str(tsdata.dtype))
    taper_func = laltaper_func_map[tsdata.dtype]
    # make a LAL TimeSeries to pass to the LALSim function
    ts_lal = tsdata.astype(tsdata.dtype).lal()
    if laltaper_map[tapermethod] is not None:
        taper_func(ts_lal.data, laltaper_map[tapermethod])
    if return_lal:
        return ts_lal
    else:
        return TimeSeries(ts_lal.data.data[:], delta_t=ts_lal.deltaT,
                          epoch=ts_lal.epoch)


#
# =============================================================================
#
#                        Generic time-domain window
#
# =============================================================================
#
class WindowBoundsError(Exception):
    """This should be raised if applying a window to a data set would result
    in all zeros being returned, e.g., if the end time of the right taper was
    before the start of the data.
    """
    pass


class TimeDomainWindow(object):
    """Provides generic windows for tapering data in the time domain.

    A left and right taper function may be specified to use different tapers,
    along with a taper duration for each. Tapers are applied such that the
    left taper ramps up to 1 and the right taper ramps down from 1, with zeros
    before/after.

    Available taper functions are `lal` or any window recognized by
    `scipy.signal.get_window`. If `lal`, the tapering function in
    LALSimulation is used. This taper does not use a time or a location.
    Instead, if the data are cyclic, it will apply the taper to the first (if
    it is the left taper; last if right) two cycles.  If the data are not
    cyclic, the taper will be applied to the first (last) half of the data.
    All other tapers require a taper duration to be specified on
    initialization, and a time to use it when being applied (see `apply_window`
    for details). For these tapers, the first half (for the left taper, second
    half for right) of the window generated by `scipy.signal.get_window` is
    used. No taper may also be specified on either side.

    Windows may be applied to unwhitened, whitened, or over-whitened data. If
    either of the later two, a psd or dictionary of psds must be provided on
    initialization to use for whitening.

    .. note::
        LAL tapering cannot be applied to whitened data, as the whitening will
        mess up the function's ability to find the first/last two cycles. For
        this reason, if lal tapering is requested on either side, whitening
        will be done after the lal taper is applied. If lal tapering is
        requested for both sides, then a ValueError will be raised if
        `window_whitened` is not False.

    Instances of this class may be called as a function; the `apply_window`
    method is called in that case.

    Parameters
    ----------
    left_taper : str or tuple, optional
        The name of the function to use for the left taper. May be either
        `lal`, or any name recognized by `scipy.signal.get_window`. Some taper
        functions require an additional parameter to be specified; if so,
        provide a tuple with the first argument the window name and the
        second the parameter. For details, see `scipy.signal.get_window`.
        Default is None, in which case no taper will be applied on the left.
    right_taper : str or tuple, optional
        Same as `left_taper`, but for the right side.
    left_taper_duration : float, optional
        The duration of the taper on the left side, in seconds. Required for
        scipy windows. A `ValueError` will be raised if `left_taper` is set
        to 'lal' and this is not None. Default is None.
    right_taper_duration : float, optional
        Same as `left_taper_duration`, but for the right side.
    window_whitened : {False/0, 1, 2}
        Optionally (over-)whiten the data before applying the window. If 1,
        data will be whitened; if 2, data will be over-whitened; if 0 or False,
        no whitening will be done. Default is False. If 1 or 2, psds must not
        be None.
    psds : (dict of) FrequencySeries, optional
        Required if `window_whitened` is 1 or 2. Either a FrequencySeries or a
        dictionary of FrequencySeries to use for whitening the data. If a
        dictionary, must be `detector name -> psd`. Default is None.

    Raises
    ------
    WindowBoundsError
        Raised when `apply_window` is called and the requested left (right)
        time of the window occurs after (before) the end (start) of the data,
        such that the entire data segment would be zeroed.

    Examples
    --------
    Create a window that uses a 5 second Hann window on the left and a 4
    second Kaiser window with `beta = 8` on the right:

    >>> from pycbc import window
    >>> windowts = window.TimeDomainWindow(left_taper='hann', right_taper=('kaiser', 8.), left_taper_duration=5., right_taper_duration=4.)

    Create a 16s time series of ones, and apply the above window to it,
    starting at 2s and ending at 15s:

    >>> from pycbc import types
    >>> ts = types.TimeSeries(numpy.ones(16), delta_t=1.)
    >>> windowts(ts, left_time=2., right_time=15.).data
    ArrayWithAligned([ 0.,  0.,  0.,  0.0954915,  0.3454915, 0.6545085,
                       0.9045085,  1.,  1.,  1., 1.,  1.,  0.78875245,
                       0.36897272,  0.08273982,  0. ])

    Over-whiten and window a CBC waveform using a lal taper on the left and
    a 10ms kaiser taper ending 20ms before coalescence on the right:

    >>> from pycbc import psd, waveform
    >>> hp, _ = waveform.get_td_waveform(approximant='SEOBNRv4', mass1=40., mass2=30., delta_t=1./2048, f_lower=15.)
    >>> aligopsd = psd.aLIGOZeroDetHighPower(4*2048/2+1, 1./4, 10.)
    >>> windowts = window.TimeDomainWindow(left_taper='lal', right_taper=('kaiser', 8), right_taper_duration=0.01, window_whitened=2, psds=aligopsd)
    >>> hp = windowts(hp, right_time=len(hp)*hp.delta_t-0.02)

    Apply the same taper to a frequency-domain waveform. Note that we have to
    supply a `break_time` indicating the start of the data series to account
    for the wrap around of the waveform:

    >>> hp, _ = waveform.get_fd_waveform(approximant='IMRPhenomD', mass1=40., mass2=30., delta_f=1./4, f_lower=15.)
    >>> break_time = 1.
    >>> hp = windowts(hp, right_time=1./hp.delta_f-break_time-0.02, break_time=break_time)

    See Also
    --------
    scipy.signal.get_window : function for generating windows
    laltaper_timeseries : function that only applies the lal taper
    """
    def __init__(self, left_taper=None, right_taper=None,
                 left_taper_duration=None, right_taper_duration=None,
                 window_whitened=False, psds=None):
        if int(window_whitened) not in [UNWHITENED, WHITENED, OVERWHITENED]:
            raise ValueError("window_whitened must be either {} (taper "
                             "before whitening), {} (taper after whitening) "
                             "or {} (taper after overwhitening)".format(
                                UNWHITENED, WHITENED, OVERWHITENED))
        if left_taper == 'lal':
            if left_taper_duration is not None:
                raise ValueError("The lal taper function does not take a "
                                 "duration")
        elif left_taper is not None and left_taper_duration is None:
            raise ValueError("Non-lal taper functions require a duration")
        self.left_taper = left_taper
        self.left_taper_duration = left_taper_duration
        if right_taper == 'lal':
            if right_taper_duration is not None:
                raise ValueError("The lal taper function does not take a "
                                 "duration")
        elif right_taper is not None and right_taper_duration is None:
            raise ValueError("Non-lal taper functions require a duration")
        self.right_taper = right_taper
        self.right_taper_duration = right_taper_duration
        self.window_whitened = window_whitened
        self.left_window = {}
        self.right_window = {}
        self._psds = None
        self._invpsds = {}
        self._invasds = {}
        self.set_psds(psds)
        if self.window_whitened:
            if psds is None:
                raise ValueError("must provide a psd if tapering "
                                "(over-)whitened waveform")
            if left_taper == 'lal' and right_taper == 'lal':
                raise ValueError("both left and right use lal tapering, but "
                                 "lal tapering cannot be done on whitened "
                                 "waveforms")

    def set_psds(self, psds):
        """Sets the psds attribute and calculates the inverse."""
        self._psds = psds
        self._invasds.clear()
        self._invpsds.clear()
        if psds is not None:
            # temporarily suppress numpy divide by 0 warning
            numpysettings = numpy.seterr(divide='ignore')
            if not isinstance(psds, dict):
                self._psds = {None: psds}
            for ifo,psd in self._psds.items():
                mask = psd.data == 0.
                invpsd = 1./psd
                invpsd.data[mask] = 0.
                self._invpsds[ifo] = invpsd
            numpy.seterr(**numpysettings)

    @property
    def psds(self):
        """Returns the psds attribute."""
        return self._psds

    def whiten_waveform(self, htilde, ifo=None, copy=False):
        """Whitens the given frequency domain data.

        If `window_whitened` = 1, the data will be divided by the ASD. If 2,
        the PSD. Otherwise, a ValueError is raised.

        Parameters
        ----------
        htilde : FrequencySeries
            The data in the frequency domain.
        ifo : str, optional
            If `psds` is a dictionary, the psd of the ifo to get.
        copy : bool, optional
            If True, the data will be copied before whitening. Otherwise, the
            data are whitened in place. Default is False.

        Returns
        -------
        FrequencySeries :
            The whitened data.
        """
        if copy:
            htilde = 1.*htilde
        if self.window_whitened == WHITENED:
            if self._invasds is None:
                self._invasds = {}
            try:
                wh = self._invasds[ifo]
            except KeyError:
                # compute the inverse asd
                wh = self._invpsds[ifo]**0.5
                self._invasds[ifo] = wh
        elif self.window_whitened == OVERWHITENED:
            wh = self._invpsds[ifo]
        else:
            raise ValueError("window_whitened set to {}".format(
                             self.window_whitened))
        kmax = len(htilde)
        # we can't whiten the waveform if it has frequencies outside of what
        # the psd has
        if kmax > len(wh):
            raise ValueError("htilde goes to higher frequency than the psd")
        htilde *= wh[:kmax]
        return htilde

    def get_left_window(self, delta_t):
        """Returns the left window to use for tapering.

        If the given `delta_t` has not previously been used, the window will
        be generated and cached to the `left_window` dict.

        Parameters
        ----------
        delta_t : float
            The dt of the time series the taper will be applied to.

        Returns
        -------
        Array :
            The window.
        """
        taper_size = int(self.left_taper_duration / delta_t)
        try:
            return self.left_window[taper_size]
        except KeyError:
            # generate the window at this dt
            win = signal.get_window(self.left_taper, 2*taper_size)
            self.left_window[taper_size] = Array(win[:taper_size])
            return self.left_window[taper_size]

    def get_right_window(self, delta_t):
        """Returns the right window to use for tapering.

        If the given `delta_t` has not previously been used, the window will
        be generated and cached to the `right_window` dict.

        Parameters
        ----------
        delta_t : float
            The dt of the time series the taper will be applied to.

        Returns
        -------
        Array :
            The window.
        """
        taper_size = int(self.right_taper_duration / delta_t)
        try:
            return self.right_window[taper_size]
        except KeyError:
            # generate the window at this dt
            win = signal.get_window(self.right_taper, 2*taper_size)
            self.right_window[taper_size] = Array(win[taper_size:])
            return self.right_window[taper_size]

    @staticmethod
    def _taper_by_idx(ht, win, startidx, endidx, side):
        """Determines how much of a window overlaps with the data and applies.
        """
        if endidx > 0 and startidx < len(ht):
            if startidx < 0:
                wstartidx = abs(startidx)
                startidx = 0
            else:
                wstartidx = 0
            if endidx > len(ht):
                wendidx = len(ht) - endidx
                endidx = len(ht)
            else:
                wendidx = len(win)
            ht[startidx:endidx] *= win[wstartidx:wendidx]
            # zero out everything before the taper
            if side == 'left' and startidx > 0:
                ht[:startidx] = 0.
            if side == 'right' and endidx < len(ht):
                ht[endidx:] = 0.

    def apply_window(self, h, left_time=None, right_time=None, break_time=0,
                     ifo=None, copy=True):
        """Applies the window at the given times.

        The provided left and right times should be the number of seconds from
        the start of the data segment to apply the tapers.  The input data,
        `h`, may be either a `TimeSeries` or `FrequencySeries`.  If the
        latter, the data will be transformed to the time domain before
        applying the tapers, then transformed back to a FrequencySeries before
        returning.

        Frequency-domain data (or time-domain data that was transformed from
        the frequency domain) treats time on loop, making it ambigous where in
        the array the start/end of the data occurs. For example, frequency
        domain waveforms will have the ringdown portion wrapped around to the
        beginning of the time series. To clear the ambiguity, a `break_time`
        option is provided to specify where in the time series the
        beginning/end is. The left and right times are measured from this
        point.

        If `window_whitened` is 1 or 2, the data will be whitened accordingly
        before applying the window. If `h` is a FrequencySeries, this means
        `h` will be whitened before transforming to the time domain. If `h`
        is a TimeSeries, the data will be transformed to the frequency domain,
        whitened, then transformed back to the time domain.

        If `left_taper` and `right_taper` are both None, this just returns
        (a copy of) the data.

        Parameters
        ----------
        h : TimeSeries or FrequencySeries
            The data to apply the window to.
        left_time : float
            The time at which to start the left taper, in number of seconds
            from the start of the data segment. This must be provided if
            `left_taper` is not None or 'lal'; otherwise, this must be None
            (the default). If the time is negative (i.e., before the start of
            the data), only the amount of the taper that overlaps the data
            will be applied. If the entire taper occurs before the start of
            the data, no left taper will be applied.  If the entire left taper
            occurs after the end of the data, a WindowBoundsError is raised.
        right_time : float
            The time at which to end the right taper, in number of seconds
            from the start of the data segment. This must be provided if
            `right_taper` is not None or 'lal'; otherwise, this must be None
            (the default). If the time is after the end of the data (if h is a
            `FrequencySeries`, this means `1/h.delta_f`; otherwise
            `len(h)*h.delta_t`), only the amount of the taper that overlaps
            the data will be applied. If the entire taper occurs after the end
            of the data, no right taper will be applied. If `right_time < 0`
            (the entire right taper occurs before the start of the data), a
            WindowBoundsError is raised.
        break_time : float, optional
            The number of seconds into the data segment to consider the start
            of the data. Default is 0.
        ifo : None, optional
            Must be provided if `window_whitened` is 1 or 2 and `psds` is a
            dictionary of psds.
        copy : bool, optional
            Whether to copy the data before applying the window/whitening. If
            False, the taper will be applied in place. Default is True.

        Returns
        -------
        TimeSeries or FrequencySeries
            The windowed data. If a FrequencySeries was provided, a
            FrequencySeries will be returned. Otherwise, a TimeSeries. If
            `window_whitened` is 1 or 2, the returned data will be
            (over-)whitened.
        """
        # check times
        if left_time is not None:
            # check that a time wasn't specified for lal tapers
            if self.left_taper == 'lal':
                raise ValueError("lal tapers do not require a time")
            # check that we won't be left with all zeros
            if isinstance(h, FrequencySeries):
                end_time = 1./h.delta_f
            else:
                end_time = len(h)*h.delta_t
            if left_time >= end_time:
                raise WindowBoundsError("start of left taper occurs after the "
                                        "end of the data")
        # ditto for right
        if right_time is not None:
            # check that a time wasn't specified for lal tapers
            if self.right_taper == 'lal':
                raise ValueError("lal tapers do not require a time")
            # check that we won't be left with all zeros
            if right_time < 0:
                raise WindowBoundsError("end of the right taper occurs before "
                                        "the start of the data")
        if copy:
            h = 1*h
        if isinstance(h, FrequencySeries):
            ht = None
            hf = h
            return_f = True
        else:
            ht = h
            hf = None
            return_f = False
        rollback = 0
        # lal taper function needs to be applied before whitening
        if self.left_taper == 'lal' or self.right_taper == 'lal':
            if ht is None:
                ht = hf.to_timeseries()
            if break_time:
                break_indx = -int(break_time/ht.delta_t)
                ht.roll(break_indx)
                rollback = -break_indx
                # set break time to 0, so don't roll again
                break_time = 0
            tmeth = ''
            if self.left_taper == 'lal':
                tmeth = 'start'
            if self.right_taper == 'lal':
                tmeth = ''.join([tmeth, 'end'])
            ht = laltaper_timeseries(ht, tapermethod=tmeth)
            hf = None
            if tmeth == 'startend':
                # just return, since there's nothing else to do
                if rollback:
                    ht.roll(rollback)
                if return_f:
                    return ht.to_frequencyseries()
                else:
                    return ht
        # check that a time is provided
        if left_time is None and not (self.left_taper is None or
                                      self.left_taper == 'lal') :
            raise ValueError("must provide a time for non-lal tapers")
        if right_time is None and not (self.right_taper is None or
                                      self.right_taper == 'lal') :
            raise ValueError("must provide a time for non-lal tapers")
        #
        #   Whiten
        #
        if self.window_whitened:
            if hf is None:
                hf = ht.to_frequencyseries(delta_f=self.psds[ifo].delta_f)
            self.whiten_waveform(hf, ifo=ifo, copy=False)
            ht = hf.to_timeseries()
        elif ht is None:
            ht = hf.to_timeseries()
        # roll if break time is not at the beginning
        if break_time:
            break_indx = -int(break_time/ht.delta_t)
            ht.roll(break_indx)
            rollback = -break_indx
        #
        #   apply left taper
        #
        if left_time is not None:
            if right_time is not None and right_time <= left_time:
                raise ValueError("right_time must be > left_time")
            win = self.get_left_window(ht.delta_t)
            startidx = int(left_time / ht.delta_t)
            endidx = startidx + len(win)
            # we don't have to worry about the startidx being > len(ht), since
            # that would have triggered the catch that the left time be before
            # the end of the data, above
            self._taper_by_idx(ht, win, startidx, endidx, 'left')
        #
        #   apply right taper
        #
        if right_time is not None:
            win = self.get_right_window(ht.delta_t)
            endidx = int(numpy.ceil(right_time / ht.delta_t))
            startidx = endidx - len(win)
            # we don't have to worry about the endidx being < 0, since
            # that would have triggered the catch that the left time be before
            # the end of the data, above
            self._taper_by_idx(ht, win, startidx, endidx, 'right')
        #
        #   Return
        #
        if rollback:
            ht.roll(rollback)
        if return_f:
            return ht.to_frequencyseries()
        else:
            return ht

    __call__ = apply_window

    @classmethod
    def from_config(cls, cp, section='window', psds=None):
        """Initializes and instance of this class from the given config file.

        If a taper requires additional parameters to be specified, they
        should be specifed as `(left|right)-taper-param`. All other arguments
        in the given section will be passed to the class as keyword arguments,
        with dashes replaced with underscores. For example, the following will
        return a window that applies a kaiser taper on the left with `beta=8`,
        and no taper on the right:

        .. code-block:: ini

            [{section}]
            left-taper = kaiser
            left-taper-param = 8
            left-taper-duration = 4

        Parameters
        ----------
        cp : ConfigParser
            Config file parser to retrieve the settings from.
        section : str, optional
            The section to retrieve the results from. Default is 'taper'.
        psds : (dict of) FrequencySeries, optional
            The psd(s) to use for whitening. Must be provided if
            `window-whitened` is in the config file and is set to 1 or 2.
        
        Returns
        -------
        TimeDomainWindow :
            Instance of this class initialized with the parameters specified
            in the config file.
        """
        opts = {}
        # parse the whitening
        if cp.has_option(section, 'window-whitened'):
            window_whitened = cp.get(section, 'window-whitened')
            try:
                window_whitened = int(window_whitened)
            except ValueError:
                raise ValueError("window-whitened must be either 0 (no "
                                 "whitening), 1 (whiten), or 2 (overwhiten)")
            opts['window_whitened'] = window_whitened
            opts['psds'] = psds
        # get everything else
        for opt in cp.options(section):
            if opt == 'window-whitened':
                continue
            val = cp.get(section, opt)
            try:
                val = float(val)
            except ValueError:
                pass
            opts[opt.replace('-', '_')] = val
        # if taper parameters were provided, add to the appropriate taper opt
        taper_param = opts.pop('left_taper_param', None)
        if taper_param is not None:
            try:
                opts['left_taper'] = (opts['left_taper'], taper_param)
            except KeyError:
                raise ValueError("left_taper_param provided, but no "
                                 "left_taper")
        taper_param = opts.pop('right_taper_param', None)
        if taper_param is not None:
            try:
                opts['right_taper'] = (opts['right_taper'], taper_param)
            except KeyError:
                raise ValueError("right_taper_param provided, but no "
                                 "right_taper")
        return cls(**opts)


__all__ = ['TimeDomainWindow', 'WindowBoundsError', 'laltaper_timeseries',
           'laltaper_map', 'laltaper_func_map']
