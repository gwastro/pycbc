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
This modules provides classes and functions for tapering waveforms.
"""

import numpy
from scipy import signal
from pycbc import pnutils
from pycbc.types import TimeSeries, FrequencySeries
from pycbc.waveform.utils import taper_timeseries, time_from_frequencyseries

class TDomainTaper(object):
    def __init__(self, left_taper=None, right_taper=None,
                 left_taper_duration=None, right_taper_duration=None,
                 left_taper_time=None, right_taper_time=None,
                 left_taper_frequency=None, right_taper_frequency=None,
                 left_taper_freqfunc=None, right_taper_freqfunc=None,
                 taper_whitened=False, psds=None):
        if taper_whitened and not (taper_whitened == 1 or taper_whitened == 2):
            raise ValueError("taper_whitened must be either False (taper "
                             "before whitening), 1 (taper after whitening) "
                             "or 2 (taper after overwhitening)")
        if left_taper == 'lal':
            if left_taper_duration is not None:
                raise ValueError("The lal taper function does not take a "
                                 "duration")
            if left_taper_time is not None or \
                    left_taper_frequency is not None or \
                    left_taper_freqfunc is not None:
                raise ValueError("The lal taper function does not take a "
                                 "start time or frequency")
        elif left_taper is not None and left_taper_duration is None:
            raise ValueError("Non-lal taper functions require a duration")
        elif left_taper is not None and (left_taper_time is None and
                                         left_taper_frequency is None and
                                         left_taper_freqfunc is None):
            raise ValueError("Non-lal taper functions require either a taper "
                             "time, taper frequency, or frequency function")
        self.left_taper = left_taper
        self.left_taper_duration = left_taper_duration
        self.left_taper_time = left_taper_time
        self.left_taper_frequency = left_taper_frequency
        self.left_taper_freqfunc = left_taper_freqfunc
        if right_taper == 'lal':
            if right_taper_duration is not None:
                raise ValueError("The lal taper function does not take a "
                                 "duration")
            if left_taper_time is not None or \
                    left_taper_frequency is not None or \
                    left_taper_freqfunc is not None:
                raise ValueError("The lal taper function does not take a "
                                 "end time or frequency")
        elif right_taper is not None and right_taper_duration is None:
            raise ValueError("Non-lal taper functions require a duration")
        elif right_taper is not None and (right_taper_time is None and
                                          right_taper_frequency is None and
                                          right_taper_freqfunc is None):
            raise ValueError("Non-lal taper functions require either a taper "
                             "time, taper frequency, or frequency function")
        self.right_taper = right_taper
        self.right_taper_duration = right_taper_duration
        self.right_taper_time = right_taper_time
        self.right_taper_frequency = right_taper_frequency
        self.right_taper_freqfunc = right_taper_freqfunc
        self.taper_whitened = taper_whitened
        self.left_window = {}
        self.right_window = {}
        self.psds = psds
        self.asds = {}
        if self.taper_whitened:
            if psds is None:
                raise ValueError("must provide a psd if tapering "
                                "(over-)whitened waveform")
            if left_taper == 'lal' and right_taper == 'lal':
                raise ValueError("both left and right use lal tapering, but "
                                 "lal tapering cannot be done on whitened "
                                 "waveforms")
            self.whkmin = {}
            self.whkmax = {}
            for ifo,psd in psds.items():
                if self.taper_whitened == 1:
                    asd = psd**0.5
                    self.asds[ifo] = asd 
                    nzidx = numpy.nonzero(asd.data)[0]
                else:
                    nzidx = numpy.nonzero(psd.data)[0]
                self.whkmin[ifo] = nzidx[0]
                self.whkmax[ifo] = nzidx[-1] + 1

    def get_left_window(self, delta_t):
        taper_size = int(self.left_taper_duration / delta_t)
        try:
            return self.left_window[taper_size]
        except KeyError:
            # generate the window at this dt
            win = signal.get_window(self.left_taper, 2*taper_size)
            self.left_window[taper_size] = win[:taper_size]
            return self.left_window[taper_size]

    def get_right_window(self, delta_t):
        taper_size = int(self.right_taper_duration / delta_t)
        try:
            return self.right_window[taper_size]
        except KeyError:
            # generate the window at this dt
            win = signal.get_window(self.right_taper, 2*taper_size)
            self.right_window[taper_size] = win[taper_size:]
            return self.right_window[taper_size]

    def whiten_waveform(htilde, ifo):
        if self.taper_whitened == 1:
            htilde[self.whkmin:self.whkmax] /= \
                self.asds[ifo][self.whkmin:self.whkmax]
        else:
            htilde[self.whkmin:self.whkmax] /= \
                self.psds[ifo][self.whkmin:self.whkmax]
        htilde.data[:self.whkmin] = 0.
        htilde.data[self.whkmax:] = 0.

    def apply_taper(self, h, copy=True,
                    ifo=None, params=None, sample_frequencies=None):
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
        # lal taper function needs to be applied before whitening
        if self.left_taper == 'lal' or \
                self.right_taper == 'lal':
            if ht is None:
                ht = hf.to_timeseries()
            tmeth = ''
            if self.left_taper == 'lal':
                tmeth = 'start'
            if self.right_taper == 'lal':
                tmeth = ''.join([tmeth, 'end'])
            ht = taper_timeseries(ht, tapermethod=tmeth)
            hf = None
            if tmeth == 'startend':
                # just return, since there's nothing else to do
                if return_f:
                    return ht.to_frequencyseries()
                else:
                    return ht
        if self.taper_whitened:
            if ifo is None:
                raise ValueError("must provide an ifo to whiten with")
            if hf is None:
                hf = ht.to_frequencyseries(delta_f=self.psds[ifo].delta_f)
            self.whiten_waveform(hf, ifo)
            ht = hf.to_timeseries()
        elif ht is None:
            ht = hf.to_timeseries()
        # 
        #   left taper
        #
        left_time = self.left_taper_time
        left_freq = self.left_taper_frequency
        if self.left_taper_freqfunc is not None:
            if params is None:
                raise ValueError("must provide waveform parameters for the "
                                 "frequency function to use for the left")
            left_freq = pnutils.named_frequency_cutoffs[
                self.left_taper_freqfunc](params)
        t_of_f = None
        if left_freq is not None:
            # need frequencyseries version to get f(t)
            if hf is None:
                hf = ht.to_frequencyseries()
            t_of_f = time_from_frequencyseries(hf,
                sample_frequencies=sample_frequencies)
            t = t_of_f[int(left_freq / hf.delta_f)]
            if left_time is None:
                left_time = t
            else:
                left_time = max(t, left_time)
        # apply
        if left_time is not None:
            win = self.get_left_window(ht.delta_t)
            endidx = min(int(numpy.ceil(left_time / ht.delta_t)), len(ht))
            startidx = max(endidx - len(win), 0)
            ht.data[startidx:endidx] *= win[-(endidx-startidx):]
            ht.data[:startidx] = 0.
        #
        #   right taper
        #
        right_time = self.right_taper_time
        right_freq = self.right_taper_frequency
        if self.right_taper_freqfunc is not None:
            if params is None:
                raise ValueError("must provide waveform parameters for the "
                                 "frequency function to use for the right")
            right_freq = pnutils.named_frequency_cutoffs[
                self.right_taper_freqfunc](params)
        if right_freq is not None:
            # need frequencyseries version to get f(t)
            if hf is None:
                hf = ht.to_frequencyseries()
            if t_of_f is None:
                t_of_f = time_from_frequencyseries(hf,
                    sample_frequencies=sample_frequencies)
            t = t_of_f[int(numpy.ceil(right_freq / hf.delta_f))]
            if right_time is None:
                right_time = t
            else:
                right_time = min(t, right_time)
        # apply
        if right_time is not None:
            win = self.get_right_window(ht.delta_t)
            endidx = min(int(right_time / ht.delta_t), len(ht))
            startidx = max(endidx - len(win), 0)
            ht.data[startidx:endidx] *= win[:endidx-startidx]
            ht.data[endidx:] = 0.
        #
        #   Return
        #
        if return_f:
            return ht.to_frequencyseries()
        else:
            return ht

    __call__ = apply_taper

    @classmethod
    def from_config(cls, cp, section='taper', psds=None):
        opts = {}
        # parse the whitening
        if cp.has_option(section, 'taper-whitened'):
            taper_whitened = cp.get(section, 'taper-whitened')
            try:
                taper_whitened = int(taper_whitened)
            except ValueError:
                raise ValueError("taper-whitened must be either 0 (no "
                                 "whitening), 1 (whiten), or 2 (overwhiten)")
            opts['taper_whitened'] = taper_whitened
            opts['psds'] = psds
        # get everything else
        for opt in cp.options(section):
            if opt == 'taper-whitened':
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


__all__ = ['TDomainTaper']
