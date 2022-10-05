# Copyright (C) 2014  Tito Dal Canton, Josh Willis, Alex Nitz
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

"""
Provides a class representing a time series.
"""
import os as _os, h5py
from pycbc.types.array import Array, _convert, complex_same_precision_as, zeros
from pycbc.types.array import _nocomplex
from pycbc.types.frequencyseries import FrequencySeries
import lal as _lal
import numpy as _numpy
from scipy.io.wavfile import write as write_wav


class TimeSeries(Array):
    """Models a time series consisting of uniformly sampled scalar values.

    Parameters
    ----------
    initial_array : array-like
        Array containing sampled data.
    delta_t : float
        Time between consecutive samples in seconds.
    epoch : {None, lal.LIGOTimeGPS}, optional
        Time of the first sample in seconds.
    dtype : {None, data-type}, optional
        Sample data type.
    copy : boolean, optional
        If True, samples are copied to a new array.
    """

    def __init__(self, initial_array, delta_t=None,
                 epoch=None, dtype=None, copy=True):
        if len(initial_array) < 1:
            raise ValueError('initial_array must contain at least one sample.')
        if delta_t is None:
            try:
                delta_t = initial_array.delta_t
            except AttributeError:
                raise TypeError('must provide either an initial_array with a delta_t attribute, or a value for delta_t')
        if not delta_t > 0:
            raise ValueError('delta_t must be a positive number')

        # Get epoch from initial_array if epoch not given (or is None)
        # If initialy array has no epoch, set epoch to 0.
        # If epoch is provided, use that.
        if not isinstance(epoch, _lal.LIGOTimeGPS):
            if epoch is None:
                if isinstance(initial_array, TimeSeries):
                    epoch = initial_array._epoch
                else:
                    epoch = _lal.LIGOTimeGPS(0)
            elif epoch is not None:
                try:
                    epoch = _lal.LIGOTimeGPS(epoch)
                except:
                    raise TypeError('epoch must be either None or a lal.LIGOTimeGPS')
        Array.__init__(self, initial_array, dtype=dtype, copy=copy)
        self._delta_t = delta_t
        self._epoch = epoch

    def to_astropy(self, name='pycbc'):
        """ Return an astropy.timeseries.TimeSeries instance
        """
        from astropy.timeseries import TimeSeries as ATimeSeries
        from astropy.time import Time
        from astropy.units import s

        start = Time(float(self.start_time), format='gps', scale='utc')
        delta = self.delta_t * s
        return ATimeSeries({name: self.numpy()},
                           time_start=start,
                           time_delta=delta,
                           n_samples=len(self))

    def epoch_close(self, other):
        """ Check if the epoch is close enough to allow operations """
        dt = abs(float(self.start_time - other.start_time))
        return dt <= 1e-7

    def sample_rate_close(self, other):
        """ Check if the sample rate is close enough to allow operations """

        # compare our delta_t either to a another time series' or
        # to a given sample rate (float)
        if isinstance(other, TimeSeries):
            odelta_t = other.delta_t
        else:
            odelta_t = 1.0/other

        if (odelta_t - self.delta_t) / self.delta_t > 1e-4:
            return False

        if abs(1 - odelta_t / self.delta_t) * len(self) > 0.5:
            return False

        return True

    def _return(self, ary):
        return TimeSeries(ary, self._delta_t, epoch=self._epoch, copy=False)

    def _typecheck(self, other):
        if isinstance(other, TimeSeries):
            if not self.sample_rate_close(other):
                raise ValueError('different delta_t, {} vs {}'.format(
                    self.delta_t, other.delta_t))
            if not self.epoch_close(other):
                raise ValueError('different epoch, {} vs {}'.format(
                    self.start_time, other.start_time))

    def _getslice(self, index):
        # Set the new epoch---note that index.start may also be None
        if index.start is None:
            new_epoch = self._epoch
        else:
            if index.start < 0:
                raise ValueError(('Negative start index ({})'
                                  ' not supported').format(index.start))
            new_epoch = self._epoch + index.start * self._delta_t

        if index.step is not None:
            new_delta_t = self._delta_t * index.step
        else:
            new_delta_t = self._delta_t

        return TimeSeries(Array._getslice(self, index), new_delta_t,
                          new_epoch, copy=False)


    def prepend_zeros(self, num):
        """Prepend num zeros onto the beginning of this TimeSeries. Update also
        epoch to include this prepending.
        """
        self.resize(len(self) + num)
        self.roll(num)
        self._epoch = self._epoch - num * self._delta_t

    def append_zeros(self, num):
        """Append num zeros onto the end of this TimeSeries.
        """
        self.resize(len(self) + num)

    def get_delta_t(self):
        """Return time between consecutive samples in seconds.
        """
        return self._delta_t
    delta_t = property(get_delta_t,
                       doc="Time between consecutive samples in seconds.")

    def get_duration(self):
        """Return duration of time series in seconds.
        """
        return len(self) * self._delta_t
    duration = property(get_duration,
                        doc="Duration of time series in seconds.")

    def get_sample_rate(self):
        """Return the sample rate of the time series.
        """
        return 1.0/self.delta_t
    sample_rate = property(get_sample_rate,
                           doc="The sample rate of the time series.")

    def time_slice(self, start, end, mode='floor'):
        """Return the slice of the time series that contains the time range
        in GPS seconds.
        """
        if start < self.start_time:
            raise ValueError('Time series does not contain a time as early as %s' % start)

        if end > self.end_time:
            raise ValueError('Time series does not contain a time as late as %s' % end)

        start_idx = float(start - self.start_time) * self.sample_rate
        end_idx = float(end - self.start_time) * self.sample_rate

        if _numpy.isclose(start_idx, round(start_idx)):
            start_idx = round(start_idx)

        if _numpy.isclose(end_idx, round(end_idx)):
            end_idx = round(end_idx)

        if mode == 'floor':
            start_idx = int(start_idx)
            end_idx = int(end_idx)
        elif mode == 'nearest':
            start_idx = int(round(start_idx))
            end_idx = int(round(end_idx))
        else:
            raise ValueError("Invalid mode: {}".format(mode))

        return self[start_idx:end_idx]

    @property
    def delta_f(self):
        """Return the delta_f this ts would have in the frequency domain
        """
        return 1.0 / self.duration

    @property
    def start_time(self):
        """Return time series start time as a LIGOTimeGPS.
        """
        return self._epoch

    @start_time.setter
    def start_time(self, time):
        """ Set the start time
        """
        self._epoch = _lal.LIGOTimeGPS(time)

    def get_end_time(self):
        """Return time series end time as a LIGOTimeGPS.
        """
        return self._epoch + self.get_duration()
    end_time = property(get_end_time,
                        doc="Time series end time as a LIGOTimeGPS.")

    def get_sample_times(self):
        """Return an Array containing the sample times.
        """
        if self._epoch is None:
            return Array(range(len(self))) * self._delta_t
        else:
            return Array(range(len(self))) * self._delta_t + float(self._epoch)
    sample_times = property(get_sample_times,
                            doc="Array containing the sample times.")

    def at_time(self, time, nearest_sample=False, interpolate=None):
        """ Return the value at the specified gps time

        Parameters
        ----------
        nearest_sample: bool
            Return the sample at the time nearest to the chosen time rather
            than rounded down.
        interpolate: str, None
            Return the interpolated value of the time series. Choices
            are simple linear or quadratic interpolation.
        """
        if interpolate == 'linear':
            i = (time - float(self.start_time))*self.sample_rate
            di = i - int(i)
            i = int(i)
            a = self[i]
            b = self[i+1]
            return a + (b - a) * di
        elif interpolate == 'quadratic':
            ir = (time - float(self.start_time))*self.sample_rate
            i = _numpy.floor(_numpy.asarray(ir)).astype(int)
            di = ir - i
            c = self.data[i]
            xr = self.data[i + 1] - c
            xl = self.data[i - 1] - c
            a = 0.5 * (xr + xl)
            b = 0.5 * (xr - xl)
            ans = a * di**2.0 + b * di + c
            return ans

        if nearest_sample:
            time += self.delta_t / 2.0
        return self[int((time-self.start_time)*self.sample_rate)]

    def at_times(self, times, nearest_sample = False):
        """ Return an array of values at the specified gps times

        Parameters
        ----------
        times: array of floats
            The times whose values are needed
        nearest_sample: bool
            Return the samples at the times nearest to the chosen times rather
            than rounded down.

        Returns
        -------
        values: array of floats
            The values of the timeseries at the given times
        """

        if nearest_sample:
            times += self.delta_t / 2.0
        elapsed_times = times - self.start_time
        return self[(elapsed_times * self.sample_rate).astype('int')]

    def __eq__(self,other):
        """
        This is the Python special method invoked whenever the '=='
        comparison is used.  It will return true if the data of two
        time series are identical, and all of the numeric meta-data
        are identical, irrespective of whether or not the two
        instances live in the same memory (for that comparison, the
        Python statement 'a is b' should be used instead).

        Thus, this method returns 'True' if the types of both 'self'
        and 'other' are identical, as well as their lengths, dtypes,
        epochs, delta_ts and the data in the arrays, element by element.
        It will always do the comparison on the CPU, but will *not* move
        either object to the CPU if it is not already there, nor change
        the scheme of either object. It is possible to compare a CPU
        object to a GPU object, and the comparison should be true if the
        data and meta-data of the two objects are the same.

        Note in particular that this function returns a single boolean,
        and not an array of booleans as Numpy does.  If the numpy
        behavior is instead desired it can be obtained using the numpy()
        method of the PyCBC type to get a numpy instance from each
        object, and invoking '==' on those two instances.

        Parameters
        ----------
        other: another Python object, that should be tested for equality
            with 'self'.

        Returns
        -------
        boolean: 'True' if the types, dtypes, lengths, epochs, delta_ts
            and data of the two objects are each identical.
        """
        if super(TimeSeries,self).__eq__(other):
            return (self._epoch == other._epoch and self._delta_t == other._delta_t)
        else:
            return False

    def almost_equal_elem(self,other,tol,relative=True,dtol=0.0):
        """
        Compare whether two time series are almost equal, element
        by element.

        If the 'relative' parameter is 'True' (the default) then the
        'tol' parameter (which must be positive) is interpreted as a
        relative tolerance, and the comparison returns 'True' only if
        abs(self[i]-other[i]) <= tol*abs(self[i])
        for all elements of the series.

        If 'relative' is 'False', then 'tol' is an absolute tolerance,
        and the comparison is true only if
        abs(self[i]-other[i]) <= tol
        for all elements of the series.

        The method also checks that self.delta_t is within 'dtol' of
        other.delta_t; if 'dtol' has its default value of 0 then exact
        equality between the two is required.

        Other meta-data (type, dtype, length, and epoch) must be exactly
        equal.  If either object's memory lives on the GPU it will be
        copied to the CPU for the comparison, which may be slow. But the
        original object itself will not have its memory relocated nor
        scheme changed.

        Parameters
        ----------
        other: another Python object, that should be tested for
            almost-equality with 'self', element-by-element.
        tol: a non-negative number, the tolerance, which is interpreted
            as either a relative tolerance (the default) or an absolute
            tolerance.
        relative: A boolean, indicating whether 'tol' should be interpreted
            as a relative tolerance (if True, the default if this argument
            is omitted) or as an absolute tolerance (if tol is False).
        dtol: a non-negative number, the tolerance for delta_t. Like 'tol',
            it is interpreted as relative or absolute based on the value of
            'relative'.  This parameter defaults to zero, enforcing exact
            equality between the delta_t values of the two TimeSeries.

        Returns
        -------
        boolean: 'True' if the data and delta_ts agree within the tolerance,
            as interpreted by the 'relative' keyword, and if the types,
            lengths, dtypes, and epochs are exactly the same.
        """
        # Check that the delta_t tolerance is non-negative; raise an exception
        # if needed.
        if (dtol < 0.0):
            raise ValueError("Tolerance in delta_t cannot be negative")
        if super(TimeSeries,self).almost_equal_elem(other,tol=tol,relative=relative):
            if relative:
                return (self._epoch == other._epoch and
                        abs(self._delta_t-other._delta_t) <= dtol*self._delta_t)
            else:
                return (self._epoch == other._epoch and
                        abs(self._delta_t-other._delta_t) <= dtol)
        else:
            return False

    def almost_equal_norm(self,other,tol,relative=True,dtol=0.0):
        """
        Compare whether two time series are almost equal, normwise.

        If the 'relative' parameter is 'True' (the default) then the
        'tol' parameter (which must be positive) is interpreted as a
        relative tolerance, and the comparison returns 'True' only if
        abs(norm(self-other)) <= tol*abs(norm(self)).

        If 'relative' is 'False', then 'tol' is an absolute tolerance,
        and the comparison is true only if
        abs(norm(self-other)) <= tol

        The method also checks that self.delta_t is within 'dtol' of
        other.delta_t; if 'dtol' has its default value of 0 then exact
        equality between the two is required.

        Other meta-data (type, dtype, length, and epoch) must be exactly
        equal.  If either object's memory lives on the GPU it will be
        copied to the CPU for the comparison, which may be slow. But the
        original object itself will not have its memory relocated nor
        scheme changed.

        Parameters
        ----------
        other: another Python object, that should be tested for
            almost-equality with 'self', based on their norms.
        tol: a non-negative number, the tolerance, which is interpreted
            as either a relative tolerance (the default) or an absolute
            tolerance.
        relative: A boolean, indicating whether 'tol' should be interpreted
            as a relative tolerance (if True, the default if this argument
            is omitted) or as an absolute tolerance (if tol is False).
        dtol: a non-negative number, the tolerance for delta_t. Like 'tol',
            it is interpreted as relative or absolute based on the value of
            'relative'.  This parameter defaults to zero, enforcing exact
            equality between the delta_t values of the two TimeSeries.

        Returns
        -------
        boolean: 'True' if the data and delta_ts agree within the tolerance,
            as interpreted by the 'relative' keyword, and if the types,
            lengths, dtypes, and epochs are exactly the same.
        """
        # Check that the delta_t tolerance is non-negative; raise an exception
        # if needed.
        if (dtol < 0.0):
            raise ValueError("Tolerance in delta_t cannot be negative")
        if super(TimeSeries,self).almost_equal_norm(other,tol=tol,relative=relative):
            if relative:
                return (self._epoch == other._epoch and
                        abs(self._delta_t-other._delta_t) <= dtol*self._delta_t)
            else:
                return (self._epoch == other._epoch and
                        abs(self._delta_t-other._delta_t) <= dtol)
        else:
            return False

    @_convert
    def lal(self):
        """Produces a LAL time series object equivalent to self.

        Returns
        -------
        lal_data : {lal.*TimeSeries}
            LAL time series object containing the same data as self.
            The actual type depends on the sample's dtype.  If the epoch of
            self is 'None', the epoch of the returned LAL object will be
            LIGOTimeGPS(0,0); otherwise, the same as that of self.

        Raises
        ------
        TypeError
            If time series is stored in GPU memory.
        """
        lal_data = None
        ep = self._epoch

        if self._data.dtype == _numpy.float32:
            lal_data = _lal.CreateREAL4TimeSeries("",ep,0,self.delta_t,_lal.SecondUnit,len(self))
        elif self._data.dtype == _numpy.float64:
            lal_data = _lal.CreateREAL8TimeSeries("",ep,0,self.delta_t,_lal.SecondUnit,len(self))
        elif self._data.dtype == _numpy.complex64:
            lal_data = _lal.CreateCOMPLEX8TimeSeries("",ep,0,self.delta_t,_lal.SecondUnit,len(self))
        elif self._data.dtype == _numpy.complex128:
            lal_data = _lal.CreateCOMPLEX16TimeSeries("",ep,0,self.delta_t,_lal.SecondUnit,len(self))

        lal_data.data.data[:] = self.numpy()

        return lal_data

    def crop(self, left, right):
        """ Remove given seconds from either end of time series

        Parameters
        ----------
        left : float
            Number of seconds of data to remove from the left of the time series.
        right : float
            Number of seconds of data to remove from the right of the time series.

        Returns
        -------
        cropped : pycbc.types.TimeSeries
            The reduced time series
        """
        if left + right > self.duration:
            raise ValueError('Cannot crop more data than we have')

        s = int(left * self.sample_rate)
        e = len(self) - int(right * self.sample_rate)
        return self[s:e]

    def save_to_wav(self, file_name):
        """ Save this time series to a wav format audio file.

        Parameters
        ----------
        file_name : string
             The output file name
        """
        scaled = _numpy.int16(self.numpy()/max(abs(self)) * 32767)
        write_wav(file_name, int(self.sample_rate), scaled)

    def psd(self, segment_duration, **kwds):
        """ Calculate the power spectral density of this time series.

        Use the `pycbc.psd.welch` method to estimate the psd of this time segment.
        For more complete options, please see that function.

        Parameters
        ----------
        segment_duration: float
            Duration in seconds to use for each sample of the spectrum.
        kwds : keywords
            Additional keyword arguments are passed on to the `pycbc.psd.welch` method.

        Returns
        -------
        psd : FrequencySeries
            Frequency series containing the estimated PSD.
        """
        from pycbc.psd import welch
        seg_len = int(round(segment_duration * self.sample_rate))
        seg_stride = int(seg_len / 2)
        return welch(self, seg_len=seg_len,
                           seg_stride=seg_stride,
                           **kwds)

    def gate(self, time, window=0.25, method='taper', copy=True,
             taper_width=0.25, invpsd=None):
        """ Gate out portion of time series

        Parameters
        ----------
        time: float
            Central time of the gate in seconds
        window: float
            Half-length in seconds to remove data around gate time.
        method: str
            Method to apply gate, options are 'hard', 'taper', and 'paint'.
        copy: bool
            If False, do operations inplace to this time series, else return
            new time series.
        taper_width: float
            Length of tapering region on either side of excized data. Only
            applies to the taper gating method.
        invpsd: pycbc.types.FrequencySeries
            The inverse PSD to use for painting method. If not given,
            a PSD is generated using default settings.

        Returns
        -------
        data: pycbc.types.TimeSeris
            Gated time series
        """
        data = self.copy() if copy else self
        if method == 'taper':
            from pycbc.strain import gate_data
            return gate_data(data, [(time, window, taper_width)])
        elif method == 'paint':
            # Uses the hole-filling method of
            # https://arxiv.org/pdf/1908.05644.pdf
            from pycbc.strain.gate import gate_and_paint
            if invpsd is None:
                # These are some bare minimum settings, normally you
                # should probably provide a psd
                invpsd = 1. / self.filter_psd(self.duration/32, self.delta_f, 0)
            lindex = int((time - window - self.start_time) / self.delta_t)
            rindex = lindex + int(2 * window / self.delta_t)
            lindex = lindex if lindex >= 0 else 0
            rindex = rindex if rindex <= len(self) else len(self)
            return gate_and_paint(data, lindex, rindex, invpsd, copy=False)
        elif method == 'hard':
            tslice = data.time_slice(time - window, time + window)
            tslice[:] = 0
            return data
        else:
            raise ValueError('Invalid method name: {}'.format(method))

    def filter_psd(self, segment_duration, delta_f, flow):
        """ Calculate the power spectral density of this time series.

        Use the `pycbc.psd.welch` method to estimate the psd of this time segment.
        The psd is then truncated in the time domain to the segment duration
        and interpolated to the requested sample frequency.

        Parameters
        ----------
        segment_duration: float
            Duration in seconds to use for each sample of the spectrum.
        delta_f : float
            Frequency spacing to return psd at.
        flow : float
            The low frequency cutoff to apply when truncating the inverse
            spectrum.

        Returns
        -------
        psd : FrequencySeries
            Frequency series containing the estimated PSD.
        """
        from pycbc.psd import interpolate, inverse_spectrum_truncation
        p = self.psd(segment_duration)
        samples = int(round(p.sample_rate * segment_duration))
        p = interpolate(p, delta_f)
        return inverse_spectrum_truncation(p, samples,
                                           low_frequency_cutoff=flow,
                                           trunc_method='hann')

    def whiten(self, segment_duration, max_filter_duration, trunc_method='hann',
                     remove_corrupted=True, low_frequency_cutoff=None,
                     return_psd=False, **kwds):
        """ Return a whitened time series

        Parameters
        ----------
        segment_duration: float
            Duration in seconds to use for each sample of the spectrum.
        max_filter_duration : int
            Maximum length of the time-domain filter in seconds.
        trunc_method : {None, 'hann'}
            Function used for truncating the time-domain filter.
            None produces a hard truncation at `max_filter_len`.
        remove_corrupted : {True, boolean}
            If True, the region of the time series corrupted by the whitening
            is excised before returning. If false, the corrupted regions
            are not excised and the full time series is returned.
        low_frequency_cutoff : {None, float}
            Low frequency cutoff to pass to the inverse spectrum truncation.
            This should be matched to a known low frequency cutoff of the
            data if there is one.
        return_psd : {False, Boolean}
            Return the estimated and conditioned PSD that was used to whiten
            the data.
        kwds : keywords
            Additional keyword arguments are passed on to the `pycbc.psd.welch` method.

        Returns
        -------
        whitened_data : TimeSeries
            The whitened time series
        """
        from pycbc.psd import inverse_spectrum_truncation, interpolate
        # Estimate the noise spectrum
        psd = self.psd(segment_duration, **kwds)
        psd = interpolate(psd, self.delta_f)
        max_filter_len = int(round(max_filter_duration * self.sample_rate))

        # Interpolate and smooth to the desired corruption length
        psd = inverse_spectrum_truncation(psd,
                   max_filter_len=max_filter_len,
                   low_frequency_cutoff=low_frequency_cutoff,
                   trunc_method=trunc_method)

        # Whiten the data by the asd
        white = (self.to_frequencyseries() / psd**0.5).to_timeseries()

        if remove_corrupted:
            white = white[int(max_filter_len/2):int(len(self)-max_filter_len/2)]

        if return_psd:
            return white, psd

        return white

    def qtransform(self, delta_t=None, delta_f=None, logfsteps=None,
                  frange=None, qrange=(4,64), mismatch=0.2, return_complex=False):
        """ Return the interpolated 2d qtransform of this data

        Parameters
        ----------
        delta_t : {self.delta_t, float}
            The time resolution to interpolate to
        delta_f : float, Optional
            The frequency resolution to interpolate to
        logfsteps : int
            Do a log interpolation (incompatible with delta_f option) and set
            the number of steps to take.
        frange : {(30, nyquist*0.8), tuple of ints}
            frequency range
        qrange : {(4, 64), tuple}
            q range
        mismatch : float
            Mismatch between frequency tiles
        return_complex: {False, bool}
            return the raw complex series instead of the normalized power.

        Returns
        -------
        times : numpy.ndarray
            The time that the qtransform is sampled.
        freqs : numpy.ndarray
            The frequencies that the qtransform is sampled.
        qplane : numpy.ndarray (2d)
            The two dimensional interpolated qtransform of this time series.
        """
        from pycbc.filter.qtransform import qtiling, qplane
        from scipy.interpolate import interp2d

        if frange is None:
            frange = (30, int(self.sample_rate / 2 * 8))

        q_base = qtiling(self, qrange, frange, mismatch)
        _, times, freqs, q_plane = qplane(q_base, self.to_frequencyseries(),
                                          return_complex=return_complex)
        if logfsteps and delta_f:
            raise ValueError("Provide only one (or none) of delta_f and logfsteps")

        # Interpolate if requested
        if delta_f or delta_t or logfsteps:
            if return_complex:
                interp_amp = interp2d(times, freqs, abs(q_plane))
                interp_phase = interp2d(times, freqs, _numpy.angle(q_plane))
            else:
                interp = interp2d(times, freqs, q_plane)

        if delta_t:
            times = _numpy.arange(float(self.start_time),
                                    float(self.end_time), delta_t)
        if delta_f:
            freqs = _numpy.arange(int(frange[0]), int(frange[1]), delta_f)
        if logfsteps:
            freqs = _numpy.logspace(_numpy.log10(frange[0]),
                                    _numpy.log10(frange[1]),
                                     logfsteps)

        if delta_f or delta_t or logfsteps:
            if return_complex:
                q_plane = _numpy.exp(1.0j * interp_phase(times, freqs))
                q_plane *= interp_amp(times, freqs)
            else:
                q_plane = interp(times, freqs)

        return times, freqs, q_plane

    def notch_fir(self, f1, f2, order, beta=5.0, remove_corrupted=True):
        """ notch filter the time series using an FIR filtered generated from
        the ideal response passed through a time-domain kaiser
        window (beta = 5.0)

        The suppression of the notch filter is related to the bandwidth and
        the number of samples in the filter length. For a few Hz bandwidth,
        a length corresponding to a few seconds is typically
        required to create significant suppression in the notched band.

        Parameters
        ----------
        Time Series: TimeSeries
            The time series to be notched.
        f1: float
            The start of the frequency suppression.
        f2: float
            The end of the frequency suppression.
        order: int
            Number of corrupted samples on each side of the time series
        beta: float
            Beta parameter of the kaiser window that sets the side lobe attenuation.
        """
        from pycbc.filter import notch_fir
        ts = notch_fir(self, f1, f2, order, beta=beta)
        if remove_corrupted:
            ts = ts[order:len(ts)-order]
        return ts

    def lowpass_fir(self, frequency, order, beta=5.0, remove_corrupted=True):
        """ Lowpass filter the time series using an FIR filtered generated from
        the ideal response passed through a kaiser window (beta = 5.0)

        Parameters
        ----------
        Time Series: TimeSeries
            The time series to be low-passed.
        frequency: float
            The frequency below which is suppressed.
        order: int
            Number of corrupted samples on each side of the time series
        beta: float
            Beta parameter of the kaiser window that sets the side lobe attenuation.
        remove_corrupted : {True, boolean}
            If True, the region of the time series corrupted by the filtering
            is excised before returning. If false, the corrupted regions
            are not excised and the full time series is returned.
        """
        from pycbc.filter import lowpass_fir
        ts = lowpass_fir(self, frequency, order, beta=beta)
        if remove_corrupted:
            ts = ts[order:len(ts)-order]
        return ts

    def highpass_fir(self, frequency, order, beta=5.0, remove_corrupted=True):
        """ Highpass filter the time series using an FIR filtered generated from
        the ideal response passed through a kaiser window (beta = 5.0)

        Parameters
        ----------
        Time Series: TimeSeries
            The time series to be high-passed.
        frequency: float
            The frequency below which is suppressed.
        order: int
            Number of corrupted samples on each side of the time series
        beta: float
            Beta parameter of the kaiser window that sets the side lobe attenuation.
        remove_corrupted : {True, boolean}
            If True, the region of the time series corrupted by the filtering
            is excised before returning. If false, the corrupted regions
            are not excised and the full time series is returned.
        """
        from pycbc.filter import highpass_fir
        ts = highpass_fir(self, frequency, order, beta=beta)
        if remove_corrupted:
            ts = ts[order:len(ts)-order]
        return ts

    def fir_zero_filter(self, coeff):
        """Filter the timeseries with a set of FIR coefficients

        Parameters
        ----------
        coeff: numpy.ndarray
            FIR coefficients. Should be and odd length and symmetric.

        Returns
        -------
        filtered_series: pycbc.types.TimeSeries
            Return the filtered timeseries, which has been properly shifted to account
        for the FIR filter delay and the corrupted regions zeroed out.
        """
        from pycbc.filter import fir_zero_filter
        return self._return(fir_zero_filter(coeff, self))

    def resample(self, delta_t):
        """ Resample this time series to the new delta_t

        Parameters
        -----------
        delta_t: float
            The time step to resample the times series to.

        Returns
        -------
        resampled_ts: pycbc.types.TimeSeries
            The resample timeseries at the new time interval delta_t.
        """
        from pycbc.filter import resample_to_delta_t
        return resample_to_delta_t(self, delta_t)

    def save(self, path, group = None):
        """
        Save time series to a Numpy .npy, hdf, or text file. The first column
        contains the sample times, the second contains the values.
        In the case of a complex time series saved as text, the imaginary
        part is written as a third column. When using hdf format, the data is stored
        as a single vector, along with relevant attributes.

        Parameters
        ----------
        path: string
            Destination file path. Must end with either .hdf, .npy or .txt.

        group: string
            Additional name for internal storage use. Ex. hdf storage uses
            this as the key value.

        Raises
        ------
        ValueError
            If path does not end in .npy or .txt.
        """

        ext = _os.path.splitext(path)[1]
        if ext == '.npy':
            output = _numpy.vstack((self.sample_times.numpy(), self.numpy())).T
            _numpy.save(path, output)
        elif ext == '.txt':
            if self.kind == 'real':
                output = _numpy.vstack((self.sample_times.numpy(),
                                        self.numpy())).T
            elif self.kind == 'complex':
                output = _numpy.vstack((self.sample_times.numpy(),
                                        self.numpy().real,
                                        self.numpy().imag)).T
            _numpy.savetxt(path, output)
        elif ext =='.hdf':
            key = 'data' if group is None else group
            with h5py.File(path, 'a') as f:
                ds = f.create_dataset(key, data=self.numpy(),
                                      compression='gzip',
                                      compression_opts=9, shuffle=True)
                ds.attrs['start_time'] = float(self.start_time)
                ds.attrs['delta_t'] = float(self.delta_t)
        else:
            raise ValueError('Path must end with .npy, .txt or .hdf')

    def to_timeseries(self):
        """ Return time series"""
        return self

    @_nocomplex
    def to_frequencyseries(self, delta_f=None):
        """ Return the Fourier transform of this time series

        Parameters
        ----------
        delta_f : {None, float}, optional
            The frequency resolution of the returned frequency series. By
        default the resolution is determined by the duration of the timeseries.

        Returns
        -------
        FrequencySeries:
            The fourier transform of this time series.
        """
        from pycbc.fft import fft
        if not delta_f:
            delta_f = 1.0 / self.duration

        # add 0.5 to round integer
        tlen  = int(1.0 / delta_f / self.delta_t + 0.5)
        flen = int(tlen / 2 + 1)

        if tlen < len(self):
            raise ValueError("The value of delta_f (%s) would be "
                             "undersampled. Maximum delta_f "
                             "is %s." % (delta_f, 1.0 / self.duration))
        if not delta_f:
            tmp = self
        else:
            tmp = TimeSeries(zeros(tlen, dtype=self.dtype),
                             delta_t=self.delta_t, epoch=self.start_time)
            tmp[:len(self)] = self[:]

        f = FrequencySeries(zeros(flen,
                           dtype=complex_same_precision_as(self)),
                           delta_f=delta_f)
        fft(tmp, f)
        f._delta_f = delta_f
        return f

    def inject(self, other, copy=True):
        """Return copy of self with other injected into it.

        The other vector will be resized and time shifted with sub-sample
        precision before adding. This assumes that one can assume zeros
        outside of the original vector range.
        """
        # only handle equal sample rate for now.
        if not self.sample_rate_close(other):
            raise ValueError('Sample rate must be the same')
        # determine if we want to inject in place or not
        if copy:
            ts = self.copy()
        else:
            ts = self
        # Other is disjoint
        if ((other.start_time >= ts.end_time) or
           (ts.start_time > other.end_time)):
            return ts

        other = other.copy()
        dt = float((other.start_time - ts.start_time) * ts.sample_rate)

        # This coaligns other to the time stepping of self
        if not dt.is_integer():
            diff = (dt - _numpy.floor(dt)) * ts.delta_t

            # insert zeros at end
            other.resize(len(other) + (len(other) + 1) % 2 + 1)

            # fd shift to the right
            other = other.cyclic_time_shift(diff)

        # get indices of other with respect to self
        # this is already an integer to floating point precission
        left = float(other.start_time - ts.start_time) * ts.sample_rate
        left = int(round(left))
        right = left + len(other)

        oleft = 0
        oright = len(other)

        # other overhangs on left so truncate
        if left < 0:
            oleft = -left
            left = 0

        # other overhangs on right so truncate
        if right > len(ts):
            oright = len(other) - (right - len(ts))
            right = len(ts)

        ts[left:right] += other[oleft:oright]
        return ts

    add_into = inject  # maintain backwards compatibility for now

    @_nocomplex
    def cyclic_time_shift(self, dt):
        """Shift the data and timestamps by a given number of seconds

        Shift the data and timestamps in the time domain a given number of
        seconds. To just change the time stamps, do ts.start_time += dt.
        The time shift may be smaller than the intrinsic sample rate of the data.
        Note that data will be cyclically rotated, so if you shift by 2
        seconds, the final 2 seconds of your data will now be at the
        beginning of the data set.

        Parameters
        ----------
        dt : float
            Amount of time to shift the vector.

        Returns
        -------
        data : pycbc.types.TimeSeries
            The time shifted time series.
        """
        # We do this in the frequency domain to allow us to do sub-sample
        # time shifts. This also results in the shift being circular. It
        # is left to a future update to do a faster impelementation in the case
        # where the time shift can be done with an exact number of samples.
        return self.to_frequencyseries().cyclic_time_shift(dt).to_timeseries()

    def match(self, other, psd=None,
              low_frequency_cutoff=None, high_frequency_cutoff=None):
        """ Return the match between the two TimeSeries or FrequencySeries.

        Return the match between two waveforms. This is equivalent to the overlap
        maximized over time and phase. By default, the other vector will be
        resized to match self. This may remove high frequency content or the
        end of the vector.

        Parameters
        ----------
        other : TimeSeries or FrequencySeries
            The input vector containing a waveform.
        psd : Frequency Series
            A power spectral density to weight the overlap.
        low_frequency_cutoff : {None, float}, optional
            The frequency to begin the match.
        high_frequency_cutoff : {None, float}, optional
            The frequency to stop the match.

        Returns
        -------
        match: float
        index: int
            The number of samples to shift to get the match.
        """
        return self.to_frequencyseries().match(other, psd=psd,
                     low_frequency_cutoff=low_frequency_cutoff,
                     high_frequency_cutoff=high_frequency_cutoff)

    def detrend(self, type='linear'):
        """ Remove linear trend from the data

        Remove a linear trend from the data to improve the approximation that
        the data is circularly convolved, this helps reduce the size of filter
        transients from a circular convolution / filter.

        Parameters
        ----------
        type: str
            The choice of detrending. The default ('linear') removes a linear
        least squares fit. 'constant' removes only the mean of the data.
        """
        from scipy.signal import detrend
        return self._return(detrend(self.numpy(), type=type))

    def plot(self, **kwds):
        """ Basic plot of this time series
        """
        from matplotlib import pyplot

        if self.kind == 'real':
            plot = pyplot.plot(self.sample_times, self, **kwds)
            return plot
        elif self.kind == 'complex':
            plot1 = pyplot.plot(self.sample_times, self.real(), **kwds)
            plot2 = pyplot.plot(self.sample_times, self.imag(), **kwds)
            return plot1, plot2

def load_timeseries(path, group=None):
    """Load a TimeSeries from an HDF5, ASCII or Numpy file. The file type is
    inferred from the file extension, which must be `.hdf`, `.txt` or `.npy`.

    For ASCII and Numpy files, the first column of the array is assumed to
    contain the sample times. If the array has two columns, a real-valued time
    series is returned. If the array has three columns, the second and third
    ones are assumed to contain the real and imaginary parts of a complex time
    series.

    For HDF files, the dataset is assumed to contain the attributes `delta_t`
    and `start_time`, which should contain respectively the sampling period in
    seconds and the start GPS time of the data.

    The default data types will be double precision floating point.

    Parameters
    ----------
    path: string
        Input file path. Must end with either `.npy`, `.txt` or `.hdf`.

    group: string
        Additional name for internal storage use. When reading HDF files, this
        is the path to the HDF dataset to read.

    Raises
    ------
    ValueError
        If path does not end in a supported extension.
        For Numpy and ASCII input files, this is also raised if the array
        does not have 2 or 3 dimensions.
    """
    ext = _os.path.splitext(path)[1]
    if ext == '.npy':
        data = _numpy.load(path)
    elif ext == '.txt':
        data = _numpy.loadtxt(path)
    elif ext == '.hdf':
        key = 'data' if group is None else group
        with h5py.File(path, 'r') as f:
            data = f[key][:]
            series = TimeSeries(data, delta_t=f[key].attrs['delta_t'],
                                epoch=f[key].attrs['start_time'])
        return series
    else:
        raise ValueError('Path must end with .npy, .hdf, or .txt')

    delta_t = (data[-1][0] - data[0][0]) / (len(data) - 1)
    epoch = _lal.LIGOTimeGPS(data[0][0])
    if data.ndim == 2:
        return TimeSeries(data[:,1], delta_t=delta_t, epoch=epoch)
    elif data.ndim == 3:
        return TimeSeries(data[:,1] + 1j*data[:,2],
                          delta_t=delta_t, epoch=epoch)

    raise ValueError('File has %s dimensions, cannot convert to TimeSeries, \
                      must be 2 (real) or 3 (complex)' % data.ndim)
