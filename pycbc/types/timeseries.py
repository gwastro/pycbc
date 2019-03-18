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
from __future__ import division
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

    Attributes
    ----------
    delta_t
    duration
    start_time
    end_time
    sample_times
    sample_rate
    """

    def __init__(self, initial_array, delta_t=None, epoch="", dtype=None, copy=True):
        if len(initial_array) < 1:
            raise ValueError('initial_array must contain at least one sample.')
        if delta_t is None:
            try:
                delta_t = initial_array.delta_t
            except AttributeError:
                raise TypeError('must provide either an initial_array with a delta_t attribute, or a value for delta_t')
        if not delta_t > 0:
            raise ValueError('delta_t must be a positive number')
        # We gave a nonsensical default value to epoch so we can test if it's been set.
        # If the user passes in an initial_array that has an 'epoch' attribute and doesn't
        # pass in a value of epoch, then our new object's epoch comes from initial_array.
        # But if the user passed in a value---even 'None'---that will take precedence over
        # anything set in initial_array.  Finally, if the user passes in something without
        # an epoch attribute *and* doesn't pass in a value of epoch, it becomes 'None'
        if not isinstance(epoch,_lal.LIGOTimeGPS):
            if epoch == "":
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

    def _return(self, ary):
        return TimeSeries(ary, self._delta_t, epoch=self._epoch, copy=False)

    def _typecheck(self, other):
        if isinstance(other, TimeSeries):
            if other._delta_t != self._delta_t:
                raise ValueError('different delta_t')
            if self._epoch != other._epoch:
                raise ValueError('different epoch')

    def _getslice(self, index):
        # Set the new epoch---note that index.start may also be None
        if self._epoch is None:
            new_epoch = None
        elif index.start is None:
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
        return int(round(1.0/self.delta_t))
    sample_rate = property(get_sample_rate,
                           doc="The sample rate of the time series.")

    def time_slice(self, start, end):
        """Return the slice of the time series that contains the time range
        in GPS seconds.
        """
        if start < self.start_time:
            raise ValueError('Time series does not contain a time as early as %s' % start)

        if end > self.end_time:
            raise ValueError('Time series does not contain a time as late as %s' % end)

        start_idx = int((start - self.start_time) * self.sample_rate)
        end_idx = int((end - self.start_time) * self.sample_rate)
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

    def at_time(self, time, nearest_sample=False):
        """ Return the value at the specified gps time
        """
        if nearest_sample:
            time += self.delta_t / 2.0
        return self[int((time-self.start_time)*self.sample_rate)]

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
        if self._epoch is None:
            ep = _lal.LIGOTimeGPS(0,0)
        else:
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
        write_wav(file_name, self.sample_rate, scaled)

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
        seg_len = int(segment_duration * self.sample_rate)
        seg_stride = int(seg_len / 2)
        return welch(self, seg_len=seg_len,
                           seg_stride=seg_stride,
                           **kwds)

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
        max_filter_len = int(max_filter_duration * self.sample_rate)

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
            f = h5py.File(path)
            ds = f.create_dataset(key, data=self.numpy(), compression='gzip',
                                  compression_opts=9, shuffle=True)
            ds.attrs['start_time'] = float(self.start_time)
            ds.attrs['delta_t'] = float(self.delta_t)
        else:
            raise ValueError('Path must end with .npy, .txt or .hdf')

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
        return f

    def add_into(self, other):
        """Return the sum of the two time series accounting for the time stamp.

        The other vector will be resized and time shifted wiht sub-sample
        precision before adding. This assumes that one can assume zeros
        outside of the original vector range.
        """
        # only handle equal sample rate for now.
        if self.sample_rate != other.sample_rate:
            raise ValueError('Sample rate must be the same')

        # Other is disjoint
        if ((other.start_time > self.end_time) or
           (self.start_time > other.end_time)):
            return self.copy()

        other = other.copy()
        dt = float((other.start_time - self.start_time) * self.sample_rate)
        if not dt.is_integer():
            diff = (dt - _numpy.floor(dt))
            other.resize(len(other) + (len(other) + 1) % 2 + 1)
            other = other.cyclic_time_shift(diff)

        ts = self.copy()
        start = max(other.start_time, self.start_time)
        end = min(other.end_time, self.end_time)
        part = ts.time_slice(start, end)
        part += other.time_slice(start, end)
        return ts

    @_nocomplex
    def cyclic_time_shift(self, dt):
        """Shift the data and timestamps by a given number of seconds

        Shift the data and timestamps in the time domain a given number of
        seconds. To just change the time stamps, do ts.start_time += dt.
        The time shift may be smaller than the intrinsic sample rate of the data.
        Note that data will be cycliclly rotated, so if you shift by 2
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

        Return the match between two waveforms. This is equivelant to the overlap
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

def load_timeseries(path, group=None):
    """
    Load a TimeSeries from a .hdf, .txt or .npy file. The
    default data types will be double precision floating point.

    Parameters
    ----------
    path: string
        source file path. Must end with either .npy or .txt.

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
        data = _numpy.load(path)
    elif ext == '.txt':
        data = _numpy.loadtxt(path)
    elif ext == '.hdf':
        key = 'data' if group is None else group
        f = h5py.File(path)
        data = f[key][:]
        series = TimeSeries(data, delta_t=f[key].attrs['delta_t'],
                                  epoch=f[key].attrs['start_time'])
        f.close()
        return series
    else:
        raise ValueError('Path must end with .npy, .hdf, or .txt')

    if data.ndim == 2:
        delta_t = (data[-1][0] - data[0][0]) / (len(data)-1)
        epoch = _lal.LIGOTimeGPS(data[0][0])
        return TimeSeries(data[:,1], delta_t=delta_t, epoch=epoch)
    elif data.ndim == 3:
        delta_t = (data[-1][0] - data[0][0]) / (len(data)-1)
        epoch = _lal.LIGOTimeGPS(data[0][0])
        return TimeSeries(data[:,1] + 1j*data[:,2],
                          delta_t=delta_t, epoch=epoch)
    else:
        raise ValueError('File has %s dimensions, cannot convert to Array, \
                          must be 2 (real) or 3 (complex)' % data.ndim)
