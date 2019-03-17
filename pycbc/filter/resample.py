# Copyright (C) 2012  Alex Nitz
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
import lal
import numpy
import scipy.signal
from pycbc.types import TimeSeries, Array, zeros, FrequencySeries, real_same_precision_as
from pycbc.types import complex_same_precision_as
from pycbc.fft import ifft, fft

_resample_func = {numpy.dtype('float32'): lal.ResampleREAL4TimeSeries,
                 numpy.dtype('float64'): lal.ResampleREAL8TimeSeries}

def lfilter(coefficients, timeseries):
    """ Apply filter coefficients to a time series

    Parameters
    ----------
    coefficients: numpy.ndarray
        Filter coefficients to apply
    timeseries: numpy.ndarray
        Time series to be filtered.

    Returns
    -------
    tseries: numpy.ndarray
        filtered array
    """
    from pycbc.filter import correlate

    # If there aren't many points just use the default scipy method
    if len(timeseries) < 2**7:
        if hasattr(timeseries, 'numpy'):
            timeseries = timeseries.numpy()
        series = scipy.signal.lfilter(coefficients, 1.0, timeseries)
        return series
    else:
        cseries = (Array(coefficients[::-1] * 1)).astype(timeseries.dtype)
        cseries.resize(len(timeseries))
        cseries.roll(len(timeseries) - len(coefficients) + 1)
        timeseries = Array(timeseries, copy=False)

        flen = len(cseries) / 2 + 1
        ftype = complex_same_precision_as(timeseries)

        cfreq = zeros(flen, dtype=ftype)
        tfreq = zeros(flen, dtype=ftype)

        fft(Array(cseries), cfreq)
        fft(Array(timeseries), tfreq)

        cout = zeros(flen, ftype)
        out = zeros(len(timeseries), dtype=timeseries)

        correlate(cfreq, tfreq, cout)
        ifft(cout, out)

        return out.numpy()  / len(out)

def fir_zero_filter(coeff, timeseries):
    """Filter the timeseries with a set of FIR coefficients

    Parameters
    ----------
    coeff: numpy.ndarray
        FIR coefficients. Should be and odd length and symmetric.
    timeseries: pycbc.types.TimeSeries
        Time series to be filtered.

    Returns
    -------
    filtered_series: pycbc.types.TimeSeries
        Return the filtered timeseries, which has been properly shifted to account
    for the FIR filter delay and the corrupted regions zeroed out.
    """
    # apply the filter
    series = lfilter(coeff, timeseries.numpy())

    # reverse the time shift caused by the filter,
    # corruption regions contain zeros
    # If the number of filter coefficients is odd, the central point *should*
    # be included in the output so we only zero out a region of len(coeff) - 1
    data = numpy.zeros(len(timeseries))
    data[len(coeff)//2:len(data)-len(coeff)//2] = series[(len(coeff) // 2) * 2:]
    return data

def resample_to_delta_t(timeseries, delta_t, method='butterworth'):
    """Resmple the time_series to delta_t

    Resamples the TimeSeries instance time_series to the given time step,
    delta_t. Only powers of two and real valued time series are supported
    at this time. Additional restrictions may apply to particular filter
    methods.

    Parameters
    ----------
    time_series: TimeSeries
        The time series to be resampled
    delta_t: float
        The desired time step

    Returns
    -------
    Time Series: TimeSeries
        A TimeSeries that has been resampled to delta_t.

    Raises
    ------
    TypeError:
        time_series is not an instance of TimeSeries.
    TypeError:
        time_series is not real valued

    Examples
    --------

    >>> h_plus_sampled = resample_to_delta_t(h_plus, 1.0/2048)
    """
    if not isinstance(timeseries,TimeSeries):
        raise TypeError("Can only resample time series")

    if timeseries.kind is not 'real':
        raise TypeError("Time series must be real")

    if timeseries.delta_t == delta_t:
        return timeseries * 1

    if method == 'butterworth':
        lal_data = timeseries.lal()
        _resample_func[timeseries.dtype](lal_data, delta_t)
        data = lal_data.data.data

    elif method == 'ldas':
        factor = int(delta_t / timeseries.delta_t)
        numtaps = factor * 20 + 1

        # The kaiser window has been testing using the LDAS implementation
        # and is in the same configuration as used in the original lalinspiral
        filter_coefficients = scipy.signal.firwin(numtaps, 1.0 / factor,
                                                  window=('kaiser', 5))

        # apply the filter and decimate
        data = fir_zero_filter(filter_coefficients, timeseries)[::factor]

    else:
        raise ValueError('Invalid resampling method: %s' % method)

    ts = TimeSeries(data, delta_t = delta_t,
                      dtype=timeseries.dtype,
                      epoch=timeseries._epoch)

    # From the construction of the LDAS FIR filter there will be 10 corrupted samples
    # explanation here http://software.ligo.org/docs/lalsuite/lal/group___resample_time_series__c.html
    ts.corrupted_samples = 10
    return ts


_highpass_func = {numpy.dtype('float32'): lal.HighPassREAL4TimeSeries,
                 numpy.dtype('float64'): lal.HighPassREAL8TimeSeries}

def notch_fir(timeseries, f1, f2, order, beta=5.0):
    """ notch filter the time series using an FIR filtered generated from
    the ideal response passed through a time-domain kaiser window (beta = 5.0)

    The suppression of the notch filter is related to the bandwidth and
    the number of samples in the filter length. For a few Hz bandwidth,
    a length corresponding to a few seconds is typically
    required to create significant suppression in the notched band.
    To achieve frequency resolution df at sampling frequency fs,
    order should be at least fs/df.

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
        (Extent of the filter on either side of zero)
    beta: float
        Beta parameter of the kaiser window that sets the side lobe attenuation.
    """
    k1 = f1 / float((int(1.0 / timeseries.delta_t) / 2))
    k2 = f2 / float((int(1.0 / timeseries.delta_t) / 2))
    coeff = scipy.signal.firwin(order * 2 + 1, [k1, k2], window=('kaiser', beta))
    data = fir_zero_filter(coeff, timeseries)
    return TimeSeries(data, epoch=timeseries.start_time, delta_t=timeseries.delta_t)

def lowpass_fir(timeseries, frequency, order, beta=5.0):
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
    """
    k = frequency / float((int(1.0 / timeseries.delta_t) / 2))
    coeff = scipy.signal.firwin(order * 2 + 1, k, window=('kaiser', beta))
    data = fir_zero_filter(coeff, timeseries)
    return TimeSeries(data, epoch=timeseries.start_time, delta_t=timeseries.delta_t)

def highpass_fir(timeseries, frequency, order, beta=5.0):
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
    """
    k = frequency / float((int(1.0 / timeseries.delta_t) / 2))
    coeff = scipy.signal.firwin(order * 2 + 1, k, window=('kaiser', beta), pass_zero=False)
    data = fir_zero_filter(coeff, timeseries)
    return TimeSeries(data, epoch=timeseries.start_time, delta_t=timeseries.delta_t)

def highpass(timeseries, frequency, filter_order=8, attenuation=0.1):
    """Return a new timeseries that is highpassed.

    Return a new time series that is highpassed above the `frequency`.

    Parameters
    ----------
    Time Series: TimeSeries
        The time series to be high-passed.
    frequency: float
        The frequency below which is suppressed.
    filter_order: {8, int}, optional
        The order of the filter to use when high-passing the time series.
    attenuation: {0.1, float}, optional
        The attenuation of the filter.

    Returns
    -------
    Time Series: TimeSeries
        A  new TimeSeries that has been high-passed.

    Raises
    ------
    TypeError:
        time_series is not an instance of TimeSeries.
    TypeError:
        time_series is not real valued

    """

    if not isinstance(timeseries, TimeSeries):
        raise TypeError("Can only resample time series")

    if timeseries.kind is not 'real':
        raise TypeError("Time series must be real")

    lal_data = timeseries.lal()
    _highpass_func[timeseries.dtype](lal_data, frequency,
                                     1-attenuation, filter_order)

    return TimeSeries(lal_data.data.data, delta_t = lal_data.deltaT,
                      dtype=timeseries.dtype, epoch=timeseries._epoch)

def interpolate_complex_frequency(series, delta_f, zeros_offset=0, side='right'):
    """Interpolate complex frequency series to desired delta_f.

    Return a new complex frequency series that has been interpolated to the
    desired delta_f.

    Parameters
    ----------
    series : FrequencySeries
        Frequency series to be interpolated.
    delta_f : float
        The desired delta_f of the output
    zeros_offset : optional, {0, int}
        Number of sample to delay the start of the zero padding
    side : optional, {'right', str}
        The side of the vector to zero pad

    Returns
    -------
    interpolated series : FrequencySeries
        A new FrequencySeries that has been interpolated.
    """
    new_n = int( (len(series)-1) * series.delta_f / delta_f + 1)
    old_N = int( (len(series)-1) * 2 )
    new_N = int( (new_n - 1) * 2 )
    time_series = TimeSeries(zeros(old_N), delta_t =1.0/(series.delta_f*old_N),
                             dtype=real_same_precision_as(series))

    ifft(series, time_series)

    time_series.roll(-zeros_offset)
    time_series.resize(new_N)

    if side == 'left':
        time_series.roll(zeros_offset + new_N - old_N)
    elif side == 'right':
        time_series.roll(zeros_offset)

    out_series = FrequencySeries(zeros(new_n), epoch=series.epoch,
                           delta_f=delta_f, dtype=series.dtype)
    fft(time_series, out_series)

    return out_series

__all__ = ['resample_to_delta_t', 'highpass', 'interpolate_complex_frequency', 'highpass_fir', 'lowpass_fir', 'notch_fir', 'fir_zero_filter']

