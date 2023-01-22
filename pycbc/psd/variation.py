""" PSD Variation """

import numpy
from numpy.fft import rfft, irfft
import scipy.signal as sig


import pycbc.psd
from pycbc.types import TimeSeries
from pycbc.filter import resample_to_delta_t


def mean_square(data, delta_t, srate, short_stride, stride):
    """ Calculate mean square of given time series once per stride

    First of all this function calculate the mean square of given time
    series once per short_stride. This is used to find and remove
    outliers due to short glitches. Here an outlier is defined as any
    element which is greater than two times the average of its closest
    neighbours. Every outlier is substituted with the average of the
    corresponding adjacent elements.
    Then, every second the function compute the mean square of the
    smoothed time series, within the stride.

    Parameters
    ----------
    data : numpy.ndarray
    delta_t : float
        Duration of the time series
    srate : int
        Sample rate of the data were it given as a TimeSeries
    short_stride : float
        Stride duration for outlier removal
    stride ; float
        Stride duration

    Returns
    -------
    m_s: List
        Mean square of given time series
    """

    # Calculate mean square of data once per short stride and replace
    # outliers
    short_ms = numpy.mean(data.reshape(-1, int(srate * short_stride)) ** 2,
                          axis=1)
    # Define an array of averages that is used to substitute outliers
    ave = 0.5 * (short_ms[2:] + short_ms[:-2])
    outliers = short_ms[1:-1] > (2. * ave)
    short_ms[1:-1][outliers] = ave[outliers]

    # Calculate mean square of data every step within a window equal to
    # stride seconds
    m_s = []
    inv_time = int(1. / short_stride)
    for index in range(int(delta_t - stride + 1)):
        m_s.append(numpy.mean(short_ms[inv_time * index:inv_time *
                                       int(index+stride)]))
    return m_s


def calc_filt_psd_variation(strain, segment, short_segment, psd_long_segment,
                            psd_duration, psd_stride, psd_avg_method, low_freq,
                            high_freq):
    """ Calculates time series of PSD variability

    This function first splits the segment up into 512 second chunks. It
    then calculates the PSD over this 512 second. The PSD is used to
    to create a filter that is the composition of three filters:
    1. Bandpass filter between f_low and f_high.
    2. Weighting filter which gives the rough response of a CBC template.
    3. Whitening filter.
    Next it makes the convolution of this filter with the stretch of data.
    This new time series is given to the "mean_square" function, which
    computes the mean square of the timeseries within an 8 seconds window,
    once per second.
    The result, which is the variance of the S/N in that stride for the
    Parseval theorem, is then stored in a timeseries.

    Parameters
    ----------
    strain : TimeSeries
        Input strain time series to estimate PSDs
    segment : {float, 8}
        Duration of the segments for the mean square estimation in seconds.
    short_segment : {float, 0.25}
        Duration of the short segments for the outliers removal.
    psd_long_segment : {float, 512}
        Duration of the long segments for PSD estimation in seconds.
    psd_duration : {float, 8}
        Duration of FFT segments for long term PSD estimation, in seconds.
    psd_stride : {float, 4}
        Separation between FFT segments for long term PSD estimation, in
        seconds.
    psd_avg_method : {string, 'median'}
        Method for averaging PSD estimation segments.
    low_freq : {float, 20}
        Minimum frequency to consider the comparison between PSDs.
    high_freq : {float, 480}
        Maximum frequency to consider the comparison between PSDs.

    Returns
    -------
    psd_var : TimeSeries
        Time series of the variability in the PSD estimation
    """
    # Calculate strain precision
    if strain.precision == 'single':
        fs_dtype = numpy.float32
    elif strain.precision == 'double':
        fs_dtype = numpy.float64

    # Convert start and end times immediately to floats
    start_time = float(strain.start_time)
    end_time = float(strain.end_time)

    # Resample the data
    strain = resample_to_delta_t(strain, 1.0 / 2048)
    srate = int(strain.sample_rate)

    # Fix the step for the PSD estimation and the time to remove at the
    # edge of the time series.
    step = 1.0
    strain_crop = 8.0

    # Find the times of the long segments
    times_long = numpy.arange(start_time, end_time,
                              psd_long_segment - 2 * strain_crop
                              - segment + step)

    # Create a bandpass filter between low_freq and high_freq
    filt = sig.firwin(4 * srate, [low_freq, high_freq], pass_zero=False,
                      window='hann', nyq=srate / 2)
    filt.resize(int(psd_duration * srate))
    # Fourier transform the filter and take the absolute value to get
    # rid of the phase.
    filt = abs(rfft(filt))

    psd_var_list = []
    for tlong in times_long:
        # Calculate PSD for long segment
        if tlong + psd_long_segment <= float(end_time):
            astrain = strain.time_slice(tlong, tlong + psd_long_segment)
            plong = pycbc.psd.welch(
                astrain,
                seg_len=int(psd_duration * strain.sample_rate),
                seg_stride=int(psd_stride * strain.sample_rate),
                avg_method=psd_avg_method)
        else:
            astrain = strain.time_slice(tlong, end_time)
            plong = pycbc.psd.welch(
                           strain.time_slice(end_time - psd_long_segment,
                                             end_time),
                           seg_len=int(psd_duration * strain.sample_rate),
                           seg_stride=int(psd_stride * strain.sample_rate),
                           avg_method=psd_avg_method)
        astrain = astrain.numpy()
        freqs = numpy.array(plong.sample_frequencies, dtype=fs_dtype)
        plong = plong.numpy()

        # Make the weighting filter - bandpass, which weight by f^-7/6,
        # and whiten. The normalization is chosen so that the variance
        # will be one if this filter is applied to white noise which
        # already has a variance of one.
        fweight = freqs ** (-7./6.) * filt / numpy.sqrt(plong)
        fweight[0] = 0.
        norm = (sum(abs(fweight) ** 2) / (len(fweight) - 1.)) ** -0.5
        fweight = norm * fweight
        fwhiten = numpy.sqrt(2. / srate) / numpy.sqrt(plong)
        fwhiten[0] = 0.
        full_filt = sig.hann(int(psd_duration * srate)) * numpy.roll(
            irfft(fwhiten * fweight), int(psd_duration / 2) * srate)
        # Convolve the filter with long segment of data
        wstrain = sig.fftconvolve(astrain, full_filt, mode='same')
        wstrain = wstrain[int(strain_crop * srate):-int(strain_crop * srate)]
        # compute the mean square of the chunk of data
        delta_t = len(wstrain) * strain.delta_t
        variation = mean_square(wstrain, delta_t, srate, short_segment, segment)
        psd_var_list.append(numpy.array(variation, dtype=wstrain.dtype))

    # Package up the time series to return
    psd_var = TimeSeries(numpy.concatenate(psd_var_list), delta_t=step,
                         epoch=start_time + strain_crop + segment)

    return psd_var


def find_trigger_value(psd_var, idx, start, sample_rate):
    """ Find the PSD variation value at a particular time with the filter
    method. If the time is outside the timeseries bound, 1. is given.

    Parameters
    ----------
    psd_var : TimeSeries
        Time series of the varaibility in the PSD estimation
    idx : numpy.ndarray
        Time indices of the triggers
    start : float
        GPS start time
    sample_rate : float
        Sample rate defined in ini file

    Returns
    -------
    vals : Array
        PSD variation value at a particular time
    """
    # Find gps time of the trigger
    time = start + idx / sample_rate
    # Extract the PSD variation at trigger time through linear
    # interpolation
    if not hasattr(psd_var, 'cached_psd_var_interpolant'):
        from scipy import interpolate
        psd_var.cached_psd_var_interpolant = \
            interpolate.interp1d(psd_var.sample_times.numpy(), psd_var.numpy(),
                                 fill_value=1.0, bounds_error=False)
    vals = psd_var.cached_psd_var_interpolant(time)

    return vals
