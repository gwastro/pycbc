""" PSD Variation """

import numpy
import pycbc.psd

from pycbc.types import TimeSeries, FrequencySeries, zeros

def calc_psd_variation(strain, psd_short_segment, psd_long_segment,
                       short_psd_duration, short_psd_stride, psd_avg_method,
                       low_freq, high_freq):
    """Calculates time series of PSD variability

    This function first splits the segment up into 512 second chunks. It
    then calculates the PSD over this 512 second period as well as in 4
    second chunks throughout each 512 second period. Next the function
    estimates how different the 4 second PSD is to the 512 second PSD and
    produces a timeseries of this variability.

    Parameters
    ----------
    strain : TimeSeries
        Input strain time series to estimate PSDs
    psd_short_segment : {float, 8}
        Duration of the short segments for PSD estimation in seconds.
    psd_long_segment : {float, 512}
        Duration of the long segments for PSD estimation in seconds.
    short_psd_duration : {float, 4}
        Duration of the segments for PSD estimation in seconds.
    short_psd_stride : {float, 2}
        Separation between PSD estimation segments in seconds.
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
    start_time = numpy.float(strain.start_time)
    end_time = numpy.float(strain.end_time)

    # Find the times of the long segments
    times_long = numpy.arange(start_time, end_time, psd_long_segment)

    # Set up the empty time series for the PSD variation estimate
    psd_var = TimeSeries(zeros(int(numpy.ceil((end_time -
                                   start_time) / psd_short_segment))),
                         delta_t=psd_short_segment, copy=False,
                         epoch=start_time)

    ind = 0
    for tlong in times_long:
        # Calculate PSD for long segment and separate the long segment in to
        # overlapping shorter segments
        if tlong + psd_long_segment <= end_time:
            psd_long = pycbc.psd.welch(
                           strain.time_slice(tlong, tlong + psd_long_segment),
                           seg_len=int(short_psd_duration * strain.sample_rate),
                           seg_stride=int(short_psd_stride *
                                          strain.sample_rate),
                           avg_method=psd_avg_method)
            times_short = numpy.arange(tlong, tlong + psd_long_segment,
                                       psd_short_segment)
        else:
            psd_long = pycbc.psd.welch(
                           strain.time_slice(end_time - psd_long_segment,
                                             end_time),
                           seg_len=int(short_psd_duration * strain.sample_rate),
                           seg_stride=int(short_psd_stride *
                                          strain.sample_rate),
                           avg_method=psd_avg_method)
            times_short = numpy.arange(tlong, end_time, psd_short_segment)

        # Calculate the PSD of the shorter segments
        psd_short = []
        for tshort in times_short:
            if tshort + psd_short_segment <= end_time:
                pshort = pycbc.psd.welch(
                            strain.time_slice(tshort, tshort +
                                              psd_short_segment),
                            seg_len=int(short_psd_duration *
                                        strain.sample_rate),
                            seg_stride=int(short_psd_stride *
                                           strain.sample_rate),
                            avg_method=psd_avg_method)
            else:
                pshort = pycbc.psd.welch(
                            strain.time_slice(tshort - psd_short_segment,
                                              end_time),
                            seg_len=int(short_psd_duration *
                                        strain.sample_rate),
                            seg_stride=int(short_psd_stride *
                                           strain.sample_rate),
                            avg_method=psd_avg_method)
            psd_short.append(pshort)

        # Estimate the range of the PSD to compare
        kmin = int(low_freq / psd_long.delta_f)
        kmax = int(high_freq / psd_long.delta_f)
        # Comapre the PSD of the short segment to the long segment
        # The weight factor gives the rough response of a cbc template across
        # the defined frequency range given the expected PSD (i.e. long PSD)
        # Then integrate the weighted ratio of the actual PSD (i.e. short PSD)
        # with the expected PSD (i.e. long PSD) over the specified frequency
        # range
        freqs = FrequencySeries(psd_long.sample_frequencies,
                                delta_f=psd_long.delta_f,
                                epoch=psd_long.epoch, dtype=fs_dtype)
        weight = numpy.array(
                     freqs[kmin:kmax]**(-7./3.) / psd_long[kmin:kmax])
        weight /= weight.sum()
        diff = numpy.array([(weight * numpy.array(p_short[kmin:kmax] /
                             psd_long[kmin:kmax])).sum()
                             for p_short in psd_short])

        # Store variation value
        for i, val in enumerate(diff):
            psd_var[ind+i] = val

        ind = ind+len(diff)

    return psd_var

def find_trigger_value(psd_var, idx, start, sample_rate):
    """ Find the PSD variation value at a particular time

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
    # Find where in the psd variation time series the trigger belongs
    ind = numpy.digitize(time, psd_var.sample_times)
    ind -= 1
    vals = psd_var[ind]
    return vals
