""" PSD Variation """

import numpy

from pycbc.types import TimeSeries, zeros

def calc_psd_variation(strain, psd_short_segment, psd_long_segment,
                       overlap, low_freq, high_freq):
    """Calculates time series of PSD variability

    This function first splits the segment up in to 512 second chunks. It
    then calculates the PSD over this 512 second period as well as in 4
    second chunks throughout each 512 second period. Next the function
    estimates how different the 4 second PSD is to the 512 second PSD and
    produces a timeseries of this variability.

    Parameters
    ----------
    strain : TimeSeries
        Input strain time series to estimate PSDs
    psd_short_segment : {float, 4}
        Duration of the short segments for PSD estimation in seconds.
    psd_long_segment : {float, 512}
        Duration of the long segments for PSD estimation in seconds.
    overlap : {float, 0.5}
        Duration in seconds to use for each sample of the PSD.
    low_freq : {float, 20}
        Minimum frequency to consider the comparison between PSDs.
    high_freq : {float, 512}
        Maximum frequency to consider the comparison between PSDs.

    Returns
    -------
    psd_var : TimeSeries
        Time series of the variability in the PSD estimation
    """

    # Find the times of the long segments
    times_long = numpy.arange(float(strain.start_time),
                              float(strain.end_time), psd_long_segment)

    # Set up the empty time series for the PSD variation estimate
    psd_var = TimeSeries(zeros(int((strain.end_time - strain.start_time) /
                  psd_short_segment)), delta_t=psd_short_segment,
                  copy=False, epoch=strain.start_time)

    ind = 0
    for tlong in times_long:
        # Calculate PSD for long segment and separate the long segment in to
        # Overlapping shorter segments
        if tlong + psd_long_segment <= float(strain.end_time):
            psd_long = strain.time_slice(tlong,
                          tlong + psd_long_segment).psd(overlap)
            times_short = numpy.arange(tlong, tlong + psd_long_segment,
                                       psd_short_segment)
        else:
            psd_long = strain.time_slice(float(strain.end_time) -
                           psd_long_segment,
                           float(strain.end_time)).psd(overlap)
            times_short = numpy.arange(tlong, float(strain.end_time),
                                       psd_short_segment)

        # Caculate the PSD of the shorter segments
        psd_short = []
        for tshort in times_short:
            if tshort + psd_short_segment*2 <= float(strain.end_time):
                pshort = strain.time_slice(tshort,
                             tshort + psd_short_segment*2).psd(overlap)
            else:
                pshort = strain.time_slice(tshort - psd_short_segment,
                             float(strain.end_time)).psd(overlap)
            psd_short.append(pshort)
        # Estimate the range of the PSD to compare
        kmin = int(low_freq / psd_long.delta_f)
        kmax = int(high_freq / psd_long.delta_f)
        # Comapre the PSD of the short segment to the long segment
        #diff = numpy.array([numpy.std((p_short[kmin:kmax] /
        #           psd_long[kmin:kmax])) for p_short in psd_short])
        diff = numpy.array([(p_short[kmin:kmax] / psd_long[kmin:kmax]).sum()
                           for p_short in psd_short])
        diff /= (kmax - kmin)

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
