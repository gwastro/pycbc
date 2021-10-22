# Copyright (C) 2017  Alex Nitz
#
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
from pycbc.types import TimeSeries, complex_same_precision_as
from numpy.random import RandomState
from numpy import log
from scipy.interpolate import interp1d
from scipy.signal import firwin2
import numpy, pycbc.psd
np = numpy

# This constant need to be constant to be able to recover identical results.
BLOCK_SAMPLES = 1638400

def block(seed, sample_rate):
    """ Return block of normal random numbers

    Parameters
    ----------
    seed : {None, int}
        The seed to generate the noise.sd
    sample_rate: float
        Sets the variance of the white noise

    Returns
    --------
    noise : numpy.ndarray
        Array of random numbers
    """
    num = BLOCK_SAMPLES
    rng = RandomState(seed % 2**32)
    variance = sample_rate / 2
    return rng.normal(size=num, scale=variance**0.5)

def normal(start, end, sample_rate=16384, seed=0):
    """ Generate data with a white Gaussian (normal) distribution

    Parameters
    ----------
    start_time : int
        Start time in GPS seconds to generate noise
    end_time : int
        End time in GPS seconds to generate noise
    sample-rate: float
        Sample rate to generate the data at. Keep constant if you want to
        ensure continuity between disjoint time spans.
    seed : {None, int}
        The seed to generate the noise.

    Returns
    --------
    noise : TimeSeries
        A TimeSeries containing gaussian noise
    """
    # This is reproduceable because we used fixed seeds from known values
    block_dur = BLOCK_SAMPLES / sample_rate
    s = int(numpy.floor(start / block_dur))
    e = int(numpy.floor(end / block_dur))

    # The data evenly divides so the last block would be superfluous
    if end % block_dur == 0:
        e -= 1

    sv = RandomState(seed).randint(-2**50, 2**50)
    data = numpy.concatenate([block(i + sv, sample_rate)
                              for i in numpy.arange(s, e + 1, 1)])
    ts = TimeSeries(data, delta_t=1.0 / sample_rate, epoch=(s * block_dur))
    return ts.time_slice(start, end)

def colored_noise(psd, start_time, end_time,
                  seed=0, sample_rate=16384,
                  low_frequency_cutoff=1.0,
                  filter_duration=128):
    """ Create noise from a PSD

    Return noise from the chosen PSD. Note that if unique noise is desired
    a unique seed should be provided.

    Parameters
    ----------
    psd : pycbc.types.FrequencySeries
        PSD to color the noise
    start_time : int
        Start time in GPS seconds to generate noise
    end_time : int
        End time in GPS seconds to generate nosie
    seed : {None, int}
        The seed to generate the noise.
    sample_rate: {16384, float}
        The sample rate of the output data. Keep constant if you want to
        ensure continuity between disjoint time spans.
    low_frequency_cutof : {1.0, float}
        The low frequency cutoff to pass to the PSD generation.
    filter_duration : {128, float}
        The duration in seconds of the coloring filter

    Returns
    --------
    noise : TimeSeries
        A TimeSeries containing gaussian noise colored by the given psd.
    """
    psd = psd.copy()

    flen = int(sample_rate / psd.delta_f) // 2 + 1
    oldlen = len(psd)
    psd.resize(flen)

    # Want to avoid zeroes in PSD.
    max_val = psd.max()
    for i in range(len(psd)):
        if i >= (oldlen-1):
            psd.data[i] = psd[oldlen - 2]
        if psd[i] == 0:
            psd.data[i] = max_val

    fil_len = int(filter_duration * sample_rate)
    wn_dur = int(end_time - start_time) + 2 * filter_duration
    if psd.delta_f >= 1. / (2.*filter_duration):
        # If the PSD is short enough, this method is less memory intensive than
        # resizing and then calling inverse_spectrum_truncation
        psd = pycbc.psd.interpolate(psd, 1.0 / (2. * filter_duration))
        # inverse_spectrum_truncation truncates the inverted PSD. To truncate
        # the non-inverted PSD we give it the inverted PSD to truncate and then
        # invert the output.
        psd = 1. / pycbc.psd.inverse_spectrum_truncation(
                                1./psd,
                                fil_len,
                                low_frequency_cutoff=low_frequency_cutoff,
                                trunc_method='hann')
        psd = psd.astype(complex_same_precision_as(psd))
        # Zero-pad the time-domain PSD to desired length. Zeroes must be added
        # in the middle, so some rolling between a resize is used.
        psd = psd.to_timeseries()
        psd.roll(fil_len)
        psd.resize(int(wn_dur * sample_rate))
        psd.roll(-fil_len)
        # As time series is still mirrored the complex frequency components are
        # 0. But convert to real by using abs as in inverse_spectrum_truncate
        psd = psd.to_frequencyseries()
    else:
        psd = pycbc.psd.interpolate(psd, 1.0 / wn_dur)
        psd = 1. / pycbc.psd.inverse_spectrum_truncation(
                                1./psd,
                                fil_len,
                                low_frequency_cutoff=low_frequency_cutoff,
                                trunc_method='hann')

    kmin = int(low_frequency_cutoff / psd.delta_f)
    psd[:kmin].clear()
    asd = (psd.squared_norm())**0.25
    del psd

    white_noise = normal(start_time - filter_duration,
                         end_time + filter_duration,
                         seed=seed,
                         sample_rate=sample_rate)
    white_noise = white_noise.to_frequencyseries()
    # Here we color. Do not want to duplicate memory here though so use '*='
    white_noise *= asd
    del asd
    colored = white_noise.to_timeseries()
    del white_noise
    return colored.time_slice(start_time, end_time)

def noise_from_string(psd_name, start_time, end_time,
                      seed=0,
                      sample_rate=16384,
                      low_frequency_cutoff=1.0,
                      filter_duration=128):
    """ Create noise from an analytic PSD

    Return noise from the chosen PSD. Note that if unique noise is desired
    a unique seed should be provided.

    Parameters
    ----------
    psd_name : str
        Name of the analytic PSD to use.
    start_time : int
        Start time in GPS seconds to generate noise
    end_time : int
        End time in GPS seconds to generate nosie
    seed : {None, int}
        The seed to generate the noise.
    sample_rate: {16384, float}
        The sample rate of the output data. Keep constant if you want to
        ensure continuity between disjoint time spans.
    low_frequency_cutof : {10.0, float}
        The low frequency cutoff to pass to the PSD generation.
    filter_duration : {128, float}
        The duration in seconds of the coloring filter

    Returns
    --------
    noise : TimeSeries
        A TimeSeries containing gaussian noise colored by the given psd.
    """
    delta_f = 1.0 / filter_duration
    flen = int(sample_rate / delta_f) // 2 + 1
    psd = pycbc.psd.from_string(psd_name, flen, delta_f, low_frequency_cutoff)
    return colored_noise(psd, start_time, end_time,
                         seed=seed,
                         sample_rate=sample_rate,
                         low_frequency_cutoff=low_frequency_cutoff,
                         filter_duration=filter_duration)


def _compute_FIR_coefficients(freqs, amps, nyquist, n_taps):
    """
    Compute the FIR filter coefficients given a filter via frequency series.

    Parameters
    ----------
    freqs : (N,) array
        Frequencies at which the filter is defined.
    amps : (N,) array
        PSD amplitudes of the filter.
    nyquist : float
        Nyquist frequency of the time series to be filtered,
        `nyquist = sampling_frequency / 2`.
    n_taps : int
        Number of "taps", i.e. FIR filter coefficients. A higher value
        leads to a more precise approximation of the desired filter.

    Returns
    -------
    coefficients : (n_taps,) array
        FIR filter coefficients as computed via `scipy.signal.firwin2`.
    """
    # not an expert on the details, but firwin2 requires the frequencies to
    #   go from 0 to 1 (nyquist scale) and the gain to be 0 at the endpoints
    f_scaled = freqs / nyquist  # firwin2 works in units of nyquist
    f_scaled[0] = 0.  # similarly, frequencies need to start at 0
    # brute force the frequency value, floating point math is annoying
    f_scaled[-1] = 1.
    amps[0] = 0.
    amps[-1] = 0.
    assert f_scaled[-1] == 1., repr(f_scaled[-1])  # last element should be 1
    coefficients = firwin2(numtaps=n_taps, freq=f_scaled, gain=amps)
    return coefficients


def variable_noise(
        psd_fn, start_time, end_time, psd_resolution, filter_duration=128.,
        seed=0, sample_rate=16384., time_resolution=-1000, progress=False):
    """ Create noise from a time dependent PSD

    Parameters
    ----------
    psd_fn : callable(array, float) -> FrequencySeries
        Function which provides the time varying PSD. Must accept a call
        signature of `psd_fn(freqs, tau)`, where `freqs` is a numpy array
        of frequencies in [Hz], and `tau` is the normalised time, i.e.
        `current_time / duration` if the start time were at 0. Must return
        a FrequencySeries with sample_frequencies=`freqs`. Ideally, the
        Frequency series covers the interval (0, nyquist), the endpoints
        will be set to 0 due to constraints in the filter construction.
    start_time : int
        Start time in GPS seconds to generate noise
    end_time : int
        End time in GPS seconds to generate nosie
    psd_resolution : int
        Number of points for which the psd_fn in evaluated before being
        passed to the filter generation. Using a large number should be
        generally harmless and without side effects. Recommended to be
        a few time the filter resolution `n_taps`.
    filter_duration : {128, float}
        Length of the filter in [sec]. From this and ´sample_rate´ the
        number of FIR filter coefficients is calculated as
        ´int(filter_duration * sample_rate)´.
    seed : {None, int}
        Seed used to generate the noise.
    sample_rate : {16384, float}
        Sample rate of the generate time series in units of [sec]. Required
        as the filter design depends on the Nyquist frequency.
    time_resolution : {-1000, float}
        The FIR filter is recomputed every `time_resolution` seconds to
        smoothly vary the noise PSD over the length of the time series.
        This directly impacts performance as filter application is much
        cheaper than the recomputing of filter coefficients.
        If the value is larger than the number of samples, a warning is
        issued if `quiet` is not set to `True`. If the value is negative,
        the filter coefficient will be recomputed at ´int(time_resolution)´
        evenly spaced break points.
    progress : {False, bool}
        If `True`, display a simple progress counter.

    Returns
    -------
    noisy_time_series : TimeSeries
        A TimeSeries containing noise with a time variant PSD.

    Notes
    -----
    Currently, the overall noise level may be off by some constant factor.
    Check this and remove the warning before serious use.

    Reference for the FIR filter application:
        https://en.wikipedia.org/wiki/Finite_impulse_response#Definition

    """
    def zero_pad_left(arr, n):
        """
        Zero pad an array to a given size. Zeros are added on the left,
        i.e. at the beginning.

        Parameters
        ----------
        arr : (N,) array
            Some array of floats.
        n : int
            Size to which the array is to be padded.

        Returns
        -------
        padded : (n,) array
            Left padded array.

        """
        if arr.size >= n:
            return arr
        padded = np.concatenate([np.zeros(n - arr.size), arr])
        return padded

    assert time_resolution != 0

    # first: compute some essential constants
    duration = end_time - start_time
    n_taps = int(sample_rate * filter_duration)
    nyquist = sample_rate / 2.  # [Hz]
    n_samples = int(duration * sample_rate)
    #   use some very small number here, f[0] will be set to 0 later anyways
    freqs = np.geomspace(1e-100, nyquist, psd_resolution)

    # generate white noise, which will then be filtered to
    #   obtain the desired PSD.
    white_time_series = normal(start_time, end_time, sample_rate, seed)
    RNG = np.random.default_rng(seed)
    white_time_series = RNG.standard_normal(size=n_samples)
    assert white_time_series.size == n_samples

    # prepare TimeSeries for filter output
    noisy_time_series = white_time_series.copy()

    # precompute the indices for re-computation of the FIR coefficients
    # first: turn the time between recomputations into a number of samples
    #   note to treat the negative time_resolution case seperately
    if time_resolution < 0:
        idx_step = int(n_samples / abs(int(time_resolution)))
    else:
        idx_step = int(duration / time_resolution)
    assert idx_step > 0
    break_indices = range(0, n_samples, idx_step)

    # This loop is a prime candidate for optimisation,
    #   e.g. pre-computing FIR coefficients and then dropping into cython
    #   for the convolution.
    for idx in range(n_samples):
        # ^ no enumerate since the samples are only used in chunks
        tau = float(idx) / n_samples
        if idx in break_indices:
            if progress:
                print('Recompute at sample idx {} of {} ({:.1%})'
                      ''.format(idx, n_samples, idx/n_samples), flush=True)
            # re-compute the coefficients
            amps = psd_fn(freqs, tau)
            assert amps.size == psd_resolution
            coeffs = _compute_FIR_coefficients(freqs, amps, nyquist, n_taps)
            # as a quirk of definition, the coefficients are time reversed
            coeffs = coeffs[::-1]
            assert coeffs.size == n_taps
        # select the chunk of the time series which is used by the filter
        #   i.e. the n_taps elements before the current one
        chunk = zero_pad_left(white_time_series[max(idx-n_taps, 0)+1:idx+1],
                              n_taps)
        #   max-term : avoid selecting before element 0
        #   idx+1 : the current sample needs to be included
        # convolution
        noisy_time_series[idx] = np.sum(chunk * coeffs)

    return noisy_time_series


# example for creating the input `psd_fn`

def generate_psd_fn_linear_in_log_transition(psd_0, psd_1,
                                             psd_fill_value=1e-100,
                                             frequency_gap=1e-3):
    """
    Generate a psd_fn according to the requirements of ´variable_noise(...)´

    psd_0 and psd_1 should start & end at the same frequencies, otherwise
    their values are transitioned to/from the psd_fill_value. The
    transition between psd_0 and psd_1 is linear in log-log space, i.e.
    psd(f) = exp(tau * psd_0(f) + (1 - tau) * psd_1(f)).

    Parameters
    ----------
    psd_0 : FrequencySeries
        PSD at the start of the noise segment.
    psd_1 : FrequencySeries
        PSD at the end of the noise segment.
    psd_fill_value : {1e-100, float}
        Filler value outside the defined frequency range.
    frequency_gap : {1e-3, float}
        Interval in [Hz] over which the PSDs are transitioned to the
        filler value at the boundaries.

    Returns
    -------
    psd_fn : FrequencySeries
        psd_fn(freqs, tau) as consumed by ´variable_noise(...)´.

    """
    freqs_0 = psd_0.sample_frequencies
    gains_0 = np.array(psd_0)
    freqs_1 = psd_1.sample_frequencies
    gains_1 = np.array(psd_1)
    # add the transition gaps
    freqs_0 = np.concatenate([[freqs_0[0] - frequency_gap], freqs_0,
                              [freqs_0[-1] + frequency_gap]])
    gains_0 = np.concatenate([[psd_fill_value], gains_0, [psd_fill_value]])
    freqs_1 = np.concatenate([[freqs_1[0] - frequency_gap], freqs_1,
                              [freqs_1[-1] + frequency_gap]])
    gains_1 = np.concatenate([[psd_fill_value], gains_1, [psd_fill_value]])
    # insert a transition from f_low and
    interpolant_0 = interp1d(x=log(psd_0.sample_frequencies),
                             y=log(np.array(psd_0)), kind='linear',
                             bounds_error=False,
                             fill_value=log(psd_fill_value))
    interpolant_1 = interp1d(x=log(psd_1.sample_frequencies),
                             y=log(np.array(psd_1)), kind='linear',
                             bounds_error=False,
                             fill_value=log(psd_fill_value))

    def psd_fn(freqs, tau):
        assert tau >= 0., repr(tau)
        assert tau <= 1., repr(tau)
        return np.exp(tau * interpolant_0(log(freqs))
                      + (1. - tau) * interpolant_1(log(freqs)))

    return psd_fn
