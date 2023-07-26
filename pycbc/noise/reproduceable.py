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
import numpy, pycbc.psd
from pycbc.types import TimeSeries, complex_same_precision_as
from numpy.random import RandomState

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
                  filter_duration=128,
                  scale=1.0):
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
    white_noise *= asd*scale
    del asd
    colored = white_noise.to_timeseries(delta_t=1.0/sample_rate)
    del white_noise
    return colored.time_slice(start_time, end_time)

def noise_from_string(psd_name, start_time, end_time,
                      seed=0,
                      sample_rate=16384,
                      low_frequency_cutoff=1.0,
                      filter_duration=128,
                      scale=1.0):
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
                         filter_duration=filter_duration,
                         scale=scale)
