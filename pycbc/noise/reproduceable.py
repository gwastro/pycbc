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
from six.moves import range
import numpy, pycbc.psd
from pycbc.types import TimeSeries, FrequencySeries, complex_same_precision_as
from numpy.random import RandomState

# These need to be constant to be able to recover identical results.
# The hope is that nobody needs a higher resolution
SAMPLE_RATE = 16384
BLOCK_SIZE = 100
FILTER_LENGTH = 128

def block(seed):
    """ Return block of normal random numbers

    Parameters
    ----------
    seed : {None, int}
        The seed to generate the noise.sd

    Returns
    --------
    noise : numpy.ndarray
        Array of random numbers
    """
    num = SAMPLE_RATE * BLOCK_SIZE
    rng = RandomState(seed % 2**32)
    variance = SAMPLE_RATE / 2
    return rng.normal(size=num, scale=variance**0.5)

def normal(start, end, seed=0):
    """ Generate data with a white Gaussian (normal) distribution

    Parameters
    ----------
    start_time : int
        Start time in GPS seconds to generate noise
    end_time : int
        End time in GPS seconds to generate nosie
    seed : {None, int}
        The seed to generate the noise.

    Returns
    --------
    noise : TimeSeries
        A TimeSeries containing gaussian noise
    """
    # This is reproduceable because we used fixed seeds from known values
    s = int(start / BLOCK_SIZE)
    e = int(end / BLOCK_SIZE)

    # The data evenly divides so the last block would be superfluous
    if end % BLOCK_SIZE == 0:
        e -= 1

    sv = RandomState(seed).randint(-2**50, 2**50)
    data = numpy.concatenate([block(i + sv) for i in numpy.arange(s, e + 1, 1)])
    ts = TimeSeries(data, delta_t=1.0 / SAMPLE_RATE, epoch=start)
    return ts.time_slice(start, end)

def colored_noise(psd, start_time, end_time, seed=0, low_frequency_cutoff=1.0):
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
    low_frequency_cutof : {1.0, float}
        The low frequency cutoff to pass to the PSD generation.

    Returns
    --------
    noise : TimeSeries
        A TimeSeries containing gaussian noise colored by the given psd.
    """
    psd = psd.copy()

    flen = int(SAMPLE_RATE / psd.delta_f) / 2 + 1
    oldlen = len(psd)
    psd.resize(flen)

    # Want to avoid zeroes in PSD.
    max_val = psd.max()
    for i in range(len(psd)):
        if i >= (oldlen-1):
            psd.data[i] = psd[oldlen - 2]
        if psd[i] == 0:
            psd.data[i] = max_val

    wn_dur = int(end_time - start_time) + 2*FILTER_LENGTH
    if psd.delta_f >= 1. / (2.*FILTER_LENGTH):
        # If the PSD is short enough, this method is less memory intensive than
        # resizing and then calling inverse_spectrum_truncation
        psd = pycbc.psd.interpolate(psd, 1.0 / (2.*FILTER_LENGTH))
        # inverse_spectrum_truncation truncates the inverted PSD. To truncate
        # the non-inverted PSD we give it the inverted PSD to truncate and then
        # invert the output.
        psd = 1. / pycbc.psd.inverse_spectrum_truncation(1./psd,
                                FILTER_LENGTH * SAMPLE_RATE,
                                low_frequency_cutoff=low_frequency_cutoff,
                                trunc_method='hann')
        psd = psd.astype(complex_same_precision_as(psd))
        # Zero-pad the time-domain PSD to desired length. Zeroes must be added
        # in the middle, so some rolling between a resize is used.
        psd = psd.to_timeseries()
        psd.roll(SAMPLE_RATE * FILTER_LENGTH)
        psd.resize(wn_dur * SAMPLE_RATE)
        psd.roll(-SAMPLE_RATE * FILTER_LENGTH)
        # As time series is still mirrored the complex frequency components are
        # 0. But convert to real by using abs as in inverse_spectrum_truncate
        psd = psd.to_frequencyseries()
    else:
        psd = pycbc.psd.interpolate(psd, 1.0 / wn_dur)
        psd = 1. / pycbc.psd.inverse_spectrum_truncation(1./psd,
                                FILTER_LENGTH * SAMPLE_RATE,
                                low_frequency_cutoff=low_frequency_cutoff,
                                trunc_method='hann')

    kmin = int(low_frequency_cutoff / psd.delta_f)
    psd[:kmin].clear()
    asd = (psd.real())**0.5
    del psd

    white_noise = normal(start_time - FILTER_LENGTH, end_time + FILTER_LENGTH,
                         seed=seed)
    white_noise = white_noise.to_frequencyseries()
    # Here we color. Do not want to duplicate memory here though so use '*='
    white_noise *= asd
    del asd
    colored = white_noise.to_timeseries()
    del white_noise
    return colored.time_slice(start_time, end_time)

def noise_from_string(psd_name, start_time, end_time, seed=0, low_frequency_cutoff=1.0):
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
    low_frequency_cutof : {10.0, float}
        The low frequency cutoff to pass to the PSD generation.

    Returns
    --------
    noise : TimeSeries
        A TimeSeries containing gaussian noise colored by the given psd.
    """
    delta_f = 1.0 / FILTER_LENGTH
    flen = int(SAMPLE_RATE / delta_f) / 2 + 1
    psd = pycbc.psd.from_string(psd_name, flen, delta_f, low_frequency_cutoff)
    return colored_noise(psd, start_time, end_time,
                         seed=seed,
                         low_frequency_cutoff=low_frequency_cutoff)
