# Copyright (C) 2012  Alex Nitz
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
"""This module contains functions to generate gaussian noise colored with a
noise spectrum.
"""

from pycbc.types import (
    zeros, complex_same_precision_as, FrequencySeries
)
import numpy.random
import pycbc.psd

def frequency_noise_from_psd(psd, seed=None):
    """ Create noise with a given psd.

    Return noise coloured with the given psd. The returned noise
    FrequencySeries has the same length and frequency step as the given psd.
    Note that if unique noise is desired a unique seed should be provided.

    Parameters
    ----------
    psd : FrequencySeries
        The noise weighting to color the noise.
    seed : {0, int} or None
        The seed to generate the noise. If None specified,
        the seed will not be reset.

    Returns
    --------
    noise : FrequencySeriesSeries
        A FrequencySeries containing gaussian noise colored by the given psd.
    """
    sigma = 0.5 * (psd / psd.delta_f) ** (0.5)
    if seed is not None:
        numpy.random.seed(seed)
    sigma = sigma.numpy()
    dtype = complex_same_precision_as(psd)

    not_zero = (sigma != 0)

    sigma_red = sigma[not_zero]
    noise_re = numpy.random.normal(0, sigma_red)
    noise_co = numpy.random.normal(0, sigma_red)
    noise_red = noise_re + 1j * noise_co

    noise = numpy.zeros(len(sigma), dtype=dtype)
    noise[not_zero] = noise_red

    return FrequencySeries(noise,
                           delta_f=psd.delta_f,
                           dtype=dtype)

def noise_from_psd(length, delta_t, psd, seed=None):
    """ Create noise with a given psd.
    Return noise with a given psd. Note that if unique noise is desired
    a unique seed should be provided.
    Parameters
    ----------
    length : int
        The length of noise to generate in samples.
    delta_t : float
        The time step of the noise.
    psd : FrequencySeries
        The noise weighting to color the noise.
    seed : {0, int}
        The seed to generate the noise.
    Returns
    --------
    noise : TimeSeries
        A TimeSeries containing gaussian noise colored by the given psd.
    """
    delta_f = 1.0 / (length * delta_t)
    flen = int(length / 2) + 1

    resampled_psd = FrequencySeries(zeros(flen), delta_f=delta_f)

    # Resample the given PSD to the frequencies we need
    for i in range(flen):
        f = i * delta_f
        if f <= psd.sample_frequencies[-1]:
            resampled_psd[i] = psd.at_frequency(f)

    # Generate frequency-domain noise
    noise_freq = frequency_noise_from_psd(resampled_psd, seed=seed)

    # Convert to time series
    # The to_timeseries() method will create a time series of the correct
    # length and delta_t
    noise_ts = noise_freq.to_timeseries()

    # The length may be off by one, so we resize if needed
    if len(noise_ts) != length:
        noise_ts.resize(length)

    return noise_ts

def noise_from_string(psd_name, length, delta_t, seed=None, low_frequency_cutoff=10.0):
    """ Create noise from an analytic PSD

    Return noise from the chosen PSD. Note that if unique noise is desired
    a unique seed should be provided.

    Parameters
    ----------
    psd_name : str
        Name of the analytic PSD to use.
    low_fr
    length : int
        The length of noise to generate in samples.
    delta_t : float
        The time step of the noise.
    seed : {None, int}
        The seed to generate the noise.
    low_frequency_cutof : {10.0, float}
        The low frequency cutoff to pass to the PSD generation.

    Returns
    --------
    noise : TimeSeries
        A TimeSeries containing gaussian noise colored by the given psd.
    """
    # We just need enough resolution to resolve lines
    delta_f = 1.0 / 8
    flen = int(.5 / delta_t / delta_f) + 1
    psd = pycbc.psd.from_string(psd_name, flen, delta_f, low_frequency_cutoff)
    return noise_from_psd(int(length), delta_t, psd, seed=seed)
