# Copyright (C) 2017  Hunter A. Gabbard, Andrew Lundgren,
#                     Duncan Macleod, Alex Nitz
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

"""
This module retrives a timeseries and then calculates
the q-transform of that time series
"""

from six.moves import range
import numpy
from numpy import ceil, log, exp
from pycbc.types.timeseries import FrequencySeries, TimeSeries
from pycbc.fft import ifft
from pycbc.types import zeros

def qplane(qplane_tile_dict, fseries, return_complex=False):
    """Performs q-transform on each tile for each q-plane and selects
       tile with the maximum energy. Q-transform can then
       be interpolated to a desired frequency and time resolution.

    Parameters
    ----------
    qplane_tile_dict:
        Dictionary containing a list of q-tile tupples for each q-plane
    fseries: 'pycbc FrequencySeries'
        frequency-series data set
    return_complex: {False, bool}
        Return the raw complex series instead of the normalized power.

    Returns
    -------
    q : float
        The q of the maximum q plane
    times : numpy.ndarray
        The time that the qtransform is sampled.
    freqs : numpy.ndarray
        The frequencies that the qtransform is samled.
    qplane : numpy.ndarray (2d)
        The two dimensional interpolated qtransform of this time series.
    """
    # store q-transforms for each q in a dict
    qplanes = {}
    max_energy, max_key = None, None
    for i, q in enumerate(qplane_tile_dict):
        energies = []
        for f0 in qplane_tile_dict[q]:
            energy = qseries(fseries, q, f0, return_complex=return_complex)
            menergy = abs(energy).max()
            energies.append(energy)

            if i == 0 or menergy > max_energy:
                max_energy = menergy
                max_key = q

        qplanes[q] = energies

    # record q-transform output for peak q
    plane = qplanes[max_key]
    frequencies = qplane_tile_dict[max_key]
    times = plane[0].sample_times.numpy()
    plane = numpy.array([v.numpy() for v in plane])
    return max_key, times, frequencies, numpy.array(plane)

def qtiling(fseries, qrange, frange, mismatch=0.2):
    """Iterable constructor of QTile tuples

    Parameters
    ----------
    fseries: 'pycbc FrequencySeries'
        frequency-series data set
    qrange:
        upper and lower bounds of q range
    frange:
        upper and lower bounds of frequency range
    mismatch:
        percentage of desired fractional mismatch

    Returns
    -------
    qplane_tile_dict: 'dict'
        dictionary containing Q-tile tuples for a set of Q-planes
    """
    qplane_tile_dict = {}
    qs = list(_iter_qs(qrange, deltam_f(mismatch)))
    for q in qs:
        qtilefreq = _iter_frequencies(q, frange, mismatch, fseries.duration)
        qplane_tile_dict[q] = numpy.array(list(qtilefreq))

    return qplane_tile_dict

def deltam_f(mismatch):
    """Fractional mismatch between neighbouring tiles

    Parameters
    ----------
    mismatch: 'float'
        percentage of desired fractional mismatch

    Returns
    -------
    :type: 'float'
    """
    return 2 * (mismatch / 3.) ** (1/2.)

def _iter_qs(qrange, deltam):
    """Iterate over the Q values

    Parameters
    ----------
    qrange:
        upper and lower bounds of q range
    deltam:
        Fractional mismatch between neighbouring tiles

    Returns
    -------
    Q-value:
        Q value for Q-tile
    """

    # work out how many Qs we need
    cumum = log(float(qrange[1]) / qrange[0]) / 2**(1/2.)
    nplanes = int(max(ceil(cumum / deltam), 1))
    dq = cumum / nplanes
    for i in range(nplanes):
        yield qrange[0] * exp(2**(1/2.) * dq * (i + .5))
    raise StopIteration()

def _iter_frequencies(q, frange, mismatch, dur):
    """Iterate over the frequencies of this 'QPlane'

    Parameters
    ----------
    q:
        q value
    frange: 'list'
        upper and lower bounds of frequency range
    mismatch:
        percentage of desired fractional mismatch
    dur:
        duration of timeseries in seconds

    Returns
    -------
    frequencies:
        Q-Tile frequency
    """
    # work out how many frequencies we need
    minf, maxf = frange
    fcum_mismatch = log(float(maxf) / minf) * (2 + q**2)**(1/2.) / 2.
    nfreq = int(max(1, ceil(fcum_mismatch / deltam_f(mismatch))))
    fstep = fcum_mismatch / nfreq
    fstepmin = 1. / dur
    # for each frequency, yield a QTile
    for i in range(nfreq):
        yield (float(minf) *
               exp(2 / (2 + q**2)**(1/2.) * (i + .5) * fstep) //
               fstepmin * fstepmin)
    raise StopIteration()

def qseries(fseries, Q, f0, return_complex=False):
    """Calculate the energy 'TimeSeries' for the given fseries

    Parameters
    ----------
    fseries: 'pycbc FrequencySeries'
        frequency-series data set
    Q:
        q value
    f0:
        central frequency
    return_complex: {False, bool}
        Return the raw complex series instead of the normalized power.

    Returns
    -------
    energy: '~pycbc.types.TimeSeries'
        A 'TimeSeries' of the normalized energy from the Q-transform of
        this tile against the data.
    """
    # normalize and generate bi-square window
    qprime = Q / 11**(1/2.)
    norm = numpy.sqrt(315. * qprime / (128. * f0))
    window_size = 2 * int(f0 / qprime * fseries.duration) + 1
    xfrequencies = numpy.linspace(-1., 1., window_size)

    start = int((f0 - (f0 / qprime)) * fseries.duration)
    end = int(start + window_size)
    center = (start + end) / 2

    windowed = fseries[start:end] * (1 - xfrequencies ** 2) ** 2 * norm

    tlen = (len(fseries)-1) * 2
    windowed.resize(tlen)
    windowed.roll(-center)

    # calculate the time series for this q -value
    windowed = FrequencySeries(windowed, delta_f=fseries.delta_f,
                            epoch=fseries.start_time)
    ctseries = TimeSeries(zeros(tlen, dtype=numpy.complex128),
                            delta_t=fseries.delta_t)
    ifft(windowed, ctseries)

    if return_complex:
        return ctseries
    else:
        energy = ctseries.squared_norm()
        medianenergy = numpy.median(energy.numpy())
        return  energy / float(medianenergy)
