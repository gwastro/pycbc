# Copyright (C) 2017  Hunter A. Gabbard, Andrew Lundgren, Duncan Macleod
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
This module retrives a timeseries and then calculates the q-transform of that time series

Example
-------
    $ python q-transform.py -s 4096 -u test -o /Users/pycbc_qtransform

"""

from math import pi
import numpy as np
from pycbc.strain  import next_power_of_2
from pycbc.types.timeseries import FrequencySeries, TimeSeries
from numpy import fft as npfft

def qtiling(fseries, qrange, frange, sampling, mismatch):
    """Iterable constructor of QTile tuples

    Parameters
    ----------
    fseries: 'LIGO gwf frame file'
        raw frequency-series data set
    qrange:
        lower and upper bound on q values    
    frange:
        lower and upper bound on frequency values
    sampling:
        sampling rate of channel
    mismatch:
        fractional mismatch percentage

    Returns
    -------
    qplane_tile_dict: 'dict'
        dictionary of tile values for a given Q-range and frequency range
    frange: 'list'
        lower and upper bound on freqeuncy range
    """

    deltam = deltam_f(mismatch)
    qrange = (float(qrange[0]), float(qrange[1]))
    frange = [float(frange[0]), float(frange[1])]
    dur = fseries.to_timeseries().duration
    qplane_tile_dict = {}

    qs = list(_iter_qs(qrange, deltam))
    if frange[0] == 0:  # set non-zero lower frequency
        frange[0] = 50 * max(qs) / (2 * pi * dur)
    if np.isinf(frange[1]):  # set non-infinite upper frequency
        frange[1] = sampling / 2 / (1 + 11**(1/2.) / min(qs))

    #lets now define the whole tiling (e.g. choosing all tiling in planes)
    for q in qs:
        qtilefreq = np.array(list(_iter_frequencies(q, frange, mismatch, dur)))
        qlst = np.empty(len(qtilefreq), dtype=float)
        qlst.fill(q)
        qtiles_array = np.vstack((qtilefreq,qlst)).T
        qplane_tiles_list = list(map(tuple,qtiles_array))
        qplane_tile_dict[q] = qplane_tiles_list

    return qplane_tile_dict, frange


def deltam_f(mismatch):
    """Fractional mismatch between neighbouring tiles

    Parameters
    ----------
    mismatch: 'float'
        fractional mismatch percentage

    Returns
    -------
    FracMismatch: 'float'
    """
    return 2 * (mismatch / 3.) ** (1/2.)

def _iter_qs(qrange, deltam):
    """Iterate over the Q values

    Parameters
    ----------
    qrange: 'list'
        lower and upper bound on q values
    deltam: 'float'
        Fractional mismatch between neighbouring tiles

    Returns
    -------
    Q's: 'list'
        A list of Q values for a given Q-plane
    """

    # work out how many Qs we need
    cumum = log(qrange[1] / qrange[0]) / 2**(1/2.)
    nplanes = int(max(ceil(cumum / deltam), 1))
    dq = cumum / nplanes
    for i in xrange(nplanes):
        yield qrange[0] * exp(2**(1/2.) * dq * (i + .5))
    raise StopIteration()

def _iter_frequencies(q, frange, mismatch, dur):
    """Iterate over the frequencies of this 'QPlane'

    Parameters
    ----------
    q: 
        q value
    frange:
        range of frequencies to iterate over
    mismatch:
        fractional mismatch percentage
    dur:
        duration of analysis period in seconds

    Returns
    -------
    frequencies: 'list'
        A list of frequency values for each QTile
    """
    # work out how many frequencies we need
    minf, maxf = frange
    fcum_mismatch = log(maxf / minf) * (2 + q**2)**(1/2.) / 2.
    nfreq = int(max(1, ceil(fcum_mismatch / deltam_f(mismatch))))
    fstep = fcum_mismatch / nfreq
    fstepmin = 1. / dur
    # for each frequency, yield a QTile
    for i in xrange(nfreq):
        yield (minf *
               exp(2 / (2 + q**2)**(1/2.) * (i + .5) * fstep) //
               fstepmin * fstepmin)
    raise StopIteration()

def qtransform(fseries, Q, f0, sampling):
    """Calculate the energy 'TimeSeries' for the given fseries

    Parameters
    ----------
    fseries: 'LIGO gwf frame file'
        raw frequency-series data set
    Q:
        Q value
    f0:
        central frequency
    sampling :
        sampling rate of channel

    Returns
    -------
    cenergy: '~pycbc.types.aligned.ArrayWithAligned'
        A 'TimeSeries' of the complex energy from the Q-transform of 
        this tile against the data.
    norm_energy: '~pycbc.types.aligned.ArrayWithAligned'
        A 'TimeSeries' of the real energy from the Q-transform of this tile against the data.
    """

    # q-transform data for each (Q, frequency) tile

    # initialize parameters
    qprime = Q / 11**(1/2.)
    dur = fseries.to_timeseries().duration

    # window fft
    window_size = 2 * int(f0 / qprime * dur) + 1

    # get start and end indices
    start = (f0 - (f0 / qprime)) * dur
    end = start + window_size

    # apply window to fft
    # normalize and generate bi-square window
    norm = np.sqrt(315. * qprime / (128. * f0))
    windowed = fseries[start:end] * (bisquare(window_size) * norm)

    # choice of output sampling rate
    output_sampling = sampling # Can lower this to highest bandwidth
    output_samples = dur * output_sampling

    # pad data, move negative frequencies to the end, and IFFT 
    padded = np.pad(windowed, padding(window_size, output_samples), mode='constant')
    wenergy = npfft.ifftshift(padded)

    # return a `TimeSeries`
    wenergy = FrequencySeries(wenergy, delta_f=1./dur)
    tdenergy = TimeSeries(zeros(output_samples, dtype=np.complex128),
                            delta_t=1./sampling)
    ifft(wenergy, tdenergy)
    cenergy = TimeSeries(tdenergy,
                         delta_t=tdenergy.delta_t, copy=False) # Normally delta_t is dur/tdenergy.size ... must figure out better way of doing this
    energy = type(cenergy)(
        cenergy.real() ** 2. + cenergy.imag() ** 2.,
        delta_t=1, copy=False)
    medianenergy = np.median(energy)
    norm_energy = energy / medianenergy

    return norm_energy, cenergy

def padding(window_size, desired_size):
    """The (left, right) padding required for the IFFT

    Parameters
    ----------
    window_size: int
        Size of window
    desired_size: int
        Desired size of window

    Returns
    -------
    tuple
       Number of values padded to the edges of each axis. 
 
    """
    pad = desired_size - window_size
    return (int(pad/2.), int((pad + 1)/2.))


def _get_indices(window_size):
    """ Windows indices for fft
    
    Parameters
    ---------
    window_size: int
        size of window

    Returns
    -------
    numpy.ndarray
        Window indices for fft using total duration of segment

    """
    half = int((windowsize - 1) / 2.)
    return np.arange(-half, half + 1)

def bisquare(size):
    """Generate the bi-square window for this row
 
    Parameters
    ---------
    size: int
        size of window

    Returns
    -------
    window : numpy.ndarray
    """
    # dimensionless frequencies
    xfrequencies = np.linspace(-1., 1., size)

    # ported from https://github.com/gwpy/gwpy/blob/master/gwpy/signal/qtransform.py
    return (1 - xfrequencies ** 2) ** 2

def n_tiles(dur, f0, Q):
    """The number of tiles in this row 
    
    Parameters
    ----------
    dur: int
        Duration of timeseries in seconds
    f0: int
        Central frequency
    Q: int
        q value

    Returns
    -------
    :type: 'int'
    """

    # ported from https://github.com/gwpy/gwpy/blob/master/gwpy/signal/qtransform.py
    tcum_mismatch = dur * 2 * pi * f0 / Q  
    return next_power_of_2(tcum_mismatch / deltam())

def deltam():
    """Fractional mismatch between neighbouring tiles

    Parameters
    ----------
    None

    Returns
    -------
    :type: 'float'
    """
    mismatch = 0.2
    return 2 * (mismatch / 3.) ** (1/2.)

