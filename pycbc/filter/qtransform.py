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


    Parameters
    ----------
    fseries: 'pycbc FrequencySeries'
        frequency-series data set
    qrange:
        upper and lower bounds of q range
    frange:
        upper and lower bounds of frequency range
    sampling:
        sampling rate of time-series
    mismatch:
        percentage of desired fractional mismatch

    Returns
    -------
    qplane_tile_dict: 'dict'
        dictionary containing Q-tile tuples for a set of Q-planes 
    frange: 'list'
        upper and lower bounds on frequency range   
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
    fseries: 'pycbc FrequencySeries'
        frequency-series data set
    Q:
        q value
    f0:
        central frequency

    Returns
    -------
    norm_energy: '~pycbc.types.aligned.ArrayWithAligned'
        A 'TimeSeries' of the normalized energy from the Q-transform of 
        this tile against the data.
    cenergy: '~pycbc.types.aligned.ArrayWithAligned'
        A 'TimeSeries' of the complex energy from the Q-transform of 
        this tile against the data.
    """

    # q-transform data for each (Q, frequency) tile

    # initialize parameters
    qprime = Q / 11**(1/2.) # ... self.qprime
    dur = fseries.to_timeseries().duration

    # window fft
    window_size = 2 * int(f0 / qprime * dur) + 1

    # get indices
    indices = _get_indices(window_size)

    # apply window to fft
    windowed = fseries[get_data_indices(dur, f0, indices)] * get_window(window_size, f0, qprime)

    # Choice of output sampling rate
    output_sampling = fseries.delta_f # Can lower this to highest bandwidth
    output_samples = dur * output_sampling

    # pad data, move negative frequencies to the end, and IFFT
    padded = np.pad(windowed, padding(window_size, output_samples), mode='constant')
    wenergy = npfft.ifftshift(padded)

    # return a 'TimeSeries'
    wenergy = FrequencySeries(wenergy, delta_f=1./dur)
    tdenergy = TimeSeries(zeros(output_samples, dtype=np.complex128),
                            delta_t=1./output_sampling)
    ifft(wenergy, tdenergy)
    cenergy = TimeSeries(tdenergy,
                         delta_t=tdenergy.delta_t, copy=False)
    if normalized:
        energy = type(cenergy)(
            cenergy.real() ** 2. + cenergy.imag() ** 2.,
            delta_t=1, copy=False)
        medianenergy = energy.numpy().median()
        result = energy / medianenergy
    else:
        result = cenergy

    return result


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

def get_data_indices(dur, f0, indices):
    """Returns the index array of interesting frequencies for this row

    Parameters
    ----------
    dur: int
        Duration of timeseries in seconds
    f0: int
        Central frequency
    indices: numpy.ndarray
        window indices for fft 

    Returns
    -------
    numpy.ndarray
        Returns index array of interesting frequencies for this row

    """
    return np.round(indices + 1 +
                       f0 * dur).astype(int)

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

def get_window(size, f0, qprime):
    """Generate the bi-square window for this row
 
    Parameters
    ---------
    size: int
        size of window
    f0: int
        Central frequency
    qprime: int
        Normalized Q (q/sqrt(11))

    Returns
    -------
    window : numpy.ndarray
    """
    # dimensionless frequencies
    xfrequencies = np.linspace(-1., 1., size)

    # ported from https://github.com/gwpy/gwpy/blob/master/gwpy/signal/qtransform.py
    norm = np.sqrt(315. * qprime / (128. * f0))
    return (1 - xfrequencies ** 2) ** 2 * norm

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
