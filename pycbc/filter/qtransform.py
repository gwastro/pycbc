# Copyright (C) 2017  Hunter A. Gabbard
# Most of this is a port of Duncan Macleod's GWPY qtransform.py script
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

from math import pi, ceil, log, exp
import numpy as np
from pycbc.types.timeseries import FrequencySeries, TimeSeries
from ../strain  import next_power_of_2
import os, sys
from pycbc.filter import highpass_fir, matched_filter
from pycbc.waveform import get_fd_waveform
from pycbc.psd import welch, interpolate
from pycbc.fft import ifft
from scipy.interpolate import (interp2d, InterpolatedUnivariateSpline)
from numpy import fft as npfft

from matplotlib import use
use('Agg')
from matplotlib import pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.pyplot import specgram

__author__ = 'Hunter Gabbard <hunter.gabbard@ligo.org>'
__credits__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

def padding(window_size, dur, f0, Q):
    """The `(left, right)` padding required for the IFFT

    Parameters
    ----------
    window_size: int
        Size of window
    dur: int
        Duration of timeseries in seconds
    f0: int
        Central frequency
    Q: int
        q value

    Returns
    -------
    tuple
       Number of values padded to the edges of each axis. 
 
    """

    pad = n_tiles(dur,f0,Q) - window_size
    return (int((pad - 1)/2.), int((pad + 1)/2.))

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

def _get_indices(dur):
    """ Windows indices for fft
    
    Paramters
    ---------
    dur: int
        Duration of timeseries in seconds

    Returns
    -------
    numpy.ndarray
        Window indices for fft using total duration of segment

    """
    half = int((int(dur) - 1.) / 2.)
    return np.arange(-half, half + 1)

def get_window(dur, indices, f0, qprime, Q, sampling):
    """Generate the bi-square window for this row
 
    Paramters
    ---------
    dur: int
        Duration of timeseries in seconds
    f0: int
        Central frequency
    qprime: int
        Normalized Q `(q/sqrt(11))
    Q: int
        q value
    sampling: int
        sampling frequency of timeseries

    Returns
    -------
    window : `numpy.ndarray`
    """
    # real frequencies
    wfrequencies = indices / dur

    # dimensionless frequencies
    xfrequencies = wfrequencies * qprime / f0

    # normalize and generate bi-square window
    norm = n_tiles(dur,f0,Q) / (dur * sampling) * (
        315 * qprime / (128 * f0)) ** (1/2.)
    return (1 - xfrequencies ** 2) ** 2 * norm

def n_tiles(dur,f0,Q):
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

    tcum_mismatch = dur * 2 * pi * f0 / Q
    return next_power_of_2(tcum_mismatch / deltam())

def deltam():
    """Fractional mismatch between neighbouring tiles

    Parameters
    ----------
    None

    Returns
    -------
    :type: `float`
    """
    mismatch = 0.2
    return 2 * (mismatch / 3.) ** (1/2.)

def main():
    return None

if __name__ == '__main__':
    main()
