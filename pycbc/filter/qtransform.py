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
"""

from math import pi, ceil, log, exp
import numpy as np
from pycbc.types.timeseries import FrequencySeries, TimeSeries
import os, sys
from pycbc.frame import read_frame
from pycbc.filter import highpass_fir, matched_filter
from pycbc.waveform import get_fd_waveform
from pycbc.psd import welch, interpolate
from pycbc.fft import ifft
import urllib
import datetime
from scipy.interpolate import (interp2d, InterpolatedUnivariateSpline)
from numpy import fft as npfft
import argparse
import datetime

from matplotlib import use
use('Agg')
from matplotlib import pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.pyplot import specgram

__author__ = 'Hunter Gabbard <hunter.gabbard@ligo.org>'
__credits__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

def get_window(dur, indices, f0, qprime, Q, sampling):
    """Generate the bi-square window for this row
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

    :type: 'int'
    """

    tcum_mismatch = dur * 2 * pi * f0 / Q
    return next_power_of_two(tcum_mismatch / deltam())

def next_power_of_two(x):
    """Return the smallest power of two greater than or equal to `x`
    """
    return 2**(ceil(log(x, 2)))

def deltam():
    """Fractional mismatch between neighbouring tiles
    :type: `float`
    """
    mismatch = 0.2
    return 2 * (mismatch / 3.) ** (1/2.)

def main():
    #Get Current time
    cur_time = datetime.datetime.now()

    #construct the argument parse and parse the arguments
    ap = argparse.ArgumentParser()
    ap.add_argument("-u", "--usertag", required=False, default=cur_time,
        help="label for given run")
    ap.add_argument("-o", "--output-dir", required=False,
        help="path to output directory")
    ap.add_argument("-n", "--normalize", required=False, default=True,
        help="normalize the energy of the output")
    ap.add_argument("-s", "--samp-freq", required=True, type=float,
        help="Sampling frequency of channel")

    args = ap.parse_args()


    #Initialize parameters
    out_dir = args.output_dir
    now = args.usertag
    os.makedirs('%s/run_%s' % (out_dir,now))  # Fail early if the dir already exists
    normalized = args.normalize # Set this as needed
    sampling = args.samp_freq #sampling frequency
    mismatch=.2
    qrange=(4,64)
    frange=(0,np.inf)

    # Read data and remove low frequency content
    fname = 'H-H1_LOSC_4_V2-1126259446-32.gwf'
    url = "https://losc.ligo.org/s/events/GW150914/" + fname
    urllib.urlretrieve(url, filename=fname)
    h1 = read_frame('H-H1_LOSC_4_V2-1126259446-32.gwf', 'H1:LOSC-STRAIN')
    h1 = TimeSeries(np.random.normal(size=64*4096), delta_t = 1. / sampling)
    h1 = highpass_fir(h1, 15, 8)

    # Calculate the noise spectrum
    psd = interpolate(welch(h1), 1.0 / 32)

    #perform Q-tiling
    Qbase, frange = qtiling(h1, qrange, frange, sampling, normalized, mismatch)

    #Choose Q-plane and plot
    qplane = Qplane(Qbase, h1, sampling, normalized, out_dir, now, frange)

    #Plot spectrogram
    plotter(qplane, out_dir, now, frange, h1, sampling)

    print 'Done!'

if __name__ == '__main__':
    main()
