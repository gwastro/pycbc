# Copyright (C) 2015  Alex Nitz
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
""" Utilities for handling frequency compressed an unequally spaced frequency
domain waveforms.
"""
import lalsimulation, lal, numpy
from pycbc import pnutils
from pycbc.types import Array
from pycbc.opt import omp_libs, omp_flags
from pycbc import WEAVE_FLAGS
from scipy.weave import inline
from pycbc.types import FrequencySeries, zeros

def rough_time_estimate(m1, m2, flow, fudge_length=1.1, fudge_min=0.02):
    """ A very rough estimate of the duration of the waveform.

    An estimate of the waveform duration starting from flow. This is intended
    to be fast but not necessarily accurate. It should be an overestimate of the
    length. It is derived from a simplification of the 0PN post-newtonian
    terms and includes a fudge factor for possible ringdown, etc.

    Parameters
    ----------
    m1: float
        mass of first component object in solar masses
    m2: float
        mass of second component object in solar masses
    flow: float
        starting frequency of the waveform
    fudge_length: optional, {1.1, float}
        Factor to multiply length estimate by to ensure it is a convservative value
    fudge_min: optional, {0.2, float}
        Minimum signal duration that can be returned. This should be long enough
    to encompass the ringdown and errors in the precise end time.             

    Returns
    -------
    time: float
        Time from flow untill the end of the waveform
    """
    m = m1 + m2
    msun = m * lal.MTSUN_SI
    t =  5.0 / 256.0 * m * m * msun / (m1 * m2) / (numpy.pi * msun * flow) **  (8.0 / 3.0)

    # fudge factoriness
    return .022 if t < 0 else (t + fudge_min) * fudge_length 

def rough_frequency_samples(m1, m2, flow, fmax, df_min):
    """ Return a quick estimate of the frequency values which are needed to reproduce
    a waveform. Results are integer multiples of the df_min.

    Parameters
    ----------
    m1: float
        mass of first component object in solar masses
    m2: float
        mass of second component object in solar masses
    flow: float
        starting frequency of the waveform
    fmax: float
        ending frequency of the waveform
    df_min: float
        The size of a frequency sample step.
    """
    kmin = int(flow / df_min)
    kmax = int(fmax / df_min)
    k = kmin
    ksamples = []
    while k < kmax:
        ksamples.append(k)
        k += int(1.0 / rough_time_estimate(m1, m2, k * df_min) / df_min)
    ksamples.append(kmax)
    return numpy.array(ksamples) 
   
