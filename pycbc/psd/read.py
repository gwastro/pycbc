#!/usr/bin/python
# Copyright (C) 2012 Alex Nitz, Tito Dal Canton
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
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
"""Utilites to read PSDs from files.
"""

import numpy
import scipy.interpolate
from pycbc.types import FrequencySeries

def from_asd_txt(filename, length, delta_f, low_freq_cutoff):
    """Read an ASCII file containing one-sided ASD data and generate
    a frequency series with the corresponding PSD. The ASD data is
    interpolated in order to match the desired resolution of the
    generated frequency series.

    Parameters
    ----------
    filename : string
        Path to a two-column ASCII file. The first column must contain
        the frequency (positive frequencies only) and the second column
        must contain the amplitude density.
    length : int
        Length of the frequency series in samples.
    delta_f : float
        Frequency resolution of the frequency series in Herz.
    low_freq_cutoff : float
        Frequencies below this value are set to zero.

    Returns
    -------
    psd : FrequencySeries
        The generated frequency series.

    Raises
    ------
    ValueError
        If the ASCII file contains negative, infinite or NaN frequencies
        or amplitude densities.
    """
    asd_data = numpy.loadtxt(filename)
    if (asd_data < 0).any() or numpy.logical_not(numpy.isfinite(asd_data)).any():
        raise ValueError('Invalid ASD data in ' + filename)

    flog = numpy.log10(asd_data[:, 0])
    slog = numpy.log10(asd_data[:, 1] ** 2)

    psd_interp = scipy.interpolate.interp1d(flog, slog)

    kmin = int(low_freq_cutoff / delta_f)
    psd = numpy.zeros(length, dtype=numpy.float64)
    for k in xrange(kmin, length):
        psd[k] = 10. ** float(psd_interp(numpy.log10(k * delta_f)))

    return FrequencySeries(psd, delta_f=delta_f)

