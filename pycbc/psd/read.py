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
""" Utilites to read PSDs from files. """

import numpy
import scipy.interpolate
from pycbc.types import FrequencySeries

def from_txt(filename, length, delta_f, low_freq_cutoff):
    """Returns a PSD from a two-column ASCII file containing
    frequency on the first column and ASD on the second.
    """
    asd_data = numpy.loadtxt(filename)
    if (asd_data < 0).any() or numpy.logical_not(numpy.isfinite(asd_data)).any():
        raise ValueError('Invalid ASD data in ' + filename)
    psd_interp = scipy.interpolate.interp1d(asd_data[:, 0], asd_data[:, 1] ** 2)

    kmin = int(low_freq_cutoff / delta_f)
    psd = numpy.zeros(length, dtype=numpy.float64)
    for k in xrange(kmin, length):
        psd[k] = float(psd_interp(k * delta_f))

    return FrequencySeries(psd, delta_f=delta_f)

