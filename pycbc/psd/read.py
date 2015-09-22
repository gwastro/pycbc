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
"""Utilities to read PSDs from files.
"""

import numpy
import scipy.interpolate
from pycbc.types import FrequencySeries

def from_txt(filename, length, delta_f, low_freq_cutoff, is_asd_file=True):
    """Read an ASCII file containing one-sided ASD or PSD  data and generate
    a frequency series with the corresponding PSD. The ASD or PSD data is
    interpolated in order to match the desired resolution of the
    generated frequency series.

    Parameters
    ----------
    filename : string
        Path to a two-column ASCII file. The first column must contain
        the frequency (positive frequencies only) and the second column
        must contain the amplitude density OR power spectral density.
    length : int
        Length of the frequency series in samples.
    delta_f : float
        Frequency resolution of the frequency series in Herz.
    low_freq_cutoff : float
        Frequencies below this value are set to zero.
    is_asd_file : Boolean
        If false assume that the second column holds power spectral density.
        If true assume that the second column holds amplitude spectral density.
        Default: True

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
    file_data = numpy.loadtxt(filename)
    if (file_data < 0).any() or \
                            numpy.logical_not(numpy.isfinite(file_data)).any():
        raise ValueError('Invalid data in ' + filename)

    freq_data = file_data[:, 0]
    noise_data = file_data[:, 1]
    
    # Only include points above the low frequency cutoff
    if freq_data[0] > low_freq_cutoff:
        raise ValueError('Lowest frequency in input file ' + filename + \
          ' is higher than requested low-frequency cutoff ' + str(low_freq_cutoff))
    
    kmin = int(low_freq_cutoff / delta_f)    
    flow = kmin * delta_f  
          
    data_start = (0 if freq_data[0]==low_freq_cutoff else numpy.searchsorted(freq_data, flow) - 1)
    
    # If the cutoff is exactly in the file, start there
    if freq_data[data_start+1] == low_freq_cutoff:
        data_start += 1
    
    freq_data = freq_data[data_start:]
    noise_data = noise_data[data_start:]    

    flog = numpy.log(freq_data)
    if is_asd_file:
        slog = numpy.log(noise_data ** 2)
    else:
        slog = numpy.log(noise_data)

    psd_interp = scipy.interpolate.interp1d(flog, slog)

    kmin = int(low_freq_cutoff / delta_f)
    psd = numpy.zeros(length, dtype=numpy.float64)

    vals = numpy.log(numpy.arange(kmin, length) * delta_f) 
    psd[kmin:] =  numpy.exp(psd_interp(vals))

    return FrequencySeries(psd, delta_f=delta_f)
