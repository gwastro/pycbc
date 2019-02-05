# Copyright (C) 2019 Miriam Cabero
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
""" Functions for removing frequency lines from real data.
"""

import numpy
from pycbc.types import TimeSeries, FrequencySeries, Array, complex128
from pycbc.types.array import complex_same_precision_as, zeros


def imag_median(complex_list):
    x = numpy.median([complex_number.real for complex_number in complex_list])
    y = numpy.median([complex_number.imag for complex_number in complex_list])
    return x + 1.j*y

def avg_inner_product(data1, data2, bin_size):
    # Calculate the time-domain inner product in each bin of duration bin_size
    # and return the median from all bins to avoid outliers due to the presence of 
    # a signal in a particular bin
    assert data1.duration == data2.duration
    assert data1.sample_rate == data2.sample_rate
    seglen = int(bin_size * data1.sample_rate)
    inner_prod = []
    for n in range(int(data1.duration / bin_size)):
        start, end = n * seglen, (n+1) * seglen
        norm = len(data1[start:end])
        bin_prod = 2 * sum(data1.data[start:end].real * \
                    numpy.conjugate(data2.data[start:end])) / norm
        inner_prod.append(bin_prod)

    inner_median = imag_median(inner_prod)
    return inner_prod, numpy.abs(inner_median), numpy.angle(inner_median)

def line_model(freq, data, tref, amp=1, phi=0):
    # Simple time-domain model for the frequency line, 
    # with the same duration as data.
    # The returned data are complex to allow measuring the amplitude and phase
    # of the strain data, for extracting purposes use only the real part
    t = data.sample_times
    freq_line = TimeSeries(zeros(len(data)), delta_t=data.delta_t, 
                                             epoch=data.start_time)

    alpha = 2 * numpy.pi * freq * (t - tref) + phi
    freq_line.data = amp * numpy.exp(1.j * alpha)

    return freq_line

def matching_line(freq, data, tref, bin_size=1):
    template_line = line_model(freq, data, tref=tref)
    # Measure amplitude and phase of the line in the data
    _, amp, phi = avg_inner_product(data, template_line,
                                    bin_size=bin_size)
    return line_model(freq, data, tref=tref, amp=amp, phi=phi)

def calibration_lines(freqs, data, tref=None):
    # Calibration lines are stable and the frequency is well know, so we use
    # a different (simpler) method to remove them. Give freqs as a list
    if tref is None:
        tref = float(data.start_time)
    for freq in freqs:
        measured_line = matching_line(freq, data, tref,
                                      bin_size=data.duration)
        data -= measured_line.data.real

    return data

def clean_data(freqs, data, chunk, avg_bin):
    # In general, noise lines are time-varying (wandering). To account for this
    # time variation, divide data into chunks of size chunk > bin_size 
    # (bin_size is the segment duration for averaging the inner product)
    tref = float(data.start_time)
    if avg_bin >= chunk:
        raise ValueError('The bin size for averaging the inner product must ' \
                         'be less than the chunk size.')
    steps = numpy.arange(0, int(data.duration/chunk)-0.5, 0.5)
    seglen = chunk * data.sample_rate

    for freq in freqs:
        for n in steps:
            start, end = int(n*seglen), int((n+1)*seglen)
            chunk_line = matching_line(freq, data[start:end], 
                                       tref, bin_size=avg_bin)
            
            # Apply hann window on sides of chunk_line to smooth boundaries
            # and avoid discontinuities
            hann_window = numpy.hanning(len(chunk_line))
            apply_hann = TimeSeries(numpy.ones(len(chunk_line)), 
                                    delta_t=chunk_line.delta_t, 
                                    epoch=chunk_line.start_time)
            if n == 0:
                apply_hann.data[len(hann_window)/2:] *= \
                                hann_window[len(hann_window)/2:]
            elif n == steps[-1]:
                apply_hann.data[:len(hann_window)/2] *= \
                                hann_window[:len(hann_window)/2]
            else:
                apply_hann.data *= hann_window
            chunk_line.data *= apply_hann.data
            data.data[start:end] -= chunk_line.data.real
            
    return data
