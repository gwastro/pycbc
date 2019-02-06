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
from pycbc.types import TimeSeries, zeros

def complex_median(complex_list):
    """ Get the median value of a list of complex numbers.

    Parameters
    ----------
    complex_list: list
        List of complex numbers to calculate the median.

    Returns
    -------
    a + 1.j*b: complex number
        The median of the real and imaginary parts.
    """
    median_real = numpy.median([complex_number.real
                     for complex_number in complex_list])
    median_imag = numpy.median([complex_number.imag
                     for complex_number in complex_list])
    return median_real + 1.j*median_imag

def avg_inner_product(data1, data2, bin_size):
    """ Calculate the time-domain inner product averaged over bins.

    Parameters
    ----------
    data1: pycbc.types.TimeSeries
        First data set.
    data2: pycbc.types.TimeSeries
        Second data set, with same duration and sample rate as data1.
    bin_size: float
        Duration of the bins the data will be divided into to calculate
        the inner product.

    Returns
    -------
    inner_prod: list
        The (complex) inner product of data1 and data2 obtained in each bin.
    amp: float
        The absolute value of the median of the inner product.
    phi: float
        The angle of the median of the inner product.
    """
    assert data1.duration == data2.duration
    assert data1.sample_rate == data2.sample_rate
    seglen = int(bin_size * data1.sample_rate)
    inner_prod = []
    for idx in range(int(data1.duration / bin_size)):
        start, end = idx * seglen, (idx+1) * seglen
        norm = len(data1[start:end])
        bin_prod = 2 * sum(data1.data[start:end].real *
                           numpy.conjugate(data2.data[start:end])) / norm
        inner_prod.append(bin_prod)

    # Get the median over all bins to avoid outliers due to the presence
    # of a signal in a particular bin.
    inner_median = complex_median(inner_prod)
    return inner_prod, numpy.abs(inner_median), numpy.angle(inner_median)

def line_model(freq, data, tref, amp=1, phi=0):
    """ Simple time-domain model for a frequency line.

    Parameters
    ----------
    freq: float
        Frequency of the line.
    data: pycbc.types.TimeSeries
        Reference data, to get delta_t, start_time, duration and sample_times.
    tref: float
        Reference time for the line model.
    amp: {1., float}, optional
        Amplitude of the frequency line.
    phi: {0. float}, optional
        Phase of the frequency line (radians).

    Returns
    -------
    freq_line: pycbc.types.TimeSeries
        A timeseries of the line model with frequency 'freq'. The returned
        data are complex to allow measuring the amplitude and phase of the
        corresponding frequency line in the strain data. For extraction, use
        only the real part of the data.
    """
    freq_line = TimeSeries(zeros(len(data)), delta_t=data.delta_t,
                           epoch=data.start_time)

    times = data.sample_times - float(tref)
    alpha = 2 * numpy.pi * freq * times + phi
    freq_line.data = amp * numpy.exp(1.j * alpha)

    return freq_line

def matching_line(freq, data, tref, bin_size=1):
    """ Find the parameter of the line with frequency 'freq' in the data.

    Parameters
    ----------
    freq: float
        Frequency of the line to find in the data.
    data: pycbc.types.TimeSeries
        Data from which the line wants to be measured.
    tref: float
        Reference time for the frequency line.
    bin_size: {1, float}, optional
        Duration of the bins the data will be divided into for averaging.

    Returns
    -------
    line_model: pycbc.types.TimeSeries
        A timeseries containing the frequency line with the amplitude
        and phase measured from the data.
    """
    template_line = line_model(freq, data, tref=tref)
    # Measure amplitude and phase of the line in the data
    _, amp, phi = avg_inner_product(data, template_line,
                                    bin_size=bin_size)
    return line_model(freq, data, tref=tref, amp=amp, phi=phi)

def calibration_lines(freqs, data, tref=None):
    """ Extract the calibration lines from strain data.

    Parameters
    ----------
    freqs: list
        List containing the frequencies of the calibration lines.
    data: pycbc.types.TimeSeries
        Strain data to extract the calibration lines from.
    tref: {None, float}, optional
        Reference time for the line. If None, will use data.start_time.

    Returns
    -------
    data: pycbc.types.TimeSeries
        The strain data with the calibration lines removed.
    """
    if tref is None:
        tref = float(data.start_time)
    for freq in freqs:
        measured_line = matching_line(freq, data, tref,
                                      bin_size=data.duration)
        data -= measured_line.data.real

    return data

def clean_data(freqs, data, chunk, avg_bin):
    """ Extract time-varying (wandering) lines from strain data.

    Parameters
    ----------
    freqs: list
        List containing the frequencies of the wandering lines.
    data: pycbc.types.TimeSeries
        Strain data to extract the wandering lines from.
    chunk: float
        Duration of the chunks the data will be divided into to account
        for the time variation of the wandering lines. Should be smaller
        than data.duration, and allow for at least a few chunks.
    avg_bin: float
        Duration of the bins each chunk will be divided into for averaging
        the inner product when measuring the parameters of the line. Should
        be smaller than chunk.

    Returns
    -------
    data: pycbc.types.TimeSeries
        The strain data with the wandering lines removed.
    """
    if avg_bin >= chunk:
        raise ValueError('The bin size for averaging the inner product '
                         'must be less than the chunk size.')
    if chunk >= data.duration:
        raise ValueError('The chunk size must be less than the '
                         'data duration.')
    steps = numpy.arange(0, int(data.duration/chunk)-0.5, 0.5)
    seglen = chunk * data.sample_rate

    tref = float(data.start_time)
    for freq in freqs:
        for step in steps:
            start, end = int(step*seglen), int((step+1)*seglen)
            chunk_line = matching_line(freq, data[start:end],
                                       tref, bin_size=avg_bin)

            # Apply hann window on sides of chunk_line to smooth boundaries
            # and avoid discontinuities
            hann_window = numpy.hanning(len(chunk_line))
            apply_hann = TimeSeries(numpy.ones(len(chunk_line)),
                                    delta_t=chunk_line.delta_t,
                                    epoch=chunk_line.start_time)
            if step == 0:
                apply_hann.data[len(hann_window)/2:] *= \
                                hann_window[len(hann_window)/2:]
            elif step == steps[-1]:
                apply_hann.data[:len(hann_window)/2] *= \
                                hann_window[:len(hann_window)/2]
            else:
                apply_hann.data *= hann_window
            chunk_line.data *= apply_hann.data
            data.data[start:end] -= chunk_line.data.real

    return data
