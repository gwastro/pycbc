# -*- coding: utf-8 -*-
# Copyright (C) 2021 Sebastian Gaebel
#
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
Example to illustrate the use of the pycbc.noise.reproduceable.variable_noise
function.

@author: Sebastian Gaebel
@email: gaebel.sebastian@gmail.com
"""

from numpy import log
from scipy.interpolate import interp1d
from pycbc.types import FrequencySeries
from pycbc.noise.reproduceable import \
    generate_psd_fn_linear_in_log_transition, variable_noise
import matplotlib.pyplot as plt
import numpy as np


if __name__ == '__main__':
    # read from some plot very roughly by hand
    rough_LIGO_noise = (np.array([10., 50., 400., 10000.]),
                        np.array([1e-22, 4e-24, 2e-24, 1e-22]))
    rough_ET_noise = (np.array([1.3, 60., 300., 10000.]),
                      np.array([7e-22, 2e-25, 2e-25, 6e-24]))
    # turn into FrequencySeries
    freqs, amps = rough_LIGO_noise
    _interp = interp1d(x=log(freqs), y=log(amps), kind='linear',
                       bounds_error=True)
    def interpolator(f):
        return np.exp(_interp(log(f)))
    ligo_freqs = np.linspace(rough_LIGO_noise[0][0],  # f_start
                             rough_LIGO_noise[0][1],  # f_end
                             1000)  # resolution, somewhat arbitrary
    delta_f = ligo_freqs[1] - ligo_freqs[0]
    ligo_amps = interpolator(ligo_freqs)
    LIGO_noise = FrequencySeries(initial_array=ligo_amps, delta_f=delta_f)
    del interpolator

    freqs, amps = rough_ET_noise
    _interp = interp1d(x=log(freqs), y=log(amps), kind='linear',
                       bounds_error=True)
    def interpolator(f):
        return np.exp(_interp(log(f)))
    et_freqs = np.linspace(rough_ET_noise[0][0],  # f_start
                           rough_ET_noise[0][1],  # f_end
                           1000)  # resolution, somewhat arbitrary
    delta_f = et_freqs[1] - et_freqs[0]
    et_amps = interpolator(et_freqs)
    ET_noise = FrequencySeries(initial_array=et_amps, delta_f=delta_f)
    del interpolator

    psd_fn = generate_psd_fn_linear_in_log_transition(LIGO_noise, ET_noise)
    sample_rate = 16834.
    start_time = 12345678
    duration = 10000
    psd_resolution = 100000
    filter_duration = 32.
    time_resolution = -1000

    noise = variable_noise(psd_fn=psd_fn,
                           start_time=start_time,
                           end_time=start_time+duration,
                           psd_resolution=psd_resolution,
                           filter_duration=filter_duration,
                           seed=0,
                           sample_rate=sample_rate,
                           time_resolution=time_resolution,
                           progress=True)
    times, values = noise.sample_times, np.array(noise)

    # plot the result

    # Note:
    #   while the slopes all look fine, the overall levels in my plots are off.
    #   however since the relative noise levels are consistent, I'll chalk
    #   that up to messing up the units or something similar in the FFT plots.

    # 1st: overall noise
    freqs = np.fft.rfftfreq(times.size, 1./sample_rate)
    psd_0 = psd_fn(freqs, 0.)
    psd_1 = psd_fn(freqs, 1.)
    spectrum = np.abs(np.fft.rfft(values)) / 10**5.5
    spectrum = np.convolve(spectrum, np.ones(100), mode='same')  # apply a bit of smoothing

    plt.figure(figsize=(12, 8))
    plt.loglog(freqs, spectrum, label='Gen. Noise')
    plt.loglog(freqs, psd_0, label='Initial PSD')
    plt.loglog(freqs, psd_1, label='Final PSD')
    plt.legend()
    plt.grid(True, which='both')
    plt.xlim(7., 1e4)
    plt.ylim(1e-25, 1e-21)
    plt.tight_layout()
    plt.savefig('test-noise-gen.png', dpi=120)

    # 2nd: noise shapes at start, middle, and end in blocks
    block_size = int(duration * sample_rate / 10)  # 10% of total
    mid_start_idx = int(values.size//2-block_size)


    freqs = np.fft.rfftfreq(block_size, 1./sample_rate)
    psd_0 = psd_fn(freqs, 0.)
    psd_0p5 = psd_fn(freqs, 0.5)
    psd_1 = psd_fn(freqs, 1.)
    ts_0 = values[:block_size]
    ts_0p5 = values[mid_start_idx:mid_start_idx+block_size]
    ts_1 = values[-block_size:]
    spec_0 = np.abs(np.fft.rfft(ts_0)) / 10**5.
    spec_0 = np.convolve(spec_0, np.ones(100), mode='same')  # apply a bit of smoothing
    spec_0p5 = np.abs(np.fft.rfft(ts_0p5)) / 10**5.
    spec_0p5 = np.convolve(spec_0p5, np.ones(100), mode='same')  # apply a bit of smoothing
    spec_1 = np.abs(np.fft.rfft(ts_1)) / 10**5.
    spec_1 = np.convolve(spec_1, np.ones(100), mode='same')  # apply a bit of smoothing

    plt.figure(figsize=(12, 8))
    plt.plot(freqs, spec_0, 'b', label='Start Noise')
    plt.plot(freqs, psd_0, 'b--', label='Start PSD')
    plt.plot(freqs, spec_0p5, 'r', label='Mid Noise')
    plt.plot(freqs, psd_0p5, 'r--', label='Mid PSD')
    plt.plot(freqs, spec_1, 'g', label='End Noise')
    plt.plot(freqs, psd_1, 'g--', label='End PSD')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.grid(True, which='both')
    plt.xlim(7., 1e4)
    plt.ylim(1e-26, 1e-22)
    plt.tight_layout()
    plt.savefig('test-noise-gen-blocks.png', dpi=120)
