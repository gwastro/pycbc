#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 12:54:12 2021

Generates colored noise for which the noise PSD may change over time.

@author: Sebastian Gaebel
@email: gaebel.sebastian@gmail.com⎄
"""

import numpy as np
from numpy import log
from scipy.interpolate import interp1d
from scipy.signal import firwin2
from warnings import warn as warn_fn


def _compute_coefficients(freqs, amps, nyquist, n_taps):
    """
    Compute the FIR filter coefficients given a filter via frequency series.

    Parameters
    ----------
    freqs : (N,) array
        Frequencies at which the filter is defined.
    amps : (N,) array
        PSD amplitudes of the filter.
    nyquist : float
        Nyquist frequency of the time series to be filtered,
        `nyquist = sampling_frequency / 2`.
    n_taps : int
        Number of "taps", i.e. FIR filter coefficients. A higher value
        leads to a more precise approximation of the desired filter.

    Returns
    -------
    coefficients : (n_taps,) array
        FIR filter coefficients as computed via `scipy.signal.firwin2`.
    """
    # not an expert on the details, but firwin2 requires the frequencies to
    #   go from 0 to 1 (nyquist scale) and the gain to be 0 at the endpoints
    f_scaled = freqs / nyquist  # firwin2 works in units of nyquist
    f_scaled[0] = 0. # similarly, frequencies need to start at 0
    f_scaled[-1] = 1.  # brute force the frequency value, floating point math sucks
    amps[0] = 0.
    amps[-1] = 0.
    assert f_scaled[-1] == 1., repr(f_scaled[-1])  # last element should be 1
    coefficients = firwin2(numtaps=n_taps, freq=f_scaled, gain=amps)
    return coefficients


def generate_time_varying_noise(
        psd_fn, sample_rate, duration, psd_resolution, n_taps,
        time_resolution, use_time_directly=False, seed=0,
        quiet=False, skip_checks=False, progress=False):
    """ Create a time series of pure noise with time dependent noise spectrum.

    Parameters
    ----------
    psd_fn : callable
        Function which provides the time varying PSD. Must accept a call
        signature of `psd_fn(freqs, tau)`, where `freqs` is a numpy array
        of frequencies in [Hz], and `tau` is the normalised time, i.e.
        `current_time / duration`. Must return a numpy array of the same
        shape as `freqs`.
    sample_rate : float
        Sample rate of the generate time series in units of [sec]. Required
        as the filter design depends on the Nyquist frequency.
    duration : float
        Duration of the resulting time series in [sec].
    psd_resolution : int
        Number of points for which the psd_fn in evaluated before being
        passed to the filter generation. Using a large number should be
        generally harmless and without side effects. Recommended to be
        a few time the filter resolution `n_taps`.
    n_taps : int
        Number of coefficients used in the FIR filter, i.e. length of the
        filter in samples.
    time_resolution : int
        The FIR filter is recomputed every `time_resolution` samples to
        smoothly vary the noise PSD over the length of the time series.
        This directly impacts performance as filter application is much
        cheaper than the recomputing of filter coefficients.
        If the value is larger than the number of samples, a warning is
        issued if `quiet` is not set to `True`.
    use_time_directly : bool
        If `True`, the `tau` parameter in `psd_fn(freq, tau)` will be
        passed as the true time within the generated time series.
        The default is False.
    seed : int, optional
        Seed used to initialise the numpy random number generator.
        The default is 0.
    quiet : bool, optional
        If `True`, warnings will be supressed. The default is False.
    skip_checks : bool, optional
        If `True`, some internal consistency checks will be skipped for
        the sake of slightly better performance.
    progress : bool, optional
        If `True`, display a simple progress counter. The default is False.

    Returns
    -------
    times, noisy_time_series : pair of (duration*sample_rate,) arrays
        Time stamps and time series of pure noise with variable
        noise properties.

    Notes
    -----
    Currently, the overall noise level may be off by a constant factor.
    Check this and remove the warning before serious use.

    Reference for the FIR filter application:
        https://en.wikipedia.org/wiki/Finite_impulse_response#Definition

    """
    def warn(*args, **kwargs):
        """
        Internal helper to issue qarning only if `quiet` is `False`.

        Parameters
        ----------
        *args and **kwargs directly passed on to `warnings.warn` iff
        `quiet` is Falsey.

        Returns
        -------
        None.

        """
        if not quiet:
            warn_fn(*args, **kwargs)
        return

    def zero_pad_left(arr, n):
        """
        Zero pad an array to a given size. Zeros are added on the left,
        i.e. at the beginning.

        Parameters
        ----------
        arr : (N,) array
            Some array of floats.
        n : int
            Size to which the array is to be padded.

        Returns
        -------
        padded : (n,) array
            Left padded array.

        """
        if arr.size >= n:
            return arr
        padded = np.concatenate([np.zeros(n - arr.size), arr])
        return padded

    # some small value to use in place of zeros to avoid trouble
    #   while working in log space
    zero_replacement = 1e-100

    # first: compute some essential constants
    nyquist = sample_rate / 2.  # [Hz]
    delta_t = 1. / sample_rate  # [sec], time between samples
    n_samples = int(duration * sample_rate)
    freqs = np.geomspace(zero_replacement, nyquist, psd_resolution)

    # generate white noise, which will then be filtered to
    #   obtain the desired PSD.
    RNG = np.random.default_rng(seed)
    white_time_series = RNG.standard_normal(size=n_samples)
    times = np.arange(0., duration, delta_t)
    assert times.size == n_samples
    assert white_time_series.size == n_samples

    # prepare array for filter output
    noisy_time_series = np.zeros_like(white_time_series)

    # precompute the break points for re-computation of the FIR coefficients
    if use_time_directly:
        # simply compute the number of re-evals and replace the value
        #   of time_resolution to avoid new code
        time_resolution = int(duration / time_resolution)
        # technically the shifts the breakspoint due to rounding
        #   but that really ought not matter

    # This loop is a prime candidate for optimisation,
    #   e.g. pre-computing FIR coefficients and then dropping into cython
    #   for the convolution.
    for idx in range(n_samples):
        tau = float(idx) / n_samples
        if use_time_directly:
            tau = idx * delta_t
        if idx % time_resolution == 0:
            if progress:
                print('Recompute at sample idx {} of {} ({:.1%})'
                      ''.format(idx, n_samples, idx/n_samples), flush=True)
            # re-compute the coefficients
            amps = psd_fn(freqs, tau)
            if not skip_checks:
                assert amps.size == psd_resolution
            coeffs = _compute_coefficients(freqs, amps, nyquist, n_taps)
            # as a quirk of definition, the coefficients are time reversed
            coeffs = coeffs[::-1]
            if not skip_checks:
                assert coeffs.size == n_taps
        # select the chunk of the time series which is used by the filter
        chunk = zero_pad_left(white_time_series[max(idx-n_taps, 0)+1 : idx+1],
                              n_taps)
        #   max-term : avoid selecting before element 0
        #   idx+1 : the current sample needs to be included
        # convolution
        noisy_time_series[idx] = np.sum(chunk * coeffs)

    return times, noisy_time_series


# example for creating the input `psd_fn`

def example_pdf_fn_ligo_to_et_generator(gap):
    """
    Example of a psd_fn which transitions from a LIGO-like PSD to an
    ET-like PSD linearly in log-log space.

    Parameters
    ----------
    gap : float
        Transition width to the fill value.

    Returns
    -------
    psd_fn : callable

    """
    # read from some plot very roughly by hand
    rough_LIGO_noise = (np.array([10., 50., 400., 10000.]),
                        np.array([1e-22, 4e-24, 2e-24, 1e-22]))
    rough_ET_noise = (np.array([1.3, 60., 300., 10000.]),
                      np.array([7e-22, 2e-25, 2e-25, 6e-24]))

    PSD_FILLER_VALUE = 1e-100

    def insert_transition(pair, gap):
        """ Transition the gain to the filler value at the edge of their
        defined interval. Using this is less sharp than a instant drop to
        the fillter value than the interpolation function would yield.

        Parameters
        ----------
        pair : frequency & gain pair
        gap : float
            Distance over which the gain is transitioned to the filler value.

        Returns
        -------
        psd_fn : callable
            PSD function as required for the use in
            `generate_time_varying_noise` above.

        """
        f, a = pair
        f = np.concatenate([[f[0] - gap], f, [f[-1] + gap]])
        a = np.concatenate([[PSD_FILLER_VALUE], a, [PSD_FILLER_VALUE]])
        return f, a

    # LIGO-like
    f, a = insert_transition(rough_LIGO_noise, 1e-3)
    psd_0 = interp1d(x=log(f), y=log(a), kind='linear', bounds_error=False,
                     fill_value=log(PSD_FILLER_VALUE))
    # ET-like
    f, a = insert_transition(rough_ET_noise, 1e-3)
    psd_1 = interp1d(x=log(f), y=log(a), kind='linear', bounds_error=False,
                     fill_value=log(PSD_FILLER_VALUE))

    def psd_fn(freqs, tau):
        # mix linearly
        a_0 = tau * psd_0(log(freqs))
        a_1 = (1. - tau) * psd_1(log(freqs))
        amps = np.exp(a_0 + a_1)
        return amps
    return psd_fn


if __name__ == '__main__':
    # some simple tests
    import matplotlib.pyplot as plt

    gap = 1e-3
    psd_fn = example_pdf_fn_ligo_to_et_generator(gap=gap)
    sample_rate = 16834.
    duration = 10000.
    psd_resolution = 6400
    n_taps = 2048
    time_resolution = 99999

    times, values = generate_time_varying_noise(psd_fn=psd_fn,
                                   sample_rate=sample_rate,
                                   duration=duration,
                                   psd_resolution=psd_resolution,
                                   n_taps=n_taps,
                                   time_resolution=time_resolution,
                                   progress=True)

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
