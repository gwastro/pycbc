""" Generation of sine-Gaussian bursty type things
"""

import pycbc.types
import numpy
import functools

@functools.lru_cache(maxsize=128)
def cached_arange(kmax, delta_f):
    return numpy.arange(0, kmax) * delta_f


def fd_sine_gaussian(amp, quality, central_frequency, fmin, fmax, delta_f):
    """ Generate a Fourier domain sine-Gaussian

    Parameters
    ----------
    amp: float
        Amplitude of the sine-Gaussian
    quality: float
        The quality factor
    central_frequency: float
        The central frequency of the sine-Gaussian
    fmin: float
        The minimum frequency to generate the sine-Gaussian. This determines
        the length of the output vector.
    fmax: float
        The maximum frequency to generate the sine-Gaussian
    delta_f: float
        The size of the frequency step

    Returns
    -------
    sg: pycbc.types.Frequencyseries
        A Fourier domain sine-Gaussian
    """
    # Optimization note: Ian has profiled and done optimization on this
    # function. If further speed up is needed caching the v vector and
    # avoiding a numpy.zeros call would be the next thing to speed up.
    # After that the numpy.exp call dominates.
    kmin = int(round(fmin / delta_f))
    kmax = int(round(fmax / delta_f))

    pi = numpy.pi
    tau = (quality / (2 * pi * central_frequency))
    quality_sq = quality**2

    f = cached_arange(kmax, delta_f)

    # exp(exp_term1) and exp(exp_term2) are often 0 (at double-precision
    # level) but still slow to compute. Want to shortcut this by not
    # computing terms at values where we don't need to. Use e**(-50) ~ 0 as
    # the point at which we no longer compute np.exp. Given that the maximum
    # value is O(1), e**(-50) / e**(-1) is 0 at double precision and this is
    # safe.

    # We first figure out which points we need to compute amplitudes for
    exp_term_cutoff = -50

    v = numpy.zeros(kmax, dtype=numpy.complex128)
    indices = numpy.zeros(kmax, dtype=bool)
    # Only consider points larger than kmin
    indices[kmin:] = 1
    # Find frequencies at which first term is equal to exp_term_cutoff
    low_freq_first_term = (
        central_frequency - (-exp_term_cutoff)**0.5 / (tau * pi)
    )
    high_freq_first_term = (
        central_frequency + (-exp_term_cutoff)**0.5 / (tau * pi)
    )
    low_freq_first_idx = max(kmin, int(low_freq_first_term//delta_f))
    high_freq_first_idx = min(kmax, int(high_freq_first_term//delta_f))
    # Find frequency at which second term drops to exp_term_cutoff
    high_freq_second_idx = (
        int(-exp_term_cutoff / quality_sq * central_frequency // delta_f)
    )

    exp_term_1 = -(
        tau * pi * 
        (f[low_freq_first_idx:high_freq_first_idx] - central_frequency)
    )**2.0

    A_term = amp * (pi**0.5) / 2 * tau

    v[low_freq_first_idx:high_freq_first_idx] = (
        A_term * numpy.exp(exp_term_1)
    )
    # If the first term is already less than e**50 don't need the second
    # term at all ... It's often the case that the second term is not needed.
    if high_freq_second_idx > kmin:
        exp_term_2 = (
            -quality_sq * f[kmin:high_freq_second_idx] / central_frequency
        )
        v[kmin:high_freq_second_idx] *= (1 + numpy.exp(exp_term_2))

    return pycbc.types.FrequencySeries(v, delta_f=delta_f, copy=False)

