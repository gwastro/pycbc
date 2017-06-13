# Module to hold recalibration stuff used during parameter estimation
# This will draw heavily from pycbc_adjust_strain and cal.py which
# are part of Chris Biwer's pycbc-cal repository

import numpy
import scipy
import pycbc

def tf_from_file(path, delimiter=" "):
    """ Convert the contents of a file with the columns
    [freq, real(h), imag(h)] to a numpy.array with columns
    [freq, real(h)+j*imag(h)].
    """
    data = numpy.loadtxt(path, delimiter=delimiter)
    freq = data[:, 0]
    h = data[:, 1] + 1.0j * data[:, 2]
    return numpy.array([freq, h]).transpose()

def update_c(fs=None, qinv=None, fc0=None, freqs=None, c_res=None):
    detuning_term = freqs**2 / (freqs**2 - 1.0j * freqs * fs * qinv + fs**2)
    return c_res * 1.0 / (1 + 1.0j * freqs / fc0) * detuning_term

def update_g(fs=None, qinv=None, a_tst0=None, a_pu0=None, fc0=None,
             freqs=None, c_res=None, d0=None):
    c = update_c(fs=fs, qinv=qinv, fc0=fc0, freqs=freqs, c_res=c_res)
    return c * d0 * (a_tst0 + a_pu0)

def update_r(fs=None, qinv=None, a_tst0=None, a_pu0=None, fc0=None,
             freqs=None, c_res=None, d0=None):
    c = update_c(fs=fs, qinv=qinv, fc0=fc0, freqs=freqs, c_res=c_res)
    g = update_g(fs=fs, qinv=qinv, a_tst0=a_tst0, a_pu0=a_pu0, fc0=fc0,
                 freqs=freqs, c_res=c_res, d0=d0)
    return (1.0 + g) / c

def adjust_strain(strain, fc0=None, c0=None, d0=None, a_tst0=None,
                  a_pu0=None, fs0=None, qinv0=None,
                  fs=None, qinv=None, freqs=None):
    """Adjust the FrequencySeries strain
    """
    g0 = c0 * d0 * (a_tst0 + a_pu0)
    r0 = (1.0 + g0) / c0
    init_detuning = freqs**2 / (freqs**2 - 1.0j * freqs * fs0 * qinv0 + fs0**2)
    c_res = c0 * (1 + 1.0j * freqs / fc0) / init_detuning
    r_adjusted = update_r(fs=fs, qinv=qinv, a_tst0=a_tst0, a_pu0=a_pu0, fc0=fc0,
                          freqs=freqs, c_res=c_res, d0=d0)

    # calculate error function
    k = r_adjusted / r0
    # decompose into amplitude and unwrapped phase
    k_amp = numpy.abs(k)
    k_phase = numpy.unwrap(numpy.angle(k))

    # convert to FrequencySeries by interpolating then resampling
    order = 1
    k_amp_off = scipy.interpolate.UnivariateSpline(freqs, k_amp, k=order,
                                                   s=0)
    k_phase_off = scipy.interpolate.UnivariateSpline(freqs, k_phase,
                                                     k=order, s=0)
    freq_even = strain.sample_frequencies.numpy()
    k_even_sample = k_amp_off(freq_even) * \
                    numpy.exp(1.0j * k_phase_off(freq_even))
    strain_adjusted = pycbc.types.FrequencySeries(
                      strain.numpy() * k_even_sample, delta_f=strain.delta_f)

    return strain_adjusted
