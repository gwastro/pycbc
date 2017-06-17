# Module to hold recalibration stuff used during parameter estimation
# This will draw heavily from pycbc_adjust_strain and cal.py which
# are part of Chris Biwer's pycbc-cal repository

import numpy
import scipy
import pycbc

def update_c(fs=None, qinv=None, fc=None, freqs=None, c_res=None,
             kappa_c=1.0):
    detuning_term = freqs**2 / (freqs**2 - 1.0j * freqs * fs * qinv + fs**2)
    return c_res * kappa_c / (1 + 1.0j * freqs / fc) * detuning_term

def update_g(fs=None, qinv=None, a_tst0=None, a_pu0=None, fc=None,
             freqs=None, c_res=None, d0=None, kappa_tst_re=1.0,
             kappa_tst_im=0.0, kappa_pu_re=1.0, kappa_pu_im=0.0,
             kappa_c=1.0):
    c = update_c(fs=fs, qinv=qinv, fc=fc, freqs=freqs, c_res=c_res,
                 kappa_c=kappa_c)
    a_tst = a_tst0 * (kappa_tst_re + 1.0j * kappa_tst_im)
    a_pu = a_pu0 * (kappa_pu_re + 1.0j * kappa_pu_im)
    return c * d0 * (a_tst + a_pu)

def update_r(fs=None, qinv=None, a_tst0=None, a_pu0=None, fc=None,
             freqs=None, c_res=None, d0=None, kappa_c=1.0,
             kappa_tst_re=1.0, kappa_tst_im=0.0, kappa_pu_re=1.0,
             kappa_pu_im=0.0):
    c = update_c(fs=fs, qinv=qinv, fc=fc, freqs=freqs, c_res=c_res,
                 kappa_c=kappa_c)
    g = update_g(fs=fs, qinv=qinv, a_tst0=a_tst0, a_pu0=a_pu0, fc=fc,
                 freqs=freqs, c_res=c_res, d0=d0, kappa_c=kappa_c,
                 kappa_tst_re=kappa_tst_re, kappa_tst_im=kappa_tst_im,
                 kappa_pu_re=kappa_pu_re, kappa_pu_im=kappa_pu_im)
    return (1.0 + g) / c

def adjust_strain(strain, fc0=None, c0=None, d0=None, a_tst0=None,
                  a_pu0=None, fs0=None, qinv0=None, fs=None,
                  qinv=None, fc=None, freqs=None, kappa_c=1.0,
                  kappa_tst_re=1.0, kappa_tst_im=0.0,
                  kappa_pu_re=1.0, kappa_pu_im=0.0):
    """Adjust the FrequencySeries strain
    """
    g0 = c0 * d0 * (a_tst0 + a_pu0)
    r0 = (1.0 + g0) / c0
    init_detuning = freqs**2 / (freqs**2 - 1.0j * freqs * fs0 * qinv0 + fs0**2)
    c_res = c0 * (1 + 1.0j * freqs / fc0) / init_detuning

    r_adjusted = update_r(fs=fs, qinv=qinv, a_tst0=a_tst0, a_pu0=a_pu0, fc=fc,
                          freqs=freqs, c_res=c_res, d0=d0, kappa_c=kappa_c,
                          kappa_tst_re=kappa_tst_re, kappa_tst_im=kappa_tst_im,
                          kappa_pu_re=kappa_pu_re, kappa_pu_im=kappa_pu_im)

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
