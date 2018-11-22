""" Utilities for introducing nonlinear tidal effects into waveform approximants
"""
import pycbc.conversions
import numpy
import lal

def nltides_fourier_phase_difference(f, delta_f, f0, amplitude, n, m1, m2):
    """Calculate the change to the Fourier phase change due
    to non-linear tides. Note that the Fourier phase Psi(f)
    is not the same as the gravitational-wave phase phi(f) and
    is computed by
    Delta Psi(f) = 2 \pi f Delta t(f) - Delta phi(f)

    Parameters
    ----------
    f: numpy.array
        Array of frequency values to calculate the fourier phase difference
    delta_f: float
        Frequency resolution of f array
    f0: float
        Frequency that NL effects switch on
    amplitude: float
        Amplitude of effect
    n: float
        Growth dependence of effect
    m1: float
        Mass of component 1
    m2: float
        Mass of component 2

    Returns
    -------
    delta_psi: numpy.array
        Fourier phase as a function of frequency
    """

    kmin = int(f0/delta_f)
    kmax = len(f)

    f_ref, t_of_f_factor, phi_of_f_factor = \
        pycbc.conversions.nltides_coefs(amplitude, n, m1, m2)

    # Fourier phase shift below f0 from \Delta \phi(f)
    delta_psi_f_le_f0 = numpy.ones(kmin)
    delta_psi_f_le_f0 *= - phi_of_f_factor * (f0/f_ref)**(n-3.)

    # Fourier phase shift above f0 from \Delta \phi(f)
    delta_psi_f_gt_f0 = - phi_of_f_factor * (f[kmin:kmax]/f_ref)**(n-3.)

    # Fourier phase shift below f0 from 2 pi f \Delta t(f)
    delta_psi_f_le_f0 += 2.0 * lal.lal.PI * f[0:kmin] * t_of_f_factor * \
        (f0/f_ref)**(n-4.)

    # Fourier phase shift above f0 from 2 pi f \Delta t(f)
    delta_psi_f_gt_f0 += 2.0 * lal.lal.PI * f[kmin:kmax] * t_of_f_factor * \
        (f[kmin:kmax]/f_ref)**(n-4.)

    # Return the shift to the Fourier phase
    return numpy.concatenate((delta_psi_f_le_f0, delta_psi_f_gt_f0), axis=0)


def nonlinear_tidal_spa(**kwds):
    """Generates a frequency-domain waveform that implements the
    TaylorF2+NL tide model described in https://arxiv.org/abs/1808.07013
    """

    from pycbc import waveform
    from pycbc.types import Array

    # We start with the standard TaylorF2 based waveform
    kwds.pop('approximant')
    hp, hc = waveform.get_fd_waveform(approximant="TaylorF2", **kwds)

    # Add the phasing difference from the nonlinear tides
    f = numpy.arange(len(hp)) * hp.delta_f
    pd =  Array(numpy.exp(-1.0j * nltides_fourier_phase_difference(f,
               hp.delta_f,
               kwds['f0'], kwds['amplitude'], kwds['n'],
               kwds['mass1'], kwds['mass2'])),
               dtype=hp.dtype)
    hp *= pd
    hc *= pd
    return hp, hc
