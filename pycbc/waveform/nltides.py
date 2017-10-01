""" Utilities for introducing nonlinear tidal effects into waveform approximants
"""
import pycbc.pnutils
import numpy

def nonlinear_phase_difference(f, f0, A, n, m1, m2):
    """Return the phase difference for nonlinear tides

    Implemenents the phase difference approximation of nonlinear
    tides in Essick, et al. https://arxiv.org/pdf/1609.06362v2.pdf

    Parameters
    ----------
    f: numpy.ndarray
        Array of frequency values to calculate the fourier phase difference.
    f0: float
        Frequency that NL effects switch on
    A: float
        Amplitude of effect
    n: float
        Growth dependence of effect
    m1: float
        Mass of component 1
    m2: float
        Mass of component 2
    """
    x0 = f0 / 100.0
    x = f / 100.0
    mc, _ = pycbc.pnutils.mass1_mass2_to_mchirp_eta(m1, m2)
    dphi = 0.4 * (mc / 1.2) ** (-10.0 / 3) * A * 10**8
    dphi *= (x0**(n - 3.0) - x**(n - 3.0)) / (n - 3.0)
    return dphi

def nonlinear_tidal_spa(**kwds):
    from . import waveform
    from pycbc.types import Array

    # We start with the standard TaylorF2 based waveform
    kwds.pop('approximant')
    hp, hc = waveform.get_fd_waveform(approximant="TaylorF2", **kwds)

    # Add the phasing difference from the nonlinear tides
    kmin = int((kwds['f0'] / hp.delta_f))
    f = numpy.arange(kmin, len(hp)) * hp.delta_f
    pd =  Array(numpy.exp(1.0j * nonlinear_phase_difference(f,
               kwds['f0'], kwds['A'], kwds['n'], kwds['mass1'], kwds['mass2'])),
               dtype=hp.dtype)
    hp[kmin:] *= pd
    hc[kmin:] *= pd
    return hp, hc
