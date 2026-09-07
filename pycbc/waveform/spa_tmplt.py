#  Adapted from code in LALSimInspiralTaylorF2.c
#
#  Copyright (C) 2007 Jolien Creighton, B.S. Sathyaprakash, Thomas Cokelaer
#  Copyright (C) 2012 Leo Singer, Alex Nitz
#
#  This program is free software you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with with program see the file COPYING. If not, write to the
#  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#  MA  02111-1307  USA

"""This module contains functions for generating common SPA template precalculated
   vectors.
"""
from math import sqrt, log
import warnings
import numpy

import pycbc.pnutils
from pycbc.scheme import schemed
from pycbc.types import FrequencySeries, Array, complex64, float32, zeros
from pycbc.waveform.utils import ceilpow2
from pycbc.constants import PI, MTSUN_SI, PC_SI, MRSUN_SI
from pycbc.libutils import import_optional

lal = import_optional('lal')
lalsimulation = import_optional('lalsimulation')

def findchirp_chirptime(m1, m2, fLower, porder=-1, s1z=0., s2z=0.):
    """Estimate the chirp time, i.e. the time from ``fLower`` to coalescence,
    of a TaylorF2 / stationary-phase-approximation waveform.

    The post-Newtonian coefficients are taken directly from LAL's
    ``SimInspiralTaylorF2AlignedPhasing`` rather than being hardcoded, so
    aligned-spin contributions to the phasing are included.

    Parameters
    ----------
    m1, m2 : float
        Component masses in solar masses.
    fLower : float or numpy.ndarray
        Lower frequency cutoff in Hz.
    porder : int, optional
        Twice the post-Newtonian order of the phasing (e.g. 7 for 3.5PN).
        The default (-1) lets LAL use the highest implemented order.
    s1z, s2z : float, optional
        Dimensionless spin components aligned with the orbital angular
        momentum. Default to zero (non-spinning).
    """
    m1 = float(m1)
    m2 = float(m2)
    m_sec = (m1 + m2) * MTSUN_SI
    eta = m1 * m2 / (m1 + m2) ** 2

    lal_pars = lal.CreateDict()
    if porder != -1:
        # otherwise LAL defaults to its highest implemented order, matching
        # the behaviour of spa_tmplt with phase_order=-1
        lalsimulation.SimInspiralWaveformParamsInsertPNPhaseOrder(
            lal_pars, porder)
    phasing = lalsimulation.SimInspiralTaylorF2AlignedPhasing(
        m1, m2, float(s1z), float(s2z), lal_pars)

    # PN expansion parameter v evaluated at the lower frequency cutoff
    v = (PI * m_sec * fLower) ** (1.0 / 3.0)
    lnv = numpy.log(v)

    pfaN = phasing.v[0]

    # In the stationary phase approximation the time-frequency relation is
    # t(f) = (1 / 2 pi) dPsi/df, so a phasing term (phi_k + phi_kl ln v) v^k
    # contributes to the chirp time tC = t_coalescence - t(fLower) with a
    # factor (5 - k) / 5 relative to the Newtonian term, plus an extra
    # -phi_kl / 5 v^k piece from differentiating the logarithm.
    series = 1.0
    for k in range(2, 8):
        phi_k = phasing.v[k] / pfaN
        phi_kl = phasing.vlogv[k] / pfaN
        series = series + v ** k * (
            (5.0 - k) / 5.0 * (phi_k + phi_kl * lnv) - phi_kl / 5.0)

    tN = 5.0 * m_sec / (256.0 * eta * v ** 8)
    return tN * series


def spa_length_in_time(**kwds):
    """
    Returns the length in time of the template,
    based on the masses, PN order, and low-frequency
    cut-off.
    """
    m1 = kwds['mass1']
    m2 = kwds['mass2']
    flow = kwds['f_lower']
    porder = int(kwds['phase_order'])
    s1z = kwds.get('spin1z') or 0.
    s2z = kwds.get('spin2z') or 0.

    return findchirp_chirptime(m1, m2, flow, porder, s1z=s1z, s2z=s2z)


def spa_amplitude_factor(**kwds):
    m1 = kwds['mass1']
    m2 = kwds['mass2']

    _, eta = pycbc.pnutils.mass1_mass2_to_mchirp_eta(m1, m2)

    FTaN = 32. * eta * eta / 5.
    dETaN = 2. * -eta / 2.

    M = m1 + m2

    m_sec = M * MTSUN_SI
    piM = PI * m_sec

    amp0 = 4. * m1 * m2 / (1e6 * PC_SI) * MRSUN_SI * MTSUN_SI * sqrt(PI / 12.)

    fac = numpy.sqrt(-dETaN / FTaN) * amp0 * (piM ** (-7./6.))
    return -fac


_prec = None
def spa_tmplt_precondition(length, delta_f, kmin=0):
    """Return the amplitude portion of the TaylorF2 approximant, used to precondition
    the strain data. The result is cached, and so should not be modified, only read.
    """
    global _prec
    if _prec is None or _prec.delta_f != delta_f or len(_prec) < length:
        v = numpy.arange(0, (kmin + length*2), 1.) * delta_f
        v = numpy.power(v[1:len(v)], -7./6.)
        _prec = FrequencySeries(v, delta_f=delta_f, dtype=float32)
    return _prec[kmin:kmin + length]


def spa_tmplt_norm(psd, length, delta_f, f_lower):
    amp = spa_tmplt_precondition(length, delta_f)
    k_min = int(f_lower / delta_f)
    sigma = (amp[k_min:length].numpy() ** 2. / psd[k_min:length].numpy())
    norm_vec = numpy.zeros(length)
    norm_vec[k_min:length] = sigma.cumsum() * 4. * delta_f
    return norm_vec


def spa_tmplt_end(**kwds):
    return pycbc.pnutils.f_SchwarzISCO(kwds['mass1'] + kwds['mass2'])


def spa_distance(psd, mass1, mass2, lower_frequency_cutoff, snr=8):
    """ Return the distance at a given snr (default=8) of the SPA TaylorF2
    template.
    """
    kend = int(spa_tmplt_end(mass1=mass1, mass2=mass2) / psd.delta_f)
    norm1 = spa_tmplt_norm(psd, len(psd), psd.delta_f, lower_frequency_cutoff)
    norm2 = spa_amplitude_factor(mass1=mass1, mass2=mass2) ** 2.0

    if kend >= len(psd):
        kend = len(psd) - 2
    return sqrt(norm1[kend] * norm2) / snr


@schemed("pycbc.waveform.spa_tmplt_")
def spa_tmplt_engine(htilde, kmin, phase_order, delta_f, piM, pfaN,
                     pfa2, pfa3, pfa4, pfa5, pfl5,
                     pfa6, pfl6, pfa7, amp_factor):
    """ Calculate the spa tmplt phase
    """
    err_msg = "This function is a stub that should be overridden using the "
    err_msg += "scheme. You shouldn't be seeing this error!"
    raise ValueError(err_msg)


def spa_tmplt(**kwds):
    """ Generate a minimal TaylorF2 approximant with optimizations for the sin/cos
    """
    distance = kwds['distance']
    mass1 = kwds['mass1']
    mass2 = kwds['mass2']
    s1z = kwds['spin1z']
    s2z = kwds['spin2z']
    phase_order = int(kwds['phase_order'])
    #amplitude_order = int(kwds['amplitude_order'])
    spin_order = int(kwds['spin_order'])

    if 'out' in kwds:
        out = kwds['out']
    else:
        out = None

    amp_factor = spa_amplitude_factor(mass1=mass1, mass2=mass2) / distance

    lal_pars = lal.CreateDict()
    if phase_order != -1:
        lalsimulation.SimInspiralWaveformParamsInsertPNPhaseOrder(
            lal_pars, phase_order)

    if spin_order != -1:
        lalsimulation.SimInspiralWaveformParamsInsertPNSpinOrder(
            lal_pars, spin_order)

    # Calculate the PN terms
    phasing = lalsimulation.SimInspiralTaylorF2AlignedPhasing(
                                    float(mass1), float(mass2),
                                    float(s1z), float(s2z),
                                    lal_pars)

    pfaN = phasing.v[0]
    pfa2 = phasing.v[2] / pfaN
    pfa3 = phasing.v[3] / pfaN
    pfa4 = phasing.v[4] / pfaN
    pfa5 = phasing.v[5] / pfaN
    pfa6 = (phasing.v[6] - phasing.vlogv[6] * log(4)) / pfaN
    pfa7 = phasing.v[7] / pfaN

    pfl5 = phasing.vlogv[5] / pfaN
    pfl6 = phasing.vlogv[6] / pfaN

    piM = PI * (mass1 + mass2) * MTSUN_SI

    if 'sample_points' not in kwds:
        f_lower = kwds['f_lower']
        delta_f = kwds['delta_f']
        kmin = int(f_lower / float(delta_f))

        # Get max frequency one way or another
        # f_final is assigned default value 0 in parameters.py
        if 'f_final' in kwds and kwds['f_final'] > 0.:
            fstop = kwds['f_final']
        elif 'f_upper' in kwds:
            fstop = kwds['f_upper']
            warnings.warn('f_upper is deprecated in favour of f_final!',
                          DeprecationWarning)
        else:
            # Schwarzschild ISCO frequency
            vISCO = 1. / sqrt(6.)
            fstop = vISCO * vISCO * vISCO / piM
        if fstop <= f_lower:
            raise ValueError("cannot generate waveform! f_lower >= f_final"
                             f" ({f_lower}, {fstop})")

        kmax = int(fstop / delta_f)
        f_max = ceilpow2(fstop)
        n = int(f_max / delta_f) + 1

        if not out:
            htilde = FrequencySeries(zeros(n, dtype=numpy.complex64), delta_f=delta_f, copy=False)
        else:
            if type(out) is not Array:
                raise TypeError("Output must be an instance of Array")
            if len(out) < kmax:
                kmax = len(out)
            if out.dtype != complex64:
                raise TypeError("Output array is the wrong dtype")
            htilde = FrequencySeries(out, delta_f=delta_f, copy=False)

        spa_tmplt_engine(htilde[kmin:kmax], kmin, phase_order,
                         delta_f, piM, pfaN,
                         pfa2, pfa3, pfa4, pfa5, pfl5,
                         pfa6, pfl6, pfa7, amp_factor)
    else:
        from .spa_tmplt_cpu import spa_tmplt_inline_sequence
        htilde = numpy.empty(len(kwds['sample_points']), dtype=numpy.complex64)
        spa_tmplt_inline_sequence(
            piM, pfaN, pfa2, pfa3, pfa4, pfa5, pfl5, pfa6, pfl6, pfa7,
            amp_factor, kwds['sample_points'], htilde)

    return htilde
