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
import numpy, lal, lalsimulation, pycbc.pnutils
from pycbc.scheme import schemed
from pycbc.types import FrequencySeries, Array, complex64, float32, zeros
from pycbc.waveform.utils import ceilpow2

def findchirp_chirptime(m1, m2, fLower, porder):
    # variables used to compute chirp time
    m1 = float(m1)
    m2 = float(m2)
    m = m1 + m2
    eta = m1 * m2 / m / m
    c0T = c2T = c3T = c4T = c5T = c6T = c6LogT = c7T = 0.

    # All implemented option
    if porder == -1:
        porder = 7

    if porder >= 7:
        c7T = lal.PI * (14809.0 * eta * eta / 378.0 - 75703.0 * eta / 756.0 - 15419335.0 / 127008.0)

    if porder >= 6:
        c6T = lal.GAMMA * 6848.0 / 105.0 - 10052469856691.0 / 23471078400.0 +\
            lal.PI * lal.PI * 128.0 / 3.0 + \
            eta * (3147553127.0 / 3048192.0 - lal.PI * lal.PI * 451.0 / 12.0) -\
            eta * eta * 15211.0 / 1728.0 + eta * eta * eta * 25565.0 / 1296.0 +\
            eta * eta * eta * 25565.0 / 1296.0 + numpy.log(4.0) * 6848.0 / 105.0
        c6LogT = 6848.0 / 105.0

    if porder >= 5:
        c5T = 13.0 * lal.PI * eta / 3.0 - 7729.0 * lal.PI / 252.0

    if porder >= 4:
        c4T = 3058673.0 / 508032.0 + eta * (5429.0 / 504.0 + eta * 617.0 / 72.0)
        c3T = -32.0 * lal.PI / 5.0
        c2T = 743.0 / 252.0 + eta * 11.0 / 3.0
        c0T = 5.0 * m * lal.MTSUN_SI / (256.0 * eta)

    # This is the PN parameter v evaluated at the lower freq. cutoff
    xT = pow (lal.PI * m * lal.MTSUN_SI * fLower, 1.0 / 3.0)
    x2T = xT * xT
    x3T = xT * x2T
    x4T = x2T * x2T
    x5T = x2T * x3T
    x6T = x3T * x3T
    x7T = x3T * x4T
    x8T = x4T * x4T

    # Computes the chirp time as tC = t(v_low)
    # tC = t(v_low) - t(v_upper) would be more
    # correct, but the difference is negligble.

    # This formula works for any PN order, because
    # higher order coeffs will be set to zero.
    return c0T * (1 + c2T * x2T + c3T * x3T + c4T * x4T + c5T * x5T + (c6T + c6LogT * numpy.log(xT)) * x6T + c7T * x7T) / x8T

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

    # For now, we call the swig-wrapped function below in
    # lalinspiral.  Eventually would be nice to replace this
    # with a function using PN coeffs from lalsimulation.
    return findchirp_chirptime(m1, m2, flow, porder)

def spa_amplitude_factor(**kwds):
    m1 = kwds['mass1']
    m2 = kwds['mass2']

    _, eta = pycbc.pnutils.mass1_mass2_to_mchirp_eta(m1, m2)

    FTaN = 32.0 * eta*eta / 5.0
    dETaN = 2 * -eta/2.0

    M = m1 + m2

    m_sec = M * lal.MTSUN_SI
    piM = lal.PI * m_sec

    amp0 = 4. * m1 * m2 / (1e6 * lal.PC_SI ) * lal.MRSUN_SI * lal.MTSUN_SI * sqrt(lal.PI/12.0)

    fac = numpy.sqrt( -dETaN / FTaN) * amp0 * (piM ** (-7.0/6.0))
    return -fac

_prec = None
def spa_tmplt_precondition(length, delta_f, kmin=0):
    """Return the amplitude portion of the TaylorF2 approximant, used to precondition
    the strain data. The result is cached, and so should not be modified only read.
    """
    global _prec
    if _prec is None or _prec.delta_f != delta_f or len(_prec) < length:
        v = numpy.arange(0, (kmin+length*2), 1.0) * delta_f
        v = numpy.power(v[1:len(v)], -7.0/6.0)
        _prec = FrequencySeries(v, delta_f=delta_f, dtype=float32)
    return _prec[kmin:kmin + length]

def spa_tmplt_norm(psd, length, delta_f, f_lower):
    amp = spa_tmplt_precondition(length, delta_f)
    k_min = int(f_lower / delta_f)
    sigma = (amp[k_min:length].numpy() ** 2.0 / psd[k_min:length].numpy())
    norm_vec = numpy.zeros(length)
    norm_vec[k_min:length] = sigma.cumsum() * 4 * delta_f
    return norm_vec

def spa_tmplt_end(**kwds):
    return pycbc.pnutils.f_SchwarzISCO(kwds['mass1']+kwds['mass2'])

def spa_distance(psd, mass1, mass2, lower_frequency_cutoff, snr=8):
    """ Return the distance at a given snr (default=8) of the SPA TaylorF2
    template.
    """
    kend = int(spa_tmplt_end(mass1=mass1, mass2=mass2) / psd.delta_f)
    norm1 = spa_tmplt_norm(psd, len(psd), psd.delta_f, lower_frequency_cutoff)
    norm2 = (spa_amplitude_factor(mass1=mass1, mass2=mass2)) ** 2.0

    if kend >= len(psd):
        kend = len(psd) - 1
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
    """ Generate a minimal TaylorF2 approximant with optimations for the sin/cos
    """
    # Pull out the input arguments
    f_lower = kwds['f_lower']
    delta_f = kwds['delta_f']
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

    #Calculate the PN terms
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

    piM = lal.PI * (mass1 + mass2) * lal.MTSUN_SI

    kmin = int(f_lower / float(delta_f))

    vISCO = 1. / sqrt(6.)
    fISCO = vISCO * vISCO * vISCO / piM
    kmax = int(fISCO / delta_f)
    f_max = ceilpow2(fISCO)
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

    spa_tmplt_engine(htilde[kmin:kmax], kmin, phase_order, delta_f, piM, pfaN,
                     pfa2, pfa3, pfa4, pfa5, pfl5,
                     pfa6, pfl6, pfa7, amp_factor)
    return htilde

