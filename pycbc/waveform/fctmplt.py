# Copyright (C) 2013  Alex Nitz
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
"""This module contains a legacy wrapper to use the FindChirpSPTemplate approximant
"""
import lalsimulation
import numpy
import lal
import pycbc
import pycbc.pnutils
from math import sqrt
from pycbc.types import Array, zeros, complex64, float32, FrequencySeries, complex128
from pycbc.waveform.spa_tmplt import spa_tmplt_precondition, spa_amplitude_factor, spa_tmplt

def findchirp_template(**p):
    import lalinspiral
    m1 = p['mass1']
    m2 = p['mass2']
    s1z = p['spin1z']
    s2z = p['spin2z']

    M = m1 + m2
    mc, et = pycbc.pnutils.mass1_mass2_to_mchirp_eta(m1, m2)

    m_sec = M * lal.MTSUN_SI;
    piM = lal.PI * m_sec;
    kmin = int(p['f_lower'] / float(p['delta_f']))
    vISCO = 1. / sqrt(6.)
    fISCO = vISCO * vISCO * vISCO / piM;
    kmax = int(fISCO / p['delta_f'])

    if 'out' in p:
        n = len(p['out'])
    else:
        n = kmax

    N = (n-1)*2
    delta_t = 1.0 / (p['delta_f'] * N)

    amp_factor = spa_amplitude_factor(mass1=m1, mass2=m2) / p['distance']

    fctmplt = lalinspiral.FindChirpTemplate()
    fctmplt.data = zeros(n, dtype=complex64).lal()

    tmplt = lalinspiral.InspiralTemplate()
    tmplt.totalMass = M
    tmplt.eta = et
    tmplt.mu = m1*m2 / M
    tmplt.spin1[2] = s1z
    tmplt.spin2[2] = s2z
    tmplt.fFinal = fISCO

    params = lalinspiral.FindChirpTmpltParams()

    vec = numpy.arange(0, n, 1)
    vec2 = numpy.zeros(n, dtype=float32)
    vec2[1:] = vec[1:] ** (-1.0/3.0)
    params.xfacVec = Array(vec2, dtype=float32).lal()

    params.approximant=lalsimulation.FindChirpSP

    # Max implemented order is 7 for fctmplt
    if p['phase_order'] == -1:
        p['phase_order'] = 7

    params.order = int(p['phase_order'])
    params.deltaT = delta_t
    params.fLow = p['f_lower']
    params.dynRange = pycbc.DYN_RANGE_FAC

    lalinspiral.FindChirpSPTemplate(fctmplt, tmplt, params)
    kfac = spa_tmplt_precondition(n, p['delta_f'])

    htilde = FrequencySeries(fctmplt.data.data, delta_f=p['delta_f'])
    htilde *= (amp_factor * kfac)

    if 'out' in p:
        p['out'][:] = htilde[:]
        htilde = FrequencySeries(p['out'], copy=False, delta_f=p['delta_f'])

    return htilde
