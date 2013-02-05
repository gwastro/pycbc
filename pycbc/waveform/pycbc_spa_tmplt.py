#  
#  Apapted from code in LALSimInpspiralTaylorF2.c 
# 
#  Copyright (C) 2007 Jolien Creighton, B.S. Sathyaprakash, Thomas Cokelaer
#  Copyright (C) 2012 Leo Singer, Alex Nitz
#  Adapted from code found in:
#    - LALSimInspiralTaylorF2.c
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with with program; see the file COPYING. If not, write to the
#  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#  MA  02111-1307  USA


import lalsimulation
import lal
import numpy
from numpy import sqrt
from math import frexp

import pycuda.tools
from pycuda.elementwise import ElementwiseKernel
from pycuda.gpuarray import to_gpu

from pycbc.setuputils import pkg_config_header_strings
from pycbc.types import FrequencySeries, zeros, Array, complex64
import pycbc.pnutils
from pycbc.filter import sigmasq_series

preamble = """
#include "stdio.h"
#include <lal/LALConstants.h>
"""
 
taylorf2_text = """
    const float f = (i + kmin ) * delta_f;
    const float v =  __powf(piM*f, 1.0/3.0);
    const float v2 = v * v;
    const float v3 = v * v * v;
    const float v4 = v2 * v2;
    const float v5 = v2 * v3;
    const float v6 = v3 * v3;
    const float v7 = v3 * v4;
    const float v10 = v5 * v5;
    float phasing = 0.;
    float shft = -LAL_TWOPI * tC;

    switch (phase_order)
    {
        case -1:
        case 7:
            phasing += pfa7 * v7;
        case 6:
            phasing += (pfa6 + pfl6 * __logf(4.*v) ) * v6;
        case 5:
            phasing += (pfa5 + pfl5 * __logf(v/v0)) * v5;
        case 4:
            phasing += pfa4 * v4;
        case 3:
            phasing += pfa3 * v3;
        case 2:
            phasing += pfa2 * v2;
        case 0:
            phasing += 1.;
            break;
        default:
            break;
    }

    phasing *= pfaN / v5;
    phasing += shft * f + phi0 + LAL_PI_4;

    float pcos;
    float psin;
    __sincosf(phasing, &psin, &pcos);


    htilde[i]._M_re = pcos;
    htilde[i]._M_im = - psin;

"""

def ceilpow2(n):
    signif,exponent = frexp(n)
    if (signif < 0):
        return 1;
    if (signif == 0.5):
        exponent -= 1;
    return (1) << exponent;

taylorf2_kernel = ElementwiseKernel("""pycuda::complex<float> *htilde, int kmin, int phase_order,
                                       float delta_f, float piM, float pfaN, 
                                       float pfa2, float pfa3, float pfa4, float pfa5, float pfl5,
                                       float pfa6, float pfl6, float pfa7, float tC, float phi0, float v0""",
                    taylorf2_text, "SPAtmplt",
                    preamble=preamble, options=pkg_config_header_strings(['lal']))

def spa_tmplt(**kwds):
    """ Return a TaylorF2 waveform using CUDA to generate the phase and amplitude
    """
    # Pull out the input arguments
    f_lower = kwds['f_lower']
    delta_f = kwds['delta_f']
    distance = kwds['distance']
    mass1 = kwds['mass1']
    mass2 = kwds['mass2']
    phase_order = int(kwds['phase_order'])
    amplitude_order = int(kwds['amplitude_order'])
    phi0 = kwds['phi0']

    if 'out' in kwds:
        out = kwds['out']
    else:
        out = None

    tC= -1.0 / delta_f 

    #Calculate the spin corrections
    beta, sigma, gamma = pycbc.pnutils.mass1_mass2_spin1z_spin2z_to_beta_sigma_gamma(
                                    mass1, mass2, kwds['spin1z'], kwds['spin2z'])

    #Calculate the PN terms #TODO: replace with functions in lalsimulation!###
    M = float(mass1) + float(mass2)
    eta = mass1 * mass2 / (M * M)
    theta = -11831./9240.;
    lambdaa = -1987./3080.0;
    pfaN = 3.0/(128.0 * eta);
    pfa2 = 5*(743.0/84 + 11.0 * eta)/9.0;
    pfa3 = -16.0*lal.LAL_PI + 4.0*beta;
    pfa4 = 5.0*(3058.673/7.056 + 5429.0/7.0 * eta + 617.0 * eta*eta)/72.0 - \
            10.0*sigma
    pfa5 = 5.0/9.0 * (7729.0/84.0 - 13.0 * eta) * lal.LAL_PI - gamma
    pfl5 = 5.0/3.0 * (7729.0/84.0 - 13.0 * eta) * lal.LAL_PI - gamma * 3
    pfa6 = (11583.231236531/4.694215680 - 640.0/3.0 * lal.LAL_PI * lal.LAL_PI- \
            6848.0/21.0*lal.LAL_GAMMA) + \
            eta * (-15335.597827/3.048192 + 2255./12. * lal.LAL_PI * \
            lal.LAL_PI - 1760./3.*theta +12320./9.*lambdaa) + \
            eta*eta * 76055.0/1728.0 - \
            eta*eta*eta*  127825.0/1296.0 
    pfl6 = -6848.0/21.0;
    pfa7 = lal.LAL_PI * 5.0/756.0 * ( 15419335.0/336.0 + 75703.0/2.0 * eta - \
            14809.0 * eta*eta)
    
    m_sec = M * lal.LAL_MTSUN_SI;
    piM = lal.LAL_PI * m_sec; 

    kmin = int(f_lower / float(delta_f))

    vISCO = 1. / sqrt(6.)
    fISCO = vISCO * vISCO * vISCO / piM;
    kmax = int(fISCO / delta_f)
    f_max = ceilpow2(fISCO);
    n = int(f_max / delta_f) + 1;

    v0 = (piM *  kmin * delta_f) ** (1.0/3.0)

    if not out:
        htilde = FrequencySeries(zeros(n,dtype=numpy.complex64), delta_f=delta_f, copy=False)
    else:
        if type(out) is not Array:
            raise TypeError("Output must be an instance of Array")
        if len(out) < kmax:
            raise TypeError("Output array is too small")
        if out.dtype != complex64:
            raise TypeError("Output array is the wrong dtype")
        htilde = FrequencySeries(out, delta_f=delta_f, copy=False)
        htilde.clear()
    
    taylorf2_kernel(htilde.data[kmin:kmax],  kmin,  phase_order,
                    delta_f,  piM,  pfaN, 
                    pfa2,  pfa3,  pfa4,  pfa5,  pfl5,
                    pfa6,  pfl6,  pfa7, tC,  phi0, v0)
    return htilde
    
def spa_tmplt_precondition(length, delta_f):
    """Return the amplitude portion of the TaylorF2 approximant, used to precondition
    the strian data.
    """
    v = numpy.arange(0, length+1, 1) * delta_f
    v = v[1:len(v)]**(-7.0/6.0)
    return FrequencySeries(v, delta_f=delta_f)
    
def spa_tmplt_norm(psd, length, delta_f, f_lower):
    amp = spa_tmplt_precondition(length, delta_f)
    k_min = int(f_lower / delta_f)
    sigma = (amp.numpy() ** 2.0 / psd.numpy())
    norm_vec = sigma * 0 
    norm_vec[k_min:length] = sigma[k_min:length].cumsum()
    return norm_vec

def spa_tmplt_end(**kwds):
    return pycbc.pnutils.schwarzschild_isco(kwds['mass1']+kwds['mass2'])
    


