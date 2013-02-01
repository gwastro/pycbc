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
    float dEnergy = 0.;
    float flux = 0.;
    float amp;
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
    switch (amplitude_order)
    {
        case -1:
        case 7:
            flux +=  FTa7 * v7;
        case 6:
            flux += ( FTa6 +  FTl6*__logf(16.*v2)) * v6;
            dEnergy +=  dETa3 * v6;
        case 5:
            flux +=  FTa5 * v5;
        case 4:
            flux +=  FTa4 * v4;
            dEnergy +=  dETa2 * v4;
        case 3:
            flux +=  FTa3 * v3;
        case 2:
            flux +=  FTa2 * v2;
            dEnergy +=  dETa1 * v2;
        case 0:
            flux += 1;
            dEnergy += 1.;
            break;
    }

    phasing *= pfaN / v5;
    flux *= FTaN * v10;
    dEnergy *= dETaN * v;
    phasing += shft * f + phi0 + LAL_PI_4;
    amp = amp0 * __powf(-dEnergy/flux, 0.5) * v;

    float pcos;
    float psin;
    __sincosf(phasing, &psin, &pcos);


    htilde[i]._M_re = amp * pcos;
    htilde[i]._M_im = - amp * psin;

"""

def ceilpow2(n):
    signif,exponent = frexp(n)
    if (signif < 0):
        return 1;
    if (signif == 0.5):
        exponent -= 1;
    return (1) << exponent;

taylorf2_kernel = ElementwiseKernel("""pycuda::complex<float> *htilde, int kmin, int phase_order,
                                       int amplitude_order, float delta_f, float piM, float pfaN, 
                                       float pfa2, float pfa3, float pfa4, float pfa5, float pfl5,
                                       float pfa6, float pfl6, float pfa7, float FTaN, float FTa2, 
                                       float FTa3, float FTa4, float FTa5, float FTa6,
                                       float FTl6, float FTa7, float dETaN, float dETa1, float dETa2, float dETa3,
                                       float amp0, float tC, float phi0, float v0""",
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

    #Calculate teh PN terms #TODO: replace with functions in lalsimulation!###
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

    FTaN = 32.0 * eta*eta / 5.0;
    FTa2 = -(12.47/3.36 + 3.5/1.2 * eta)
    FTa3 = 4.0 * lal.LAL_PI
    FTa4 = -(44.711/9.072 - 92.71/5.04 * eta - 6.5/1.8 * eta*eta)
    FTa5 = -(81.91/6.72 + 58.3/2.4 * eta) * lal.LAL_PI
    FTa6 = (664.3739519/6.9854400 + 16.0/3.0 * lal.LAL_PI*lal.LAL_PI - 
            17.12/1.05 * lal.LAL_GAMMA + 
		 (4.1/4.8 * lal.LAL_PI*lal.LAL_PI - 134.543/7.776) * eta -
		 94.403/3.024 * eta*eta - 7.75/3.24 * eta*eta*eta)
    FTl6 = -8.56/1.05
    FTa7 = -(162.85/5.04 - 214.745/1.728 * eta - 193.385/3.024 * eta*eta) \
            * lal.LAL_PI

    dETaN = 2 * -eta/2.0;
    dETa1 = 2 * -(3.0/4.0 + 1.0/12.0 * eta)
    dETa2 = 3 * -(27.0/8.0 - 19.0/8.0 * eta + 1./24.0 * eta*eta)
    dETa3 = 4 * -(67.5/6.4 - (344.45/5.76 - 20.5/9.6 * lal.LAL_PI*lal.LAL_PI) *
                             eta + 15.5/9.6 * eta*eta + 3.5/518.4 * eta*eta*eta)
  
    amp0 = 4. * mass1 * mass2 / (1.0e+03 * float(distance) * lal.LAL_PC_SI )* \
                    lal.LAL_MRSUN_SI * lal.LAL_MTSUN_SI * sqrt(lal.LAL_PI/12.0)    
    
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
                    amplitude_order,  delta_f,  piM,  pfaN, 
                    pfa2,  pfa3,  pfa4,  pfa5,  pfl5,
                    pfa6,  pfl6,  pfa7,  FTaN,  FTa2, 
                    FTa3,  FTa4,  FTa5,  FTa6,
                    FTl6,  FTa7,  dETaN, dETa1, dETa2,  dETa3,
                    amp0,  tC,  phi0, v0)
    return htilde
    


