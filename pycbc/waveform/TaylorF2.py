#  Adapted from code in LALSimInspiralTaylorF2.c 
# 
#  Copyright (C) 2007 Jolien Creighton, B.S. Sathyaprakash, Thomas Cokelaer
#  Copyright (C) 2012 Leo Singer, Alex Nitz
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


import lal
import numpy
from numpy import sqrt

from pycuda.elementwise import ElementwiseKernel
from pycbc.libutils import pkg_config_header_strings
from pycbc.types import FrequencySeries, zeros
import pycbc.pnutils

preamble = """
#include "stdio.h"
#include <lal/LALConstants.h>
"""
 
taylorf2_text = """
    const double f = (i + kmin ) * delta_f;
    const double v0 = cbrt(piM *  kmin * delta_f);
    const double v = cbrt(piM*f);
    const double v2 = v * v;
    const double v3 = v * v2;
    const double v4 = v * v3;
    const double v5 = v * v4;
    const double v6 = v * v5;
    const double v7 = v * v6;
    const double v8 = v * v7;
    const double v9 = v * v8;
    const double v10 = v * v9;
    double phasing = 0.;
    double dEnergy = 0.;
    double flux = 0.;
    double amp;
    double shft = -LAL_TWOPI * tC;

    switch (phase_order)
    {
        case -1:
        case 7:
            phasing += pfa7 * v7;
        case 6:
            phasing += (pfa6 + pfl6 * log(4.*v) ) * v6;
        case 5:
            phasing += (pfa5 + pfl5 * log(v/v0)) * v5;
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
            flux += ( FTa6 +  FTl6*log(16.*v2)) * v6;
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

    phasing += shft * f + phi0;
    amp = amp0 * sqrt(-dEnergy/flux) * v;

    htilde[i]._M_re = amp * cos(phasing - LAL_PI_4) ;
    htilde[i]._M_im = - amp * sin(phasing - LAL_PI_4);

"""

taylorf2_kernel = ElementwiseKernel("""pycuda::complex<double> *htilde, int kmin, int phase_order,
                                       int amplitude_order, double delta_f, double piM, double pfaN, 
                                       double pfa2, double pfa3, double pfa4, double pfa5, double pfl5,
                                       double pfa6, double pfl6, double pfa7, double FTaN, double FTa2, 
                                       double FTa3, double FTa4, double FTa5, double FTa6,
                                       double FTl6, double FTa7, double dETaN, double dETa1, double dETa2, double dETa3,
                                       double amp0, double tC, double phi0""",
                    taylorf2_text, "taylorf2_kernel",
                    preamble=preamble, options=pkg_config_header_strings(['lal']))

def taylorf2(**kwds):
    """ Return a TaylorF2 waveform using CUDA to generate the phase and amplitude
    """
    # Pull out the input arguments
    delta_f = kwds['delta_f']
    distance = kwds['distance']
    mass1 = kwds['mass1']
    mass2 = kwds['mass2']
    phase_order = int(kwds['phase_order'])
    amplitude_order = int(kwds['amplitude_order'])
    phi0 = kwds['phi0']

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
    pfa3 = -16.0*lal.PI + 4.0*beta;
    pfa4 = 5.0*(3058.673/7.056 + 5429.0/7.0 * eta + 617.0 * eta*eta)/72.0 - \
            10.0*sigma
    pfa5 = 5.0/9.0 * (7729.0/84.0 - 13.0 * eta) * lal.PI - gamma
    pfl5 = 5.0/3.0 * (7729.0/84.0 - 13.0 * eta) * lal.PI - gamma * 3
    pfa6 = (11583.231236531/4.694215680 - 640.0/3.0 * lal.PI * lal.PI- \
            6848.0/21.0*lal.GAMMA) + \
            eta * (-15335.597827/3.048192 + 2255./12. * lal.PI * \
            lal.PI - 1760./3.*theta +12320./9.*lambdaa) + \
            eta*eta * 76055.0/1728.0 - \
            eta*eta*eta*  127825.0/1296.0 
    pfl6 = -6848.0/21.0;
    pfa7 = lal.PI * 5.0/756.0 * ( 15419335.0/336.0 + 75703.0/2.0 * eta - \
            14809.0 * eta*eta)

    FTaN = 32.0 * eta*eta / 5.0;
    FTa2 = -(12.47/3.36 + 3.5/1.2 * eta)
    FTa3 = 4.0 * lal.PI
    FTa4 = -(44.711/9.072 - 92.71/5.04 * eta - 6.5/1.8 * eta*eta)
    FTa5 = -(81.91/6.72 + 58.3/2.4 * eta) * lal.PI
    FTa6 = (664.3739519/6.9854400 + 16.0/3.0 * lal.PI*lal.PI -
            17.12/1.05 * lal.GAMMA +
            (4.1/4.8 * lal.PI*lal.PI - 134.543/7.776) * eta -
            94.403/3.024 * eta*eta - 7.75/3.24 * eta*eta*eta)
    FTl6 = -8.56/1.05
    FTa7 = -(162.85/5.04 - 214.745/1.728 * eta - 193.385/3.024 * eta*eta) \
            * lal.PI

    dETaN = 2 * -eta/2.0;
    dETa1 = 2 * -(3.0/4.0 + 1.0/12.0 * eta)
    dETa2 = 3 * -(27.0/8.0 - 19.0/8.0 * eta + 1./24.0 * eta*eta)
    dETa3 = 4 * -(67.5/6.4 - (344.45/5.76 - 20.5/9.6 * lal.PI*lal.PI) *
                             eta + 15.5/9.6 * eta*eta + 3.5/518.4 * eta*eta*eta)
  
    amp0 = -4. * mass1 * mass2 / (1.0e+06 * float(distance) * lal.PC_SI )* \
                    lal.MRSUN_SI * lal.MTSUN_SI * sqrt(lal.PI/12.0)
    
    m_sec = M * lal.MTSUN_SI;
    piM = lal.PI * m_sec;

    kmin = int(kwds['f_lower'] / float(delta_f))

    vISCO = 1. / sqrt(6.)
    fISCO = vISCO * vISCO * vISCO / piM;
    kmax = int(fISCO / delta_f)
    f_max = fISCO
    n = int(f_max / delta_f) + 1;

    htilde = FrequencySeries(zeros(n,dtype=numpy.complex128), delta_f=delta_f, copy=False)
    taylorf2_kernel(htilde.data[kmin:kmax],  kmin,  phase_order,
                    amplitude_order,  delta_f,  piM,  pfaN, 
                    pfa2,  pfa3,  pfa4,  pfa5,  pfl5,
                    pfa6,  pfl6,  pfa7,  FTaN,  FTa2, 
                    FTa3,  FTa4,  FTa5,  FTa6,
                    FTl6,  FTa7,  dETaN, dETa1, dETa2,  dETa3,
                    amp0,  tC,  phi0)
                    
    hp = htilde 
    hc = htilde * 1j
    return hp, hc
    


