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
#



import lalsimulation
import lal
import numpy
import pycuda.tools
from pycuda.elementwise import ElementwiseKernel
from pycbc.types import FrequencySeries,zeros
from pycuda.gpuarray import to_gpu
from numpy import sqrt
from math import frexp
pntype = numpy.dtype([('pfaN', numpy.float64),
                    ('pfa2', numpy.float64),
                    ('pfa3', numpy.float64),
                    ('pfa4', numpy.float64),
                    ('pfa5', numpy.float64),
                    ('pfl5', numpy.float64),
                    ('pfa6', numpy.float64),
                    ('pfl6', numpy.float64),
                    ('pfa7', numpy.float64),
                    ('delta_f',numpy.float64),
                    ('FTaN',numpy.float64),
                    ('FTa2',numpy.float64),
                    ('FTa3',numpy.float64),
                    ('FTa4',numpy.float64),
                    ('FTa5',numpy.float64),
                    ('FTa6',numpy.float64),
                    ('FTl6',numpy.float64),
                    ('FTa7',numpy.float64),
                    ('piM',numpy.float64),
                    ('amp0',numpy.float64),
                    ('dEtaN',numpy.float64),
                    ('dEta1',numpy.float64),
                    ('dEta2',numpy.float64),
                    ('dEta3',numpy.float64),
                    ('kmin',numpy.float64),
                    ('phi0',numpy.float64),
                    ('tC',numpy.float64),
                    ('phase0', numpy.int64),
                    ('amplitudeO', numpy.int64),
                     ])
pycuda.tools.register_dtype(pntype,"pntype")

def get_taylorf2_pn_coefficients(mass1,mass2,distance,beta =0 , sigma = 0):
    # FIXME when lalsimulation has all the coeffiecients wrapped    

    M = float(mass1) + float(mass2)
    eta = mass1 * mass2 / (M * M)
    theta = -11831./9240.;
    lambdaa = -1987./3080.0;
    pfaN = 3.0/(128.0 * eta);
    pfa2 = 5*(743.0/84 + 11.0 * eta)/9.0;
    pfa3 = -16.0*lal.LAL_PI + 4.0*beta;
    pfa4 = 5.0*(3058.673/7.056 + 5429.0/7.0 * eta + 617.0 * eta*eta)/72.0 - \
            10.0*sigma
    pfa5 = 5.0/9.0 * (7729.0/84.0 - 13.0 * eta) * lal.LAL_PI
    pfl5 = 5.0/3.0 * (7729.0/84.0 - 13.0 * eta) * lal.LAL_PI
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
  
    amp0 = 4. * mass1 * mass2 / (float(distance) * lal.LAL_PC_SI )* \
                    lal.LAL_MRSUN_SI * lal.LAL_MTSUN_SI * sqrt(lal.LAL_PI/12.0)    
    
    m_sec = M * lal.LAL_MTSUN_SI;
    piM = lal.LAL_PI * m_sec; 

    pn_coefficients = numpy.array((pfaN,pfa2,pfa3,pfa4,pfa5,pfl5,pfa6,pfl6,pfa7,
                                    0,FTaN,FTa2,FTa3,FTa4,FTa5,FTa6,FTl6,FTa7,
                                    piM,amp0,dETaN,dETa1,dETa2,dETa3,0,0,0,0,0),
                                    dtype=pntype)

    return pn_coefficients

preamble = """
#include <lal/LALConstants.h>

struct pntype{
double pfaN;
double pfa2;
double pfa3;
double pfa4;
double pfa5;
double pfl5;
double pfa6;
double pfl6;
double pfa7;
double delta_f;
double FTaN;
double FTa2;
double FTa3;
double FTa4;
double FTa5;
double FTa6;
double FTl6;
double FTa7;
double piM;
double amp0;
double dETaN;
double dETa1;
double dETa2;
double dETa3;
double kmin;
double phi0;
double tC;
long phaseO;
long amplitudeO;
};
"""
 
taylorf2_text = """

    const double f = (i + cf->kmin ) * cf->delta_f;
    const double v0 = cbrt(cf->piM * cf-> kmin * cf->delta_f);
    const double v = cbrt(cf->piM*f);
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
    double shft = -LAL_TWOPI * cf->tC;
    double phi0 = cf->phi0;

    switch (cf->phaseO)
    {
        case -1:
        case 7:
            phasing += cf->pfa7 * v7;
        case 6:
            phasing += (cf->pfa6 + cf->pfl6 * log(4.*v) ) * v6;
        case 5:
            phasing += (cf->pfa5 + cf->pfl5 * log(v/v0)) * v5;
        case 4:
            phasing += cf->pfa4 * v4;
        case 3:
            phasing += cf->pfa3 * v3;
        case 2:
            phasing += cf->pfa2 * v2;
        case 0:
            phasing += 1.;
            break;
        default:
            break;
    }
    switch (cf -> amplitudeO)
    {
        case -1:
        case 7:
            flux += cf -> FTa7 * v7;
        case 6:
            flux += (cf -> FTa6 + cf -> FTl6*log(16.*v2)) * v6;
            dEnergy += cf -> dETa3 * v6;
        case 5:
            flux += cf -> FTa5 * v5;
        case 4:
            flux += cf -> FTa4 * v4;
            dEnergy += cf -> dETa2 * v4;
        case 3:
            flux += cf -> FTa3 * v3;
        case 2:
            flux += cf -> FTa2 * v2;
            dEnergy += cf -> dETa1 * v2;
        case 0:
            flux += 1.;
            dEnergy += 1.;
            break;
    }

    flux+=1.0;
    dEnergy+=1.0;

    phasing *= cf->pfaN / v5;
    flux *= cf->FTaN * v10;
    dEnergy *= cf->dETaN * v;

    phasing += shft * f + phi0;
    amp = cf->amp0 * sqrt(-dEnergy/flux) * v;
    htilde[i]._M_re = amp * cos(phasing + LAL_PI_4) ;
    htilde[i]._M_im = - amp * sin(phasing + LAL_PI_4);

"""

taylorf2_kernel = ElementwiseKernel("pycuda::complex<double> *htilde, pntype *cf",
                    taylorf2_text, "taylorf2_kernel",preamble=preamble)

def ceilpow2(n):
    signif,exponent = frexp(n)
    if (signif < 0):
        return 1;
    if (signif == 0.5):
        exponent -= 1;
    return (1) << exponent;


def taylorf2(tC = None, beta =0, sigma = 0,**kwds):
    pn_const = get_taylorf2_pn_coefficients(kwds['mass1'],kwds['mass2'],
                                    kwds['distance'],beta,sigma)

    pn_const['phase0'] = int(kwds['phase_order'])
    pn_const['amplitudeO'] = int(kwds['phase_order'])
    delta_f = kwds['delta_f']
    if tC is None:
        tC = -1.0 / delta_f 
    pn_const['delta_f'] = delta_f
    pn_const['phi0']=kwds['phi0']
    pn_const['tC']= tC

    kmin = int(kwds['f_lower'] / float(delta_f))

    piM = pn_const['piM']
    vISCO = 1. / sqrt(6.)
    fISCO = vISCO * vISCO * vISCO / piM;
    kmax = int(fISCO / delta_f)
    f_max = ceilpow2(fISCO);
    n = int(f_max / delta_f) + 1;
    pn_const['kmin'] = kmin

    htilde = FrequencySeries(zeros(n,dtype=numpy.complex128),delta_f = delta_f, 
                                                                    copy=False)
    taylorf2_kernel(htilde.data[kmin:kmax],to_gpu(pn_const))
    return htilde
    


