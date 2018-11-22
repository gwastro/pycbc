#  Copyright (C) 2013 Haris K
#  Ported from LALSimulation's LALSimInspiralSpinTaylorF2.c
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
from numpy import sqrt, double, complex128
from math import pow, log, cos, sin, acos, atan2

from pycuda.elementwise import ElementwiseKernel

from pycbc.libutils import pkg_config_header_strings
from pycbc.types import FrequencySeries, zeros
from pycbc.waveform.utils import ceilpow2

preamble = """
#include <lal/LALConstants.h>
#include <cuComplex.h>
"""

spintaylorf2_text = """
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
    double RE_prec_facP;
    double IM_prec_facP;
    double RE_prec_facC;
    double IM_prec_facC;

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

    const double gam = gamma0*v;
    const double sqrtfac = sqrt(1. + 2.*kappa*gam + gam*gam);
    const double logv = log(v);
    const double logfac1 = log(1. + kappa*gam + sqrtfac);
    const double logfac2 = log(kappa + gam + sqrtfac);
    const double kappa2 = kappa * kappa;
    const double kappa3 = kappa2 * kappa;
    const double gamma02 = gamma0 * gamma0;
    const double gamma03 = gamma02 *gamma0;

    const double alpha = prec_fac0*(  logfac2 *( dtdv2*gamma0 + dtdv3*kappa -
        dtdv5*kappa/(2.*gamma02) + dtdv4/(2.*gamma0) -
        dtdv4*kappa2/(2.*gamma0) + (dtdv5*kappa3)/(2.*gamma02) )  +
        logfac1*( - dtdv2*gamma0*kappa - dtdv3 + kappa*gamma03/2. -
        gamma03*kappa3/2. ) + logv *( dtdv2*gamma0*kappa + dtdv3 -
        kappa*gamma03/2. + gamma03*kappa3/2. ) + sqrtfac *( dtdv3 +
        dtdv4*v/2. + dtdv5/gamma02/3. + dtdv4*kappa/(2.*gamma0) +
        dtdv5*kappa*v/(6.*gamma0) - dtdv5*kappa2/(2.*gamma02) - 1/(3.*v3) -
        gamma0*kappa/(6.*v2) - dtdv2/v - gamma02/(3.*v) +
        gamma02*kappa2/(2.*v) + dtdv5*v2/3. )) - alpha_ref;

    const double beta = acos((1. + kappa*gamma0*v)/sqrt(1. + 2.*kappa*gamma0*v + gamma0*gamma0*v*v));

    const double zeta = prec_fac0*( dtdv3*gamma0*kappa*v + dtdv4*v +
        logfac2 *(-dtdv2*gamma0 - dtdv3*kappa + dtdv5*kappa/(2.*gamma02) -
        dtdv4/(2.*gamma0) + dtdv4*kappa2/(2.*gamma0) -
        dtdv5*kappa3/(2.*gamma02) ) + logv *( kappa*gamma03/2. -
        gamma03*kappa3/2. ) + logfac1 *( dtdv2*gamma0*kappa + dtdv3 -
        kappa*gamma03/2. + gamma03*kappa3/2. ) - 1/(3.*v3) -
        gamma0*kappa/(2.*v2) - dtdv2/v + dtdv4*gamma0*kappa*v2/2. +
        dtdv5*v2/2. + sqrtfac *( -dtdv3 - dtdv4*v/2. - dtdv5/(3.*gamma02) -
        dtdv4*kappa/(2.*gamma0) - dtdv5*kappa*v/(6.*gamma0) +
        dtdv5*kappa2/(2.*gamma02) + 1/(3.*v3) + gamma0*kappa/(6.*v2) +
        dtdv2/v + gamma02/(3.*v) - gamma02*kappa2/(2.*v) - dtdv5*v2/3. ) +
        dtdv5*gamma0*kappa*v3/3. ) - zeta_ref;

    double CBeta;
    double SBeta;
    double SAlpha1;
    double SAlpha2;
    double SAlpha3;
    double SAlpha4;
    double CAlpha1;
    double CAlpha2;
    double CAlpha3;
    double CAlpha4;
    sincos(beta/2.,&SBeta,&CBeta);
    sincos(-alpha,&SAlpha1,&CAlpha1);
    sincos(-2.*alpha,&SAlpha2,&CAlpha2);
    sincos(-3.*alpha,&SAlpha3,&CAlpha3);
    sincos(-4.*alpha,&SAlpha4,&CAlpha4);
    const double CBeta2 = CBeta * CBeta;
    const double CBeta3 = CBeta * CBeta2;
    const double CBeta4 = CBeta * CBeta3;
    const double SBeta2 = SBeta * SBeta;
    const double SBeta3 = SBeta * SBeta2;
    const double SBeta4 = SBeta * SBeta3;

    RE_prec_facP = ( cos(2.*psiJ_P) *
            ( SBeta4 * RE_SBfac4 * CAlpha4
            + CBeta * SBeta3 * RE_SBfac3 * CAlpha3
            + CBeta2 * SBeta2 * RE_SBfac2 * CAlpha2
            + CBeta3 * SBeta * RE_SBfac1 * CAlpha1
            + CBeta4 * RE_SBfac0 )
                 - sin(2.*psiJ_P) *
            ( SBeta4 * IM_SBfac4 * SAlpha4
            + CBeta * SBeta3 * IM_SBfac3 * SAlpha3
            + CBeta2 * SBeta2 * IM_SBfac2 * SAlpha2
            + CBeta3 * SBeta * IM_SBfac1 * SAlpha1
            + CBeta4 * IM_SBfac0 * 0 ));

    IM_prec_facP = ( cos(2.*psiJ_P) *
            (  SBeta4 * RE_SBfac4 * SAlpha4
            + CBeta * SBeta3 * RE_SBfac3 * SAlpha3
            + CBeta2 * SBeta2 * RE_SBfac2 * SAlpha2
            + CBeta3 * SBeta * RE_SBfac1 * SAlpha1
            + CBeta4 * RE_SBfac0 * 0 )
                  + sin(2.*psiJ_P) *
            ( SBeta4 * IM_SBfac4 * CAlpha4
            + CBeta * SBeta3 * IM_SBfac3 * CAlpha3
            + CBeta2 * SBeta2 * IM_SBfac2 * CAlpha2
            + CBeta3 * SBeta * IM_SBfac1 * CAlpha1
            + CBeta4 * IM_SBfac0 ));

    RE_prec_facC = ( cos(2.*psiJ_C) *
            ( SBeta4 * RE_SBfac4 * CAlpha4
            + CBeta * SBeta3 * RE_SBfac3 * CAlpha3
            + CBeta2 * SBeta2 * RE_SBfac2 * CAlpha2
            + CBeta3 * SBeta * RE_SBfac1 * CAlpha1
            + CBeta4 * RE_SBfac0 )
                 - sin(2.*psiJ_C) *
            ( SBeta4 * IM_SBfac4 * SAlpha4
            + CBeta * SBeta3 * IM_SBfac3 * SAlpha3
            + CBeta2 * SBeta2 * IM_SBfac2 * SAlpha2
            + CBeta3 * SBeta * IM_SBfac1 * SAlpha1
            + CBeta4 * IM_SBfac0 * 0 ));

    IM_prec_facC = ( cos(2.*psiJ_C) *
            (  SBeta4 * RE_SBfac4 * SAlpha4
            + CBeta * SBeta3 * RE_SBfac3 * SAlpha3
            + CBeta2 * SBeta2 * RE_SBfac2 * SAlpha2
            + CBeta3 * SBeta * RE_SBfac1 * SAlpha1
            + CBeta4 * RE_SBfac0 * 0 )
                  + sin(2.*psiJ_C) *
            ( SBeta4 * IM_SBfac4 * CAlpha4
            + CBeta * SBeta3 * IM_SBfac3 * CAlpha3
            + CBeta2 * SBeta2 * IM_SBfac2 * CAlpha2
            + CBeta3 * SBeta * IM_SBfac1 * CAlpha1
            + CBeta4 * IM_SBfac0 ));


    phasing += shft * f - 2. * phi0; // FIXME:: Sign of phi0?
    phasing += 2.*zeta;
    amp = amp0 * sqrt(-dEnergy/flux) * v;

    const double CPhasing = amp * cos(phasing - LAL_PI_4);
    const double SPhasing = amp * sin(phasing - LAL_PI_4);
    htildeP[i]._M_re = RE_prec_facP * CPhasing + IM_prec_facP * SPhasing ;
    htildeP[i]._M_im = IM_prec_facP * CPhasing - RE_prec_facP * SPhasing ;
    htildeC[i]._M_re = RE_prec_facC * CPhasing + IM_prec_facC * SPhasing ;
    htildeC[i]._M_im = IM_prec_facC * CPhasing - RE_prec_facC * SPhasing ;

"""

spintaylorf2_kernel = ElementwiseKernel("""pycuda::complex<double> *htildeP,
                                           pycuda::complex<double> *htildeC,
                                           int kmin, int phase_order,
                                           int amplitude_order, double delta_f,
                                           double piM, double pfaN,
                                           double pfa2, double pfa3,
                                           double pfa4, double pfa5,
                                           double pfl5, double pfa6,
                                           double pfl6, double pfa7,
                                           double FTaN, double FTa2,
                                           double FTa3, double FTa4,
                                           double FTa5, double FTa6,
                                           double FTl6, double FTa7,
                                           double dETaN, double dETa1,
                                           double dETa2, double dETa3,
                                           double amp0, double tC, double phi0,
                                           double kappa, double prec_fac0,
                                           double alpha_ref, double zeta_ref,
                                           double dtdv2, double dtdv3,
                                           double dtdv4, double dtdv5,
                                           double RE_SBfac0, double RE_SBfac1,
                                           double RE_SBfac2, double RE_SBfac3,
                                           double RE_SBfac4, double IM_SBfac0,
                                           double IM_SBfac1, double IM_SBfac2,
                                           double IM_SBfac3, double IM_SBfac4,
                                           double psiJ_P, double psiJ_C,
                                           double gamma0""",
                    spintaylorf2_text, "spintaylorf2_kernel",
                    preamble=preamble, options=pkg_config_header_strings(['lal']))

def spintaylorf2(**kwds):
    """ Return a SpinTaylorF2 waveform using CUDA to generate the phase and amplitude
    """
    #####Pull out the input arguments#####
    f_lower = double(kwds['f_lower'])
    delta_f = double(kwds['delta_f'])
    distance = double(kwds['distance'])
    mass1 = double(kwds['mass1'])
    mass2 = double(kwds['mass2'])
    spin1x = double(kwds['spin1x'])
    spin1y = double(kwds['spin1y'])
    spin1z = double(kwds['spin1z'])
    phi0 = double(kwds['coa_phase'])               #Orbital Phase at coalescence
    phase_order = int(kwds['phase_order'])
    amplitude_order = int(kwds['amplitude_order'])
    inclination = double(kwds['inclination'])
    lnhatx = sin(inclination)
    lnhaty = 0.
    lnhatz = cos(inclination)
    psi = 0.

    tC= -1.0 / delta_f
    M = mass1 + mass2
    eta = mass1 * mass2 / (M * M)
    m_sec = M * lal.MTSUN_SI
    piM = lal.PI * m_sec

    vISCO = 1. / sqrt(6.)
    fISCO = vISCO * vISCO * vISCO / piM
    f_max = ceilpow2(fISCO)
    n = int(f_max / delta_f + 1)
    kmax = int(fISCO / delta_f)
    kmin = int(numpy.ceil(f_lower / delta_f))
    kmax = kmax if (kmax<n) else n

    #####Calculate the Orientation#####
    v0 = pow(piM *  kmin * delta_f,1./3)
    chi = sqrt(spin1x**2+spin1y**2+spin1z**2)
    kappa = (lnhatx*spin1x+lnhaty*spin1y+lnhatz*spin1z)/chi if (chi > 0.)  else 1.
    Jx0 = mass1*mass2*lnhatx/v0 + mass1*mass1*spin1x
    Jy0 = mass1*mass2*lnhaty/v0 + mass1*mass1*spin1y
    Jz0 = mass1*mass2*lnhatz/v0 + mass1*mass1*spin1z
    thetaJ = acos(Jz0 / sqrt(Jx0**2+Jy0**2+Jz0**2))
    psiJ = atan2(Jy0, -Jx0) # FIXME: check that Jy0 and Jx0 are not both 0
    # Rotate Lnhat back to frame where J is along z, to figure out initial alpha
    rotLx = lnhatx*cos(thetaJ)*cos(psiJ) - lnhaty*cos(thetaJ)*sin(psiJ) + lnhatz*sin(thetaJ)
    rotLy = lnhatx*sin(psiJ) + lnhaty*cos(psiJ)
    alpha0 = atan2(rotLy, rotLx) # FIXME: check that rotLy and rotLx are not both 0
    psiJ_P =psiJ + psi
    psiJ_C =psiJ + psi + lal.PI/4.

    #####Calculate the Coefficients#####
    #quadparam = 1.
    gamma0 = mass1*chi/mass2
    #Calculate the spin corrections
    # FIXME should use pycbc's function, but sigma has different expression
    # in Andy's code, double check
    # pn_beta, pn_sigma, pn_gamma = pycbc.pnutils.mass1_mass2_spin1z_spin2z_to_beta_sigma_gamma(
    #                               mass1, mass2, chi*kappa, 0) # FIXME: spin2 is taken to be 0
    pn_beta = (113.*mass1/(12.*M) - 19.*eta/6.)*chi*kappa
    pn_sigma = (  (5.*(3.*kappa*kappa-1.)/2.) + (7. - kappa*kappa)/96.  ) * (mass1*mass1*chi*chi/M/M)
    pn_gamma = (5.*(146597. + 7056.*eta)*mass1/(2268.*M) - 10.*eta*(1276. + 153.*eta)/81.)*chi*kappa
    prec_fac0 = 5.*(4. + 3.*mass2/mass1)/64.
    dtdv2 = 743./336. + 11.*eta/4.
    dtdv3 = -4.*lal.PI + pn_beta
    dtdv4 = 3058673./1016064. + 5429.*eta/1008. + 617.*eta*eta/144. - pn_sigma
    dtdv5 = (-7729./672.+13.*eta/8.)*lal.PI + 9.*pn_gamma/40.

    #####Calculate the Initial Euler Angles alpha_ref, beta_ref=0 and zeta_ref#####
    gam = gamma0*v0
    sqrtfac = sqrt(1. + 2.*kappa*gam + gam*gam)
    logv0 = log(v0)
    logfac1 = log(1. + kappa*gam + sqrtfac)
    logfac2 = log(kappa + gam + sqrtfac)
    v02 = v0 * v0
    v03 = v0 * v02
    kappa2 = kappa * kappa
    kappa3 = kappa2 * kappa
    gamma02 = gamma0 * gamma0
    gamma03 = gamma02 *gamma0

    alpha_ref = prec_fac0*(  logfac2 *( dtdv2*gamma0 + dtdv3*kappa - dtdv5*kappa/(2.*gamma02) + dtdv4/(2.*gamma0) - dtdv4*kappa2/(2.*gamma0) + (dtdv5*kappa3)/(2.*gamma02) )  +  logfac1*( - dtdv2*gamma0*kappa - dtdv3 + kappa*gamma03/2. - gamma03*kappa3/2. ) + logv0 *( dtdv2*gamma0*kappa + dtdv3 - kappa*gamma03/2. + gamma03*kappa3/2. ) + sqrtfac *( dtdv3 + dtdv4*v0/2. + dtdv5/gamma02/3. + dtdv4*kappa/(2.*gamma0) + dtdv5*kappa*v0/(6.*gamma0) - dtdv5*kappa2/(2.*gamma02) - 1/(3.*v03) - gamma0*kappa/(6.*v02) - dtdv2/v0 - gamma02/(3.*v0) + gamma02*kappa2/(2.*v0) + dtdv5*v02/3. ))  - alpha0

    zeta_ref = prec_fac0*( dtdv3*gamma0*kappa*v0 + dtdv4*v0 + logfac2 *(-dtdv2*gamma0 - dtdv3*kappa + dtdv5*kappa/(2.*gamma02) - dtdv4/(2.*gamma0) + dtdv4*kappa2/(2.*gamma0) - dtdv5*kappa3/(2.*gamma02) ) + logv0 *( kappa*gamma03/2. - gamma03*kappa3/2. ) + logfac1 *( dtdv2*gamma0*kappa + dtdv3 - kappa*gamma03/2. + gamma03*kappa3/2. ) - 1/(3.*v03) - gamma0*kappa/(2.*v02) - dtdv2/v0 + dtdv4*gamma0*kappa*v02/2. + dtdv5*v02/2. + sqrtfac *( -dtdv3 - dtdv4*v0/2. - dtdv5/(3.*gamma02) - dtdv4*kappa/(2.*gamma0) - dtdv5*kappa*v0/(6.*gamma0) + dtdv5*kappa2/(2.*gamma02) + 1/(3.*v03) + gamma0*kappa/(6.*v02) + dtdv2/v0 + gamma02/(3.*v0) - gamma02*kappa2/(2.*v0) - dtdv5*v02/3. ) + dtdv5*gamma0*kappa*v03/3. )

    #####Calculate the Complex sideband factors, mm=2 is first entry#####
    RE_SBfac0= (1.+cos(thetaJ)**2)/2.
    RE_SBfac1= sin(2.*thetaJ)
    RE_SBfac2= 3.*sin(thetaJ)**2
    RE_SBfac3= -sin(2.*thetaJ)
    RE_SBfac4= (1.+cos(thetaJ)**2)/2.
    IM_SBfac0= -cos(thetaJ)
    IM_SBfac1= -2.*sin(thetaJ)
    IM_SBfac2= 0.
    IM_SBfac3= -2.*sin(thetaJ)
    IM_SBfac4= cos(thetaJ)

    #####Calculate the PN terms # FIXME replace with functions in lalsimulation #####
    theta = -11831./9240.
    lambdaa = -1987./3080.0
    pfaN = 3.0/(128.0 * eta)
    pfa2 = 5.0*(743.0/84 + 11.0 * eta)/9.0
    pfa3 = -16.0*lal.PI + 4.0*pn_beta
    pfa4 = 5.0*(3058.673/7.056 + 5429.0/7.0 * eta + 617.0 * eta*eta)/72.0 - \
            10.0*pn_sigma
    pfa5 = 5.0/9.0 * (7729.0/84.0 - 13.0 * eta) * lal.PI - pn_gamma
    pfl5 = 5.0/3.0 * (7729.0/84.0 - 13.0 * eta) * lal.PI - pn_gamma * 3
    pfa6 = (11583.231236531/4.694215680 - 640.0/3.0 * lal.PI * lal.PI- \
            6848.0/21.0*lal.GAMMA) + \
            eta * (-15335.597827/3.048192 + 2255./12. * lal.PI * \
            lal.PI - 1760./3.*theta +12320./9.*lambdaa) + \
            eta*eta * 76055.0/1728.0 - \
            eta*eta*eta*  127825.0/1296.0
    pfl6 = -6848.0/21.0
    pfa7 = lal.PI * 5.0/756.0 * ( 15419335.0/336.0 + 75703.0/2.0 * eta - \
            14809.0 * eta*eta)

    FTaN = 32.0 * eta*eta / 5.0
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

    dETaN = 2 * -eta/2.0
    dETa1 = 2 * -(3.0/4.0 + 1.0/12.0 * eta)
    dETa2 = 3 * -(27.0/8.0 - 19.0/8.0 * eta + 1./24.0 * eta*eta)
    dETa3 = 4 * -(67.5/6.4 - (344.45/5.76 - 20.5/9.6 * lal.PI*lal.PI) *
                             eta + 15.5/9.6 * eta*eta + 3.5/518.4 * eta*eta*eta)

    amp0 = -4. * mass1 * mass2 / (1.0e+06 * distance * lal.PC_SI ) * \
                    lal.MRSUN_SI * lal.MTSUN_SI * sqrt(lal.PI/12.0)

    htildeP = FrequencySeries(zeros(n,dtype=complex128), delta_f=delta_f, copy=False)
    htildeC = FrequencySeries(zeros(n,dtype=complex128), delta_f=delta_f, copy=False)
    spintaylorf2_kernel(htildeP.data[kmin:kmax], htildeC.data[kmin:kmax],
                        kmin, phase_order, amplitude_order, delta_f, piM, pfaN,
                        pfa2, pfa3, pfa4, pfa5, pfl5,
                        pfa6, pfl6, pfa7, FTaN, FTa2,
                        FTa3, FTa4, FTa5, FTa6,
                        FTl6, FTa7, dETaN, dETa1, dETa2,  dETa3,
                        amp0, tC, phi0,
                        kappa, prec_fac0, alpha_ref, zeta_ref,
                        dtdv2, dtdv3, dtdv4, dtdv5,
                        RE_SBfac0, RE_SBfac1, RE_SBfac2, RE_SBfac3, RE_SBfac4,
                        IM_SBfac0, IM_SBfac1, IM_SBfac2, IM_SBfac3, IM_SBfac4,
                        psiJ_P, psiJ_C, gamma0)

    return htildeP, htildeC
