#  Copyright (C) 2012 Prayush Kumar
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
from numpy import sqrt, log, float128

from pycuda.elementwise import ElementwiseKernel

from pycbc.libutils import pkg_config_header_strings
from pycbc.types import FrequencySeries, zeros, Array, complex64

preamble = """
#include <lal/LALConstants.h>
"""

phenomC_text = """
    /* ********* Main paper : Phys Rev D82, 064016 (2010) ********* */

    const double f = (double) (i + kmin ) * delta_f;
    const double fd = (double) m_sec * f;
    const double v = (double) cbrt(piM*f);
    const double v2 = v * v;
    const double v3 = v * v * v;
    const double v4 = v2 * v2;
    const double v5 = v2 * v3;
    const double v6 = v3 * v3;
    const double v7 = v3 * v4;
    const double w = (double) cbrt( m_sec * f );
    const double w3 = (double) w * w * w;
   
    /* ******************************************************* */
    /* *********************** Phasing *********************** */
    /* This is defined in Eq 5.1 - 5.9, 3.13 of the main paper */
    /* ******************************************************* */

    double phSPA = 1. + pfa2 * v2 + pfa3 * v3 + pfa4 * v4 + 
          (1. + log(v3)) * pfa5 * v5 + (pfa6  + pfa6log * log(v3))*v6 + 
          pfa7 * v7;
    phSPA *= (pfaN / v5);
    phSPA -= (LAL_PI/4.0);

    double phPM = (a1/(w3 * w * w)) + (a2/w3) + (a3/w) + a4 + (a5 * w * w) +(a6 * w3);
    phPM /= eta;

    double phRD = b1 + b2*fd;

    double wPlusf1 = 0.5*(1. + tanh( (4*(fd - Mf1)/d1) ));
    double wMinusf1 = 0.5*(1. - tanh( (4*(fd - Mf1)/d1) ));
    double wPlusf2 = 0.5*(1. + tanh( (4*(fd - Mf2)/d2) ));
    double wMinusf2 = 0.5*(1. - tanh( (4*(fd - Mf2)/d2) ));

    double phasing = (phSPA * ((double) wMinusf1)) + (phPM * ((double) wPlusf1 * wMinusf2)) + 
                  (phRD * ((double) wPlusf2));

    /* ******************************************************* */
    /* ********************** Amplitude **************** */
    /* *** This is defined in Eq 5.11 - 5.13, 3.10, 3.6 ****** */
    /* ******************************************************* */

    double xdot = 1. + xdota2 * v2 + xdota3 * v3 + xdota4 * v4 + xdota5 * v5 +
        (xdota6 + xdota6log * log(v2)) * v6 + xdota7 * v7;
    xdot *= (xdotaN * v5 * v5);
    double omgdot = 0.0, ampfac = 0.0;
    double ampSPA = 0.0, ampSPAre = 0.0, ampSPAim = 0.0;
    
    /* If xdot becomes negative, take ampSPA = 0.0 */
    /* This is valid because it becomes negative much after ISCO */

    if( xdot > 0.0 )
    {
      omgdot = 1.5 * v * xdot;
      ampfac = sqrt( LAL_PI / omgdot );
      ampSPAre = ampfac * AN * v2 * (1. + A2 * v2 + A3 * v3 + A4 * v4 + 
              A5 * v5 + (A6 + A6log * log(v2)) * v6);
      ampSPAim = ampfac * AN * v2 * (A5imag * v5 + A6imag * v6);
      ampSPA = sqrt( ampSPAre * ampSPAre + ampSPAim * ampSPAim );
    }

    double ampPM = ampSPA + (g1 * pow(fd, 5./6.));

    const double sig = Mfrd * del2 / Q;
    double sig2 = sig * sig;
    double L = sig2 / ((fd - Mfrd) * (fd - Mfrd) + sig2/4.);
    double ampRD = del1 * L * pow( fd, -7./6.);

    double wPlusf0 = 0.5*(1. + tanh( (4*(fd - Mf0)/d0) ));
    double wMinusf0 = 0.5*(1. - tanh( (4*(fd - Mf0)/d0) ));

    double amplitude = (ampPM * ((double) wMinusf0)) + (ampRD * ((double) wPlusf0));
    amplitude /= distance;

    /* ************** htilde **************** */
    htilde[i]._M_re = amplitude * cos( phasing );
    htilde[i]._M_im = -1.0 * amplitude * sin( phasing );

"""

phenomC_kernel = ElementwiseKernel("""pycuda::complex<double> *htilde, int kmin, double delta_f, 
                                       double eta, double Xi, double distance,
                                       double m_sec, double piM, double Mfrd,
                                       double pfaN, double pfa2, double pfa3, double pfa4, 
                                       double pfa5, double pfa6, double pfa6log, double pfa7,
                                       double a1, double a2, double a3, double a4,
                                       double a5, double a6, double b1, double b2, 
                                       double Mf1, double Mf2, double Mf0, 
                                       double d1, double d2, double d0, 
                                       double xdota2, double xdota3, double xdota4, 
                                       double xdota5, double xdota6, double xdota6log, 
                                       double xdota7, double xdotaN, double AN,
                                       double A2, double A3, double A4, double A5,
                                       double A5imag, double A6, double A6log, double A6imag,
                                       double g1, double del1, double del2, double Q""",
                    phenomC_text, "phenomC_kernel",
                    preamble=preamble, options=pkg_config_header_strings(['lal']))


def FinalSpin( Xi, eta ):
  """Computes the spin of the final BH that gets formed after merger. This is done usingn Eq 5-6 of arXiv:0710.3345"""
  s4 = -0.129
  s5 = -0.384
  t0 = -2.686
  t2 = -3.454
  t3 = 2.353
  etaXi = eta * Xi
  eta2 = eta*eta
  finspin = (Xi + s4*Xi*etaXi + s5*etaXi*eta + t0*etaXi + 2.*(3.**0.5)*eta + t2*eta2 + t3*eta2*eta)
  if finspin > 1.0:
    raise ValueError("Value of final spin > 1.0. Aborting")
  else:
    return finspin

def fRD( a, M):
  """Calculate the ring-down frequency for the final Kerr BH. Using Eq. 5.5 of Main paper"""
  f = (lal.C_SI**3.0 / (2.0*lal.PI*lal.G_SI*M*lal.MSUN_SI)) * (1.5251 - 1.1568*(1.0-a)**0.1292)
  return f

def Qa( a ):
  """Calculate the quality factor of ring-down, using Eq 5.6 of Main paper"""
  return (0.7 + 1.4187*(1.0-a)**-0.4990)

#Functions to calculate the Tanh window, defined in Eq 5.8 of the main paper
def imrphenomc_tmplt(**kwds):
    """ Return an IMRPhenomC waveform using CUDA to generate the phase and amplitude
      Main Paper: arXiv:1005.3306
    """
    # Pull out the input arguments
    f_min = float128(kwds['f_lower'])
    f_max = float128(kwds['f_final'])
    delta_f = float128(kwds['delta_f'])
    distance = float128(kwds['distance'])
    mass1 = float128(kwds['mass1'])
    mass2 = float128(kwds['mass2'])
    spin1z = float128(kwds['spin1z'])
    spin2z = float128(kwds['spin2z'])

    # phi0, tC are taken to be 0 in the paper, sec V-A, first paragraph.
    psi0 = 0. #float128(kwds['phi0'])
    tC= 0. #-1.0 / delta_f 

    if 'out' in kwds:
        out = kwds['out']
    else:
        out = None

    # Calculate binary parameters
    M = mass1 + mass2
    eta = mass1 * mass2 / (M * M)
    Xi = (mass1 * spin1z / M) + (mass2 * spin2z / M)
    Xisum = 2.*Xi
    Xiprod = Xi*Xi
    Xi2 = Xi*Xi

    m_sec = M * lal.MTSUN_SI;
    piM = lal.PI * m_sec;

    ## The units of distance given as input is taken to pe Mpc. Converting to SI
    distance *= (1.0e6 * lal.PC_SI / (2. * sqrt(5. / (64.*lal.PI)) * M * lal.MRSUN_SI * M * lal.MTSUN_SI))

    # Check if the value of f_max is correctly given, else replace with the fCut
    # used in the PhenomB code in lalsimulation. The various coefficients come
    # from Eq.(4.18) of http://arxiv.org/pdf/0710.2335 and 
    # Table I of http://arxiv.org/pdf/0712.0343
    if not f_max:
      f_max = (1.7086 * eta * eta - 0.26592 * eta + 0.28236) / piM

    # Transform the eta, chi to Lambda parameters, using Eq 5.14, Table II of Main
    # paper.
    z101 = -2.417e-03
    z102 = -1.093e-03
    z111 = -1.917e-02
    z110 = 7.267e-02
    z120 = -2.504e-01

    z201 = 5.962e-01
    z202 = -5.600e-02
    z211 = 1.520e-01
    z210 = -2.970e+00
    z220 = 1.312e+01

    z301 = -3.283e+01
    z302 = 8.859e+00
    z311 = 2.931e+01
    z310 = 7.954e+01
    z320 = -4.349e+02

    z401 = 1.619e+02
    z402 = -4.702e+01
    z411 = -1.751e+02
    z410 = -3.225e+02
    z420 = 1.587e+03

    z501 = -6.320e+02
    z502 = 2.463e+02
    z511 = 1.048e+03
    z510 = 3.355e+02
    z520 = -5.115e+03

    z601 = -4.809e+01
    z602 = -3.643e+02
    z611 = -5.215e+02
    z610 = 1.870e+03
    z620 = 7.354e+02

    z701 = 4.149e+00
    z702 = -4.070e+00
    z711 = -8.752e+01
    z710 = -4.897e+01
    z720 = 6.665e+02

    z801 = -5.472e-02
    z802 = 2.094e-02
    z811 = 3.554e-01
    z810 = 1.151e-01
    z820 = 9.640e-01

    z901 = -1.235e+00
    z902 = 3.423e-01
    z911 = 6.062e+00
    z910 = 5.949e+00
    z920 = -1.069e+01

    eta2 = eta*eta
    Xi2 = Xiprod

    # Calculate alphas, gamma, deltas from Table II and Eq 5.14 of Main paper
    a1 = z101 * Xi + z102 * Xi2 + z111 * eta * Xi + z110 * eta + z120 * eta2
    a2 = z201 * Xi + z202 * Xi2 + z211 * eta * Xi + z210 * eta + z220 * eta2
    a3 = z301 * Xi + z302 * Xi2 + z311 * eta * Xi + z310 * eta + z320 * eta2
    a4 = z401 * Xi + z402 * Xi2 + z411 * eta * Xi + z410 * eta + z420 * eta2
    a5 = z501 * Xi + z502 * Xi2 + z511 * eta * Xi + z510 * eta + z520 * eta2
    a6 = z601 * Xi + z602 * Xi2 + z611 * eta * Xi + z610 * eta + z620 * eta2

    g1 = z701 * Xi + z702 * Xi2 + z711 * eta * Xi + z710 * eta + z720 * eta2

    del1 = z801 * Xi + z802 * Xi2 + z811 * eta * Xi + z810 * eta + z820 * eta2
    del2 = z901 * Xi + z902 * Xi2 + z911 * eta * Xi + z910 * eta + z920 * eta2

    # Get the spin of the final BH
    afin = FinalSpin( Xi, eta )
    Q = Qa( abs(afin) )

    # Get the fRD
    frd = fRD( abs(afin), M)
    Mfrd = frd * m_sec

    # Define the frequencies where SPA->PM->RD
    f1 = 0.1 * frd
    Mf1 = m_sec * f1
    f2 = frd
    Mf2 = m_sec * f2
    d1 = 0.005
    d2 = 0.005
    f0 = 0.98 * frd
    Mf0 = m_sec * f0
    d0 = 0.015

    # Now use this frequency for calculation of betas
    # calculate beta1 and beta2, that appear in Eq 5.7 in the main paper.
    b2 = ((-5./3.)* a1 * pow(Mfrd,(-8./3.)) - a2/(Mfrd*Mfrd) - \
      (a3/3.)*pow(Mfrd,(-4./3.)) + (2./3.)* a5 * pow(Mfrd,(-1./3.)) + a6)/eta

    psiPMrd = (a1 * pow(Mfrd,(-5./3.)) + a2/Mfrd + a3 * pow(Mfrd,(-1./3.)) + \
      a4 + a5 * pow(Mfrd,(2./3.)) + a6 * Mfrd)/eta
    b1 = psiPMrd - (b2 * Mfrd)

    ### Calculate the PN coefficients, Eq A3 - A5 of main paper ###
    pfaN = 3.0/(128.0 * eta)
    pfa2 = (3715./756.) + (55.*eta/9.0)
    pfa3 = -16.0*lal.PI + (113./3.)*Xi - 38.*eta*Xisum/3.
    pfa4 = (152.93365/5.08032) - 50.*Xi2 + eta*(271.45/5.04 + 1.25*Xiprod) + \
        3085.*eta2/72.
    pfa5 = lal.PI*(386.45/7.56 - 65.*eta/9.) - \
        Xi*(735.505/2.268 + 130.*eta/9.) + Xisum*(1285.0*eta/8.1 + 170.*eta2/9.) - \
        10.*Xi2*Xi/3. + 10.*eta*Xi*Xiprod
    pfa6 = 11583.231236531/4.694215680 - 640.0*lal.PI*lal.PI/3. - \
        6848.0*lal.GAMMA/21. - 684.8*log(64.)/6.3 + \
        eta*(2255.*lal.PI*lal.PI/12. - 15737.765635/3.048192) + \
        76.055*eta2/1.728 - (127.825*eta2*eta/1.296) + \
        2920.*lal.PI*Xi/3. - (175. - 1490.*eta)*Xi2/3. - \
        (1120.*lal.PI/3. - 1085.*Xi/3.)*eta*Xisum + \
        (269.45*eta/3.36 - 2365.*eta2/6.)*Xiprod

    pfa6log = -6848./63.

    pfa7 = lal.PI*(770.96675/2.54016 + 378.515*eta/1.512 - 740.45*eta2/7.56) - \
        Xi*(20373.952415/3.048192 + 1509.35*eta/2.24 - 5786.95*eta2/4.32) + \
        Xisum*(4862.041225*eta/1.524096 + 1189.775*eta2/1.008 - 717.05*eta2*eta/2.16 - 830.*eta*Xi2/3. + 35.*eta2*Xiprod/3.) - \
        560.*lal.PI*Xi2 + 20.*lal.PI*eta*Xiprod + \
        Xi2*Xi*(945.55/1.68 - 85.*eta) + Xi*Xiprod*(396.65*eta/1.68 + 255.*eta2)


    xdotaN = 64.*eta/5.
    xdota2 = -7.43/3.36 - 11.*eta/4.
    xdota3 = 4.*lal.PI - 11.3*Xi/1.2 + 19.*eta*Xisum/6.
    xdota4 = 3.4103/1.8144 + 5*Xi2 + eta*(13.661/2.016 - Xiprod/8.) + 5.9*eta2/1.8
    xdota5 = -lal.PI*(41.59/6.72 + 189.*eta/8.) - Xi*(31.571/1.008 - 116.5*eta/2.4) + \
          Xisum*(21.863*eta/1.008 - 79.*eta2/6.) - 3*Xi*Xi2/4. + \
          9.*eta*Xi*Xiprod/4.
    xdota6 = 164.47322263/1.39708800 - 17.12*lal.GAMMA/1.05 + \
          16.*lal.PI*lal.PI/3 - 8.56*log(16.)/1.05 + \
          eta*(45.1*lal.PI*lal.PI/4.8 - 561.98689/2.17728) + \
          5.41*eta2/8.96 - 5.605*eta*eta2/2.592 - 80.*lal.PI*Xi/3. + \
          eta*Xisum*(20.*lal.PI/3. - 113.5*Xi/3.6) + \
          Xi2*(64.153/1.008 - 45.7*eta/3.6) - \
          Xiprod*(7.87*eta/1.44 - 30.37*eta2/1.44)

    xdota6log = -856./105.

    xdota7 = -lal.PI*(4.415/4.032 - 358.675*eta/6.048 - 91.495*eta2/1.512) - \
          Xi*(252.9407/2.7216 - 845.827*eta/6.048 + 415.51*eta2/8.64) + \
          Xisum*(158.0239*eta/5.4432 - 451.597*eta2/6.048 + 20.45*eta2*eta/4.32 + 107.*eta*Xi2/6. - 5.*eta2*Xiprod/24.) + \
          12.*lal.PI*Xi2 - Xi2*Xi*(150.5/2.4 + eta/8.) + \
          Xi*Xiprod*(10.1*eta/2.4 + 3.*eta2/8.)


    AN = 8.*eta*sqrt(lal.PI/5.)
    A2 = (-107. + 55.*eta)/42.
    A3 = 2.*lal.PI - 4.*Xi/3. + 2.*eta*Xisum/3.
    A4 = -2.173/1.512 - eta*(10.69/2.16 - 2.*Xiprod) + 2.047*eta2/1.512
    A5 = -10.7*lal.PI/2.1 + eta*(3.4*lal.PI/2.1)

    A5imag = -24.*eta

    A6 = 270.27409/6.46800 - 8.56*lal.GAMMA/1.05 + \
      2.*lal.PI*lal.PI/3. + \
      eta*(4.1*lal.PI*lal.PI/9.6 - 27.8185/3.3264) - \
      20.261*eta2/2.772 + 11.4635*eta*eta2/9.9792 - \
      4.28*log(16.)/1.05

    A6log = -428./105.

    A6imag = 4.28*lal.PI/1.05

    ### Define other parameters needed by waveform generation ###
    kmin = int(f_min / delta_f)
    kmax = int(f_max / delta_f)
    n = kmax + 1;

    if not out:
        htilde = FrequencySeries(zeros(n,dtype=numpy.complex128), delta_f=delta_f, copy=False)
    else:
        if type(out) is not Array:
            raise TypeError("Output must be an instance of Array")
        if len(out) < kmax:
            raise TypeError("Output array is too small")
        if out.dtype != complex64:
            raise TypeError("Output array is the wrong dtype")
        htilde = FrequencySeries(out, delta_f=delta_f, copy=False)

    phenomC_kernel(htilde.data[kmin:kmax], kmin, delta_f, eta, Xi, distance,
                                       m_sec,  piM,  Mfrd,
                                       pfaN,  pfa2,  pfa3,  pfa4, pfa5,  pfa6,  pfa6log,  pfa7,
                                       a1,  a2,  a3,  a4, a5,  a6,  b1,  b2, 
                                       Mf1,  Mf2,  Mf0, d1,  d2,  d0, 
                                       xdota2,  xdota3,  xdota4, xdota5,  xdota6,  xdota6log, 
                                       xdota7,  xdotaN,  AN, A2,  A3,  A4,  A5,
                                       A5imag,  A6,  A6log,  A6imag,
                                       g1,  del1,  del2,  Q )
    hp = htilde
    hc = htilde * 1j
    return hp, hc

