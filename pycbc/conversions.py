# Copyright (C) 2017  Collin Capano, Christopher M. Biwer
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
"""
This modules provides a library of functions that calculate waveform parameters
from other parameters. All exposed functions in this module's namespace return
one parameter given a set of inputs.
"""

import copy
import numpy
import lal
from pycbc.detector import Detector

#
# =============================================================================
#
#                           Helper functions
#
# =============================================================================
#
def _ensurearray(arg):
    """Ensures that the given argument is a numpy array. If it is not, the
    argument is converted to an array.
    """
    if not isinstance(arg, numpy.ndarray):
        arr = numpy.array(arg)
        arg = arr
    return arg


def _formatreturn(arg):
    """If the given argument is a numpy array with shape (1,), just returns
    that value."""
    if arg.size == 1:
        arg = arg.item()
    return arg

#
# =============================================================================
#
#                           CBC mass functions
#
# =============================================================================
#
def primary_mass(mass1, mass2):
    """Returns the larger of mass1 and mass2 (p = primary)."""
    mass1 = _ensurearray(mass1)
    mass2 = _ensurearray(mass2)
    if mass1.shape != mass2.shape:
        raise ValueError("mass1 and mass2 must have same shape")
    mp = copy.copy(mass1)
    mask = mass1 < mass2
    mp[mask] = mass2[mask]
    return _formatreturn(mp)


def secondary_mass(mass1, mass2):
    """Returns the smaller of mass1 and mass2 (s = secondary)."""
    mass1 = _ensurearray(mass1)
    mass2 = _ensurearray(mass2)
    if mass1.shape != mass2.shape:
        raise ValueError("mass1 and mass2 must have same shape")
    ms = copy.copy(mass2)
    mask = mass1 < mass2
    ms[mask] = mass1[mask]
    return _formatreturn(ms)


def mtotal_from_mass1_mass2(mass1, mass2):
    """Returns the total mass from mass1 and mass2."""
    return mass1 + mass2


def q_from_mass1_mass2(mass1, mass2):
    """Returns the mass ratio m1/m2, where m1 >= m2."""
    return primary_mass(mass1, mass2) / secondary_mass(mass1, mass2)


def invq_from_mass1_mass2(mass1, mass2):
    """Returns the inverse mass ratio m2/m1, where m1 >= m2."""
    return secondary_mass(mass1, mass2) / primary_mass(mass1, mass2)


def eta_from_mass1_mass2(mass1, mass2):
    """Returns the symmetric mass ratio from mass1 and mass2."""
    return mass1*mass2 / (mass1+mass2)**2.


def mchirp_from_mass1_mass2(mass1, mass2):
    """Returns the chirp mass from mass1 and mass2."""
    return eta_from_mass1_mass2(mass1, mass2)**(3./5) * (mass1+mass2)


def mass1_from_mtotal_q(mtotal, q):
    """Returns a component mass from the given total mass and mass ratio.

    If the mass ratio q is >= 1, the returned mass will be the primary
    (heavier) mass. If q < 1, the returned mass will be the secondary
    (lighter) mass.
    """
    return q*mtotal / (1.+q)


def mass2_from_mtotal_q(mtotal, q):
    """Returns a component mass from the given total mass and mass ratio.

    If the mass ratio q is >= 1, the returned mass will be the secondary
    (lighter) mass. If q < 1, the returned mass will be the primary (heavier)
    mass.
    """
    return mtotal / (1.+q)


def mass1_from_mtotal_eta(mtotal, eta):
    """Returns the primary mass from the total mass and symmetric mass
    ratio.
    """
    return 0.5 * mtotal * (1.0 + (1.0 - 4.0 * eta)**0.5)


def mass2_from_mtotal_eta(mtotal, eta):
    """Returns the secondary mass from the total mass and symmetric mass
    ratio.
    """
    return 0.5 * mtotal * (1.0 - (1.0 - 4.0 * eta)**0.5)


def mtotal_from_mchirp_eta(mchirp, eta):
    """Returns the total mass from the chirp mass and symmetric mass ratio.
    """
    return mchirp / (eta**(3./5.))


def mass1_from_mchirp_eta(mchirp, eta):
    """Returns the primary mass from the chirp mass and symmetric mass ratio.
    """
    mtotal = mtotal_from_mchirp_eta(mchirp, eta)
    return mass1_from_mtotal_eta(mtotal, eta)


def mass2_from_mchirp_eta(mchirp, eta):
    """Returns the primary mass from the chirp mass and symmetric mass ratio.
    """
    mtotal = mtotal_from_mchirp_eta(mchirp, eta)
    return mass2_from_mtotal_eta(mtotal, eta)


def _mass2_from_mchirp_mass1(mchirp, mass1):
    r"""Returns the secondary mass from the chirp mass and primary mass.

    As this is a cubic equation this requires finding the roots and returning
    the one that is real. Basically it can be shown that:

    .. math::
        m_2^3 - a(m_2 + m_1) = 0,
 
    where

    .. math::
        a = \frac{\mathcal{M}^5}{m_1^3}.

    This has 3 solutions but only one will be real.
    """
    a = mchirp**5 / mass1**3
    roots = numpy.roots([1,0,-a,-a*mass1])
    # Find the real one
    real_root = roots[(abs(roots - roots.real)).argmin()]
    return real_root.real

mass2_from_mchirp_mass1 = numpy.vectorize(_mass2_from_mchirp_mass1)


def _mass_from_knownmass_eta(known_mass, eta, known_is_secondary=False,
                            force_real=True):
    r"""Returns the other component mass given one of the component masses
    and the symmetric mass ratio.

    This requires finding the roots of the quadratic equation:

    .. math::
        \eta m_2^2 + (2\eta - 1)m_1 m_2 + \eta m_1^2 = 0.

    This has two solutions which correspond to :math:`m_1` being the heavier
    mass or it being the lighter mass. By default, `known_mass` is assumed to
    be the heavier (primary) mass, and the smaller solution is returned. Use
    the `other_is_secondary` to invert.

    Parameters
    ----------
    known_mass : float
        The known component mass.
    eta : float
        The symmetric mass ratio.
    known_is_secondary : {False, bool}
        Whether the known component mass is the primary or the secondary. If
        True, `known_mass` is assumed to be the secondary (lighter) mass and
        the larger solution is returned. Otherwise, the smaller solution is
        returned. Default is False.
    force_real : {True, bool}
        Force the returned mass to be real.

    Returns
    -------
    float
        The other component mass.
    """
    roots = numpy.roots([eta, (2*eta - 1)*known_mass, eta*known_mass**2.])
    if force_real:
        roots = numpy.real(roots)
    if known_is_secondary:
        return roots[roots.argmax()]
    else:
        return roots[roots.argmin()]

mass_from_knownmass_eta = numpy.vectorize(_mass_from_knownmass_eta)


def mass2_from_mass1_eta(mass1, eta, force_real=True):
    """Returns the secondary mass from the primary mass and symmetric mass
    ratio.
    """
    return mass_from_knownmass_eta(mass1, eta, known_is_secondary=False,
                                   force_real=force_real)


def mass1_from_mass2_eta(mass2, eta, force_real=True):
    """Returns the primary mass from the secondary mass and symmetric mass
    ratio.
    """
    return mass_from_knownmass_eta(mass2, eta, known_is_secondary=True,
                                   force_real=force_real)


def eta_from_q(q):
    r"""Returns the symmetric mass ratio from the given mass ratio.

    This is given by:

    .. math::
        \eta = \frac{q}{(1+q)^2}.

    Note that the mass ratio may be either < 1 or > 1.
    """
    return q / (1.+q)**2


def mass1_from_mchirp_q(mchirp, q):
    """Returns the primary mass from the given chirp mass and mass ratio."""
    return mass1_from_mchirp_eta(mchirp, eta_from_q(q))


def mass2_from_mchirp_q(mchirp, q):
    """Returns the secondary mass from the given chirp mass and mass ratio."""
    return mass2_from_mchirp_eta(mchirp, eta_from_q(q))


def _a0(f_lower):
    """Used in calculating chirp times: see Cokelaer, arxiv.org:0706.4437
       appendix 1, also lalinspiral/python/sbank/tau0tau3.py.
    """
    return 5. / (256. * (numpy.pi * f_lower)**(8./3.))

def _a3(f_lower):
    """Another parameter used for chirp times"""
    return numpy.pi / (8. * (numpy.pi * f_lower)**(5./3.))
  

def tau0_from_mtotal_eta(mtotal, eta, f_lower):
    r"""Returns :math:`\tau_0` from the total mass, symmetric mass ratio, and
    the given frequency.
    """
    # convert to seconds
    mtotal = mtotal * lal.MTSUN_SI
    # formulae from arxiv.org:0706.4437
    return _a0(f_lower) / (mtotal**(5./3.) * eta)


def tau3_from_mtotal_eta(mtotal, eta, f_lower):
    r"""Returns :math:`\tau_0` from the total mass, symmetric mass ratio, and
    the given frequency.
    """
    # convert to seconds
    mtotal = mtotal * lal.MTSUN_SI
    # formulae from arxiv.org:0706.4437
    return _a3(f_lower) / (mtotal**(2./3.) * eta)


def tau0_from_mass1_mass2(mass1, mass2, f_lower):
    r"""Returns :math:`\tau_0` from the component masses and given frequency.
    """
    mtotal = mass1 + mass2
    eta = eta_from_mass1_mass2(mass1, mass2)
    return tau0_from_mtotal_eta(mtotal, eta, f_lower)


def tau3_from_mass1_mass2(mass1, mass2, f_lower):
    r"""Returns :math:`\tau_3` from the component masses and given frequency.
    """
    mtotal = mass1 + mass2
    eta = eta_from_mass1_mass2(mass1, mass2)
    return tau3_from_mtotal_eta(mtotal, eta, f_lower)


def mtotal_from_tau0_tau3(tau0, tau3, f_lower,
                          in_seconds=False):
    r"""Returns total mass from :math:`\tau_0, \tau_3`."""
    mtotal = (tau3 / _a3(f_lower)) / (tau0 / _a0(f_lower))
    if not in_seconds:
        # convert back to solar mass units
        mtotal /= lal.MTSUN_SI
    return mtotal


def eta_from_tau0_tau3(tau0, tau3, f_lower):
    r"""Returns symmetric mass ratio from :math:`\tau_0, \tau_3`."""
    mtotal = mtotal_from_tau0_tau3(tau0, tau3, f_lower,
                                   in_seconds=True)
    eta = mtotal**(-2./3.) * (_a3(f_lower) / tau3)
    return eta
    

def mass1_from_tau0_tau3(tau0, tau3, f_lower):
    r"""Returns the primary mass from the given :math:`\tau_0, \tau_3`."""
    mtotal = mtotal_from_tau0_tau3(tau0, tau3, f_lower)
    eta = eta_from_tau0_tau3(tau0, tau3, f_lower)
    return mass1_from_mtotal_eta(mtotal, eta)


def mass2_from_tau0_tau3(tau0, tau3, f_lower):
    r"""Returns the secondary mass from the given :math:`\tau_0, \tau_3`."""
    mtotal = mtotal_from_tau0_tau3(tau0, tau3, f_lower)
    eta = eta_from_tau0_tau3(tau0, tau3, f_lower)
    return mass2_from_mtotal_eta(mtotal, eta)


#
# =============================================================================
#
#                           CBC spin functions
#
# =============================================================================
#
def chi_eff(mass1, mass2, spin1z, spin2z):
    """Returns the effective spin from mass1, mass2, spin1z, and spin2z."""
    return (spin1z * mass1 + spin2z * mass2) / (mass1 + mass2)

def chi_a(mass1, mass2, spin1z, spin2z):
    """ Returns the aligned mass-weighted spin difference from mass1, mass2,
    spin1z, and spin2z.
    """
    m_p = primary_mass(mass1, mass2)
    spin_p = primary_spin(mass1, mass2, spin1z, spin2z)
    m_s = secondary_mass(mass1, mass2)
    spin_s = secondary_spin(mass1, mass2, spin1z, spin2z)
    return (spin_s * m_s - spin_p * m_p) / (m_s + m_p)

def chi_p(mass1, mass2, spin1x, spin1y, spin2x, spin2y):
    """Returns the effective precession spin from mass1, mass2, spin1x,
    spin1y, spin2x, and spin2y.
    """
    xi1 = secondary_xi(mass1, mass2, spin1x, spin1y, spin2x, spin2y)
    xi2 = primary_xi(mass1, mass2, spin1x, spin1y, spin2x, spin2y)
    return chi_p_from_xi1_xi2(xi1, xi2)

def phi_a(mass1, mass2, spin1x, spin1y, spin2x, spin2y):
    """ Returns the angle between the in-plane perpendicular spins."""
    phi1 = phi_from_spinx_spiny(primary_spin(mass1, mass2, spin1x, spin2x),
                                primary_spin(mass1, mass2, spin1y, spin2y))
    phi2 = phi_from_spinx_spiny(secondary_spin(mass1, mass2, spin1x, spin2x),
                                secondary_spin(mass1, mass2, spin1y, spin2y))
    return phi1 - phi2

def phi_s(spin1x, spin1y, spin2x, spin2y):
    """ Returns the sum of the in-plane perpendicular spins."""
    phi1 = phi_from_spinx_spiny(spin1x, spin1y)
    phi2 = phi_from_spinx_spiny(spin2x, spin2y)
    return phi1 + phi2

def primary_spin(mass1, mass2, spin1, spin2):
    """Returns the dimensionless spin of the primary mass."""
    mass1 = _ensurearray(mass1)
    mass2 = _ensurearray(mass2)
    spin1 = _ensurearray(spin1)
    spin2 = _ensurearray(spin2)
    if (mass1.shape != mass2.shape) or (mass1.shape != spin1.shape) or (
        mass1.shape != spin2.shape):
        raise ValueError("mass1, mass2, spin1, spin2 must have same shape")
    sp = copy.copy(spin1)
    mask = mass1 < mass2
    sp[mask] = spin2[mask]
    return _formatreturn(sp)

def secondary_spin(mass1, mass2, spin1, spin2):
    """Returns the dimensionless spin of the secondary mass."""
    mass1 = _ensurearray(mass1)
    mass2 = _ensurearray(mass2)
    spin1 = _ensurearray(spin1)
    spin2 = _ensurearray(spin2)
    if (mass1.shape != mass2.shape) or (mass1.shape != spin1.shape) or (
        mass1.shape != spin2.shape):
        raise ValueError("mass1, mass2, spin1, spin2 must have same shape")
    ss = copy.copy(spin2)
    mask = mass1 < mass2
    ss[mask] = spin1[mask]
    return _formatreturn(ss)

def primary_xi(mass1, mass2, spin1x, spin1y, spin2x, spin2y):
    """Returns the effective precession spin argument for the larger mass.
    """
    spinx = primary_spin(mass1, mass2, spin1x, spin2x)
    spiny = primary_spin(mass1, mass2, spin1y, spin2y)
    return chi_perp_from_spinx_spiny(spinx, spiny)

def secondary_xi(mass1, mass2, spin1x, spin1y, spin2x, spin2y):
    """Returns the effective precession spin argument for the smaller mass.
    """
    spinx = secondary_spin(mass1, mass2, spin1x, spin2x)
    spiny = secondary_spin(mass1, mass2, spin1y, spin2y)
    return xi2_from_mass1_mass2_spin2x_spin2y(mass1, mass2, spinx, spiny)

def xi1_from_spin1x_spin1y(spin1x, spin1y):
    """Returns the effective precession spin argument for the larger mass.
    This function assumes it's given spins of the primary mass.
    """
    return chi_perp_from_spinx_spiny(spin1x, spin1y)

def xi2_from_mass1_mass2_spin2x_spin2y(mass1, mass2, spin2x, spin2y):
    """Returns the effective precession spin argument for the smaller mass.
    This function assumes it's given spins of the secondary mass.
    """
    q = q_from_mass1_mass2(mass1, mass2)
    a1 = 2 + 3 * q / 2
    a2 = 2 + 3 / (2 * q)
    return a1 / (q**2 * a2) * chi_perp_from_spinx_spiny(spin2x, spin2y)

def chi_perp_from_spinx_spiny(spinx, spiny):
    """Returns the in-plane spin from the x/y components of the spin.
    """
    return numpy.sqrt(spinx**2 + spiny**2)

def chi_perp_from_mass1_mass2_xi2(mass1, mass2, xi2):
    """Returns the in-plane spin from mass1, mass2, and xi2 for the
    secondary mass.
    """
    q = q_from_mass1_mass2(mass1, mass2)
    a1 = 2 + 3 * q / 2
    a2 = 2 + 3 / (2 * q)
    return q**2 * a2 / a1 * xi2

def chi_p_from_xi1_xi2(xi1, xi2):
    """Returns effective precession spin from xi1 and xi2.
    """
    xi1 = _ensurearray(xi1)
    xi2 = _ensurearray(xi2)
    if xi1.shape != xi2.shape:
        raise ValueError("xi1, xi2 must have same shape")
    chi_p = copy.copy(xi1)
    mask = xi1 < xi2
    chi_p[mask] = xi2[mask]
    return _formatreturn(chi_p)

def phi1_from_phi_a_phi_s(phi_a, phi_s):
    """Returns the angle between the x-component axis and the in-plane
    spin for the primary mass from phi_s and phi_a.
    """
    return (phi_s + phi_a) / 2.0

def phi2_from_phi_a_phi_s(phi_a, phi_s):
    """Returns the angle between the x-component axis and the in-plane
    spin for the secondary mass from phi_s and phi_a.
    """
    return (phi_s - phi_a) / 2.0

def phi_from_spinx_spiny(spinx, spiny):
    """Returns the angle between the x-component axis and the in-plane spin.
    """
    return numpy.arctan(spiny / spinx)

def spin1z_from_mass1_mass2_chi_eff_chi_a(mass1, mass2, chi_eff, chi_a):
    """Returns spin1z.
    """
    return (mass1 + mass2) / (2 * mass1) * (chi_eff - chi_a)

def spin2z_from_mass1_mass2_chi_eff_chi_a(mass1, mass2, chi_eff, chi_a):
    """Returns spin2z.
    """
    return (mass1 + mass2) / (2 * mass2) * (chi_eff + chi_a)

def spin1x_from_xi1_phi_a_phi_s(xi1, phi_a, phi_s):
    """Returns x-component spin for primary mass.
    """
    phi1 = phi1_from_phi_a_phi_s(phi_a, phi_s)
    return xi1 * numpy.cos(phi1)

def spin1y_from_xi1_phi_a_phi_s(xi1, phi_a, phi_s):
    """Returns y-component spin for primary mass.
    """
    phi1 = phi1_from_phi_a_phi_s(phi_s, phi_a)
    return xi1 * numpy.sin(phi1)

def spin2x_from_mass1_mass2_xi2_phi_a_phi_s(mass1, mass2, xi2, phi_a, phi_s):
    """Returns x-component spin for secondary mass.
    """
    chi_perp = chi_perp_from_mass1_mass2_xi2(mass1, mass2, xi2)
    phi2 = phi2_from_phi_a_phi_s(phi_a, phi_s)
    return chi_perp * numpy.cos(phi2)

def spin2y_from_mass1_mass2_xi2_phi_a_phi_s(mass1, mass2, xi2, phi_a, phi_s):
    """Returns y-component spin for secondary mass.
    """
    chi_perp = chi_perp_from_mass1_mass2_xi2(mass1, mass2, xi2)
    phi2 = phi2_from_phi_a_phi_s(phi_a, phi_s)
    return chi_perp * numpy.sin(phi2)

#
# =============================================================================
#
#                         Extrinsic parameter functions
#
# =============================================================================
#
def chirp_distance(dist, mchirp, ref_mass=1.4):
    """Returns the chirp distance given a distance and chirp mass.
    """
    return dist * (2.**(-1./5) * ref_mass / mchirp)**(5./6)


def _det_tc(detector_name, ra, dec, tc, ref_frame='geocentric'):
    """Returns the coalescence time of a signal in the given detector.
    
    Parameters
    ----------
    detector_name : string
        The name of the detector, e.g., 'H1'.
    ra : float
        The right ascension of the signal, in radians.
    dec : float
        The declination of the signal, in radians.
    tc : float
        The GPS time of the coalescence of the signal in the `ref_frame`.
    ref_frame : {'geocentric', string}
        The reference frame that the given coalescence time is defined in.
        May specify 'geocentric', or a detector name; default is 'geocentric'.

    Returns
    -------
    float :
        The GPS time of the coalescence in detector `detector_name`.
    """
    if ref_frame == detector_name:
        return tc
    detector = Detector(detector_name)
    if ref_frame == 'geocentric':
        return tc + detector.time_delay_from_earth_center(ra, dec, tc)
    else:
        other = Detector(ref_frame)
        return tc + detector.time_delay_from_detector(other, ra, dec, tc)

det_tc = numpy.vectorize(_det_tc)

#
# =============================================================================
#
#                         Likelihood statistic parameter functions
#
# =============================================================================
#
def snr_from_loglr(loglr):
    """Returns SNR computed from the given log likelihood ratio(s). This is
    defined as `sqrt(2*loglr)`.If the log likelihood ratio is < 0, returns 0.

    Parameters
    ----------
    loglr : array or float
        The log likelihood ratio(s) to evaluate.

    Returns
    -------
    array or float
        The SNRs computed from the log likelihood ratios.
    """
    singleval = isinstance(loglr, float)
    if singleval:
        loglr = numpy.array([loglr])
    # temporarily quiet sqrt(-1) warnings
    numpysettings = numpy.seterr(invalid='ignore')
    snrs = numpy.sqrt(2*loglr)
    numpy.seterr(**numpysettings)
    snrs[numpy.isnan(snrs)] = 0.
    if singleval:
        snrs = snrs[0]
    return snrs


__all__ = ['primary_mass', 'secondary_mass', 'mtotal_from_mass1_mass2',
           'q_from_mass1_mass2', 'invq_from_mass1_mass2',
           'eta_from_mass1_mass2', 'mchirp_from_mass1_mass2',
           'mass1_from_mtotal_q', 'mass2_from_mtotal_q',
           'mass1_from_mtotal_eta', 'mass2_from_mtotal_eta',
           'mtotal_from_mchirp_eta', 'mass1_from_mchirp_eta',
           'mass2_from_mchirp_eta', 'mass2_from_mchirp_mass1',
           'mass_from_knownmass_eta', 'mass2_from_mass1_eta',
           'mass1_from_mass2_eta', 'eta_from_q', 'mass1_from_mchirp_q',
           'mass2_from_mchirp_q', 'tau0_from_mtotal_eta',
           'tau3_from_mtotal_eta', 'tau0_from_mass1_mass2',
           'tau3_from_mass1_mass2', 'mtotal_from_tau0_tau3',
           'eta_from_tau0_tau3', 'mass1_from_tau0_tau3',
           'mass2_from_tau0_tau3', 'primary_spin', 'secondary_spin',
           'chi_eff', 'chi_a', 'chi_p', 'phi_a', 'phi_s',
           'primary_xi', 'secondary_xi',
           'xi1_from_spin1x_spin1y', 'xi2_from_mass1_mass2_spin2x_spin2y',
           'chi_perp_from_spinx_spiny', 'chi_perp_from_mass1_mass2_xi2',
           'chi_p_from_xi1_xi2', 'phi_from_spinx_spiny',
           'phi1_from_phi_a_phi_s', 'phi2_from_phi_a_phi_s',
           'spin1z_from_mass1_mass2_chi_eff_chi_a',
           'spin2z_from_mass1_mass2_chi_eff_chi_a',
           'spin1x_from_xi1_phi_a_phi_s', 'spin1y_from_xi1_phi_a_phi_s',
           'spin2x_from_mass1_mass2_xi2_phi_a_phi_s',
           'spin2y_from_mass1_mass2_xi2_phi_a_phi_s',
           'chirp_distance', 'det_tc', 'snr_from_loglr',
          ]
