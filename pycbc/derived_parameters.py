# Copyright (C) 2017  Collin Capano
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
from other parameters.
"""

import copy
import numpy
import lal
from pycbc import coordinates
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

# The function library provides the functions that can be used.
function_library = {}

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

function_library[primary_mass.func_name] = primary_mass


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

function_library[secondary_mass.func_name] = secondary_mass


def mtotal_from_mass1_mass2(mass1, mass2):
    """Returns the total mass from mass1 and mass2."""
    return mass1 + mass2

function_library[mtotal_from_mass1_mass2.func_name] = mtotal_from_mass1_mass2


def q_from_mass1_mass2(mass1, mass2):
    """Returns the mass ratio m1/m2, where m1 >= m2."""
    return primary_mass(mass1, mass2) / secondary_mass(mass1, mass2)

function_library[q_from_mass1_mass2.func_name] = q_from_mass1_mass2


def invq_from_mass1_mass2(mass1, mass2):
    """Returns the inverse mass ratio m2/m1, where m1 >= m2."""
    return secondary_mass(mass1, mass2) / primary_mass(mass1, mass2)

function_library[invq_from_mass1_mass2.func_name] = invq_from_mass1_mass2


def eta_from_mass1_mass2(mass1, mass2):
    """Returns the symmetric mass ratio from mass1 and mass2."""
    return mass1*mass2 / (mass1+mass2)**2.

function_library[eta_from_mass1_mass2.func_name] = eta_from_mass1_mass2


def mchirp_from_mass1_mass2(mass1, mass2):
    """Returns the chirp mass from mass1 and mass2."""
    return eta_from_mass1_mass2(mass1, mass2)**(3./5) * (mass1+mass2)

function_library[mchirp_from_mass1_mass2.func_name] = mchirp_from_mass1_mass2


def mass1_from_mtotal_q(mtotal, q):
    """Returns a component mass from the given total mass and mass ratio.
    
    If the mass ratio q is >= 1, the returned mass will be the primary
    (heavier) mass. If q < 1, the returned mass will be the secondary
    (lighter) mass.
    """
    return q*mtotal / (1.+q)

function_library[mass1_from_mtotal_q.func_name] = mass1_from_mtotal_q


def mass2_from_mtotal_q(mtotal, q):
    """Returns a component mass from the given total mass and mass ratio.

    If the mass ratio q is >= 1, the returned mass will be the secondary
    (lighter) mass. If q < 1, the returned mass will be the primary (heavier)
    mass.
    """
    return mtotal / (1.+q)

function_library[mass2_from_mtotal_q.func_name] = mass2_from_mtotal_q


def mass1_from_mtotal_eta(mtotal, eta):
    """Returns the primary mass from the total mass and symmetric mass
    ratio.
    """
    return 0.5 * mtotal * (1.0 + (1.0 - 4.0 * eta)**0.5)

function_library[mass1_from_mtotal_eta.func_name] = mass1_from_mtotal_eta


def mass2_from_mtotal_eta(mtotal, eta):
    """Returns the secondary mass from the total mass and symmetric mass
    ratio.
    """
    return 0.5 * mtotal * (1.0 - (1.0 - 4.0 * eta)**0.5)

function_library[mass2_from_mtotal_eta.func_name] = mass2_from_mtotal_eta


def mtotal_from_mchirp_eta(mchirp, eta):
    """Returns the total mass from the chirp mass and symmetric mass ratio.
    """
    return mchirp / (eta**(3./5.))

function_library[mtotal_from_mchirp_eta.func_name] = mtotal_from_mchirp_eta


def mass1_from_mchirp_eta(mchirp, eta):
    """Returns the primary mass from the chirp mass and symmetric mass ratio.
    """
    mtotal = mtotal_from_mchirp_eta(mchirp, eta)
    return mass1_from_mtotal_eta(mtotal, eta)

function_library[mass1_from_mchirp_eta.func_name] = mass1_from_mchirp_eta


def mass2_from_mchirp_eta(mchirp, eta):
    """Returns the primary mass from the chirp mass and symmetric mass ratio.
    """
    mtotal = mtotal_from_mchirp_eta(mchirp, eta)
    return mass2_from_mtotal_eta(mtotal, eta)

function_library[mass2_from_mchirp_eta.func_name] = mass2_from_mchirp_eta


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
function_library['mass2_from_mchirp_mass1'] = mass2_from_mchirp_mass1

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
function_library['mass_from_knownmass_eta'] = mass_from_knownmass_eta


def mass2_from_mass1_eta(mass1, eta, force_real=True):
    """Returns the secondary mass from the primary mass and symmetric mass
    ratio.
    """
    return mass_from_knownmass_eta(mass1, eta, known_is_secondary=False,
                                   force_real=force_real)
function_library[mass2_from_mass1_eta.func_name] = mass2_from_mass1_eta


def mass1_from_mass2_eta(mass2, eta, force_real=True):
    """Returns the primary mass from the secondary mass and symmetric mass
    ratio.
    """
    return mass_from_knownmass_eta(mass2, eta, known_is_secondary=True,
                                   force_real=force_real)

function_library[mass1_from_mass2_eta.func_name] = mass1_from_mass2_eta


def eta_from_q(q):
    r"""Returns the symmetric mass ratio from the given mass ratio.

    This is given by:

    .. math::
        \eta = \frac{q}{(1+q)^2}.

    Note that the mass ratio may be either < 1 or > 1.
    """
    return q / (1.+q)**2

function_library[eta_from_q.func_name] = eta_from_q


def mass1_from_mchirp_q(mchirp, q):
    """Returns the primary mass from the given chirp mass and mass ratio."""
    return mass1_from_mchirp_eta(mchirp, eta_from_q(q))

function_library[mass1_from_mchirp_q.func_name] = mass1_from_mchirp_q


def mass1_from_mchirp_q(mchirp, q):
    """Returns the secondary mass from the given chirp mass and mass ratio."""
    return mass2_from_mchirp_eta(mchirp, eta_from_q(q))

function_library[mass1_from_mchirp_q.func_name] = mass1_from_mchirp_q


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

function_library[tau0_from_mtotal_eta.func_name] = tau0_from_mtotal_eta


def tau3_from_mtotal_eta(mtotal, eta, f_lower):
    r"""Returns :math:`\tau_0` from the total mass, symmetric mass ratio, and
    the given frequency.
    """
    # convert to seconds
    mtotal = mtotal * lal.MTSUN_SI
    # formulae from arxiv.org:0706.4437
    return _a3(f_lower) / (mtotal**(2./3.) * eta)

function_library[tau3_from_mtotal_eta.func_name] = tau3_from_mtotal_eta


def tau0_from_mass1_mass2(mass1, mass2, f_lower):
    r"""Returns :math:`\tau_0` from the component masses and given frequency.
    """
    mtotal = mass1 + mass2
    eta = eta_from_mass1_mass2(mass1, mass2)
    return tau0_from_mtotal_eta(mtotal, eta)

function_library[tau0_from_mass1_mass2.func_name] = tau0_from_mass1_mass2


def tau3_from_mass1_mass2(mass1, mass2, f_lower):
    r"""Returns :math:`\tau_3` from the component masses and given frequency.
    """
    mtotal = mass1 + mass2
    eta = eta_from_mass1_mass2(mass1, mass2)
    return tau3_from_mtotal_eta(mtotal, eta)

function_library[tau3_from_mass1_mass2.func_name] = tau3_from_mass1_mass2


def mtotal_from_tau0_tau3(tau0, tau3, f_lower):
    r"""Returns total mass from :math:`\tau_0, \tau_3`."""
    mtotal = (tau3 / _a3(f_lower)) / (tau0 / _a0(f_lower))
    # convert back to solar mass units
    return mtotal/lal.MTSUN_SI

function_library[mtotal_from_tau0_tau3.func_name] = mtotal_from_tau0_tau3


def eta_from_tau0_tau3(tau0, tau3, f_lower):
    r"""Returns symmetric mass ratio from :math:`\tau_0, \tau_3`."""
    mtotal = mtotal_from_tau0_tau3(tau0, tau3, f_lower)
    eta = mtotal**(-2./3.) * (_a3(f_lower) / tau3)
    return eta
    
function_library[eta_from_tau0_tau3.func_name] = eta_from_tau0_tau3


def mass1_from_tau0_tau3(tau0, tau3, f_lower):
    r"""Returns the primary mass from the given :math:`\tau_0, \tau_3`."""
    mtotal = mtotal_from_tau0_tau3(tau0, tau3, f_lower)
    eta = eta_from_tau0_tau3(tau0, tau3, f_lower)
    return mass1_from_mtotal_eta(mtotal, eta)

function_library[mass1_from_tau0_tau3.func_name] = mass1_from_tau0_tau3


def mass2_from_tau0_tau3(tau0, tau3, f_lower):
    r"""Returns the secondary mass from the given :math:`\tau_0, \tau_3`."""
    mtotal = mtotal_from_tau0_tau3(tau0, tau3, f_lower)
    eta = eta_from_tau0_tau3(tau0, tau3, f_lower)
    return mass2_from_mtotal_eta(mtotal, eta)

function_library[mass2_from_tau0_tau3.func_name] = mass2_from_tau0_tau3


#
# =============================================================================
#
#                           CBC spin functions
#
# =============================================================================
#
def chi_eff(mass1, mass2, spin1z, spin2z):
    """Returns the effective spin from mass1, mass2, spin1z, and spin2z."""
    return (spin1z * mass1 + spin2z * mass2) / (mass1+mass2)

function_library[chi_eff.func_name] = chi_eff


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

function_library[primary_spin.func_name] = primary_spin


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

function_library[secondary_spin.func_name] = secondary_spin


def primary_spinx(mass1, mass2, spin1x, spin2x):
    """Returns the x-component of the dimensionless spin of the primary mass.
    """
    return primary_spin(mass1, mass2, spin1x, spin2x)


def primary_spiny(mass1, mass2, spin1y, spin2y):
    """Returns the y-component of the dimensionless spin of the primary mass.
    """
    return primary_spin(mass1, mass2, spin1y, spin2y)


def primary_spinz(mass1, mass2, spin1z, spin2z):
    """Returns the z-component of the dimensionless spin of the primary mass.
    """
    return primary_spin(mass1, mass2, spin1z, spin2z)


def secondary_spinx(mass1, mass2, spin1x, spin2x):
    """Returns the x-component of the dimensionless spin of the secondary mass.
    """
    return secondary_spin(mass1, mass2, spin1x, spin2x)


def secondary_spiny(mass1, mass2, spin1y, spin2y):
    """Returns the y-component of the dimensionless spin of the secondary mass.
    """
    return secondary_spin(mass1, mass2, spin1y, spin2y)


def secondary_spinz(mass1, mass2, spin1z, spin2z):
    """Returns the z-component of the dimensionless spin of the secondary mass.
    """
    return secondary_spin(mass1, mass2, spin1z, spin2z)


def spin_a(spinx, spiny, spinz):
    """Returns the dimensionless spin magnitude from the given spin components.
    """
    return coordinates.cartesian_to_spherical_rho(spinx, spiny, spinz)

function_library[spin_a.func_name] = spin_a


def spin_azimuthal(spinx, spinz):
    """Returns the azimuthal spin angle."""
    # do not need to normalize by mass because it cancels
    return coordinates.cartesian_to_spherical_azimuthal(spinx, spiny)

function_library[spin_azimuthal.func_name] = spin_azimuthal


def spin_polar(spinx, spiny, spinz):
    """Returns the polar spin angle."""
    # do not need to normalize by mass because it cancels
    return coordinates.cartesian_to_spherical_polar(spin1x, spin1y, spin1z)

function_library[spin_polar.func_name] = spin_polar


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

function_library[chirp_distance.func_name] = chirp_distance


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
        return tc + detector.time_delay_from_center(ra, dec, tc)
    else:
        other = Detector(ref_frame)
        return tc + detector.time_delay_from_detector(other, ra, dec, tc) 

det_tc = numpy.vectorize(_det_tc)
function_library['det_tc'] = det_tc

__all__ = function_library.keys()
