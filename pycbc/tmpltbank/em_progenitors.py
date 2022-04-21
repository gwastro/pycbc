# Copyright (C) 2015 Francesco Pannarale
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

import math
import numpy as np
import pycbc
import os.path
import scipy.optimize
import scipy.interpolate
import logging
import sys

##############################################################################
# Innermost Stable Spherical Orbit (ISSO) solver in the Perez-Giz (PG)       #
# formalism [see Stone, Loeb, Berger, PRD 87, 084053 (2013)].                #
##############################################################################

# Equation that determines the ISCO radius (in BH mass units)
def ISCO_eq(r, chi):
    """
    Polynomial that enables the calculation of the Kerr
    inntermost stable circular orbit (ISCO) radius via its
    roots.

    Parameters
    -----------
    r: float
        the radial coordinate in BH mass units
    chi: float
        the BH dimensionless spin parameter

    Returns
    ----------
    float
        (r*(r-6))**2-chi**2*(2*r*(3*r+14)-9*chi**2)
    """
    return (r*(r-6))**2-chi**2*(2*r*(3*r+14)-9*chi**2)

# Equation that determines the ISSO radius (in BH mass units) at one of the
# poles
def ISSO_eq_at_pole(r, chi):
    """
    Polynomial that enables the calculation of the Kerr polar
    (inclination = +/- pi/2) innermost stable spherical orbit
    (ISSO) radius via its roots.  Physical solutions are
    between 6 and 1+sqrt[3]+sqrt[3+2sqrt[3]].

    Parameters
    ----------
    r: float
        the radial coordinate in BH mass units
    chi: float
        the BH dimensionless spin parameter

    Returns
    -------
    float
        r**3*(r**2*(r-6)+chi**2*(3*r+4))+chi**4*(3*r*(r-2)+chi**2)
    """
    return r**3*(r**2*(r-6)+chi**2*(3*r+4))+chi**4*(3*r*(r-2)+chi**2)

# Equation that determines the ISSO radius (in BH mass units) for a generic
# orbital inclination
def PG_ISSO_eq(r, chi, incl):
    """Polynomial that enables the calculation of a generic innermost
    stable spherical orbit (ISSO) radius via its roots.  Physical
    solutions are between the equatorial ISSO (aka the ISCO) radius
    and the polar ISSO radius. See Stone, Loeb, Berger, PRD 87, 084053 (2013).

    Parameters
    ----------
    r: float
        the radial coordinate in BH mass units
    chi: float
        the BH dimensionless spin parameter
    incl: float
        inclination angle between the BH spin and the orbital angular
        momentum in radians

    Returns
    -------
    float
        ``r**8*Z+chi**2*(1-cos_incl**2)*(chi**2*(1-cos_incl**2)*Y-2*r**4*X)``
        where
        ``X=chi**2*(chi**2*(3*chi**2+4*r*(2*r-3))+r**2*(15*r*(r-4)+28))-6*r**4*(r**2-4)``
        ``Y=chi**4*(chi**4+r**2*(7*r*(3*r-4)+36))+6*r*(r-2)*(chi**6+2*r**3*(chi**2*(3*r+2)+3*r**2*(r-2)))``
        ``Z=ISCO_eq(r,chi)``
    """
    chi2 = chi*chi
    chi4 = chi2*chi2
    r2 = r*r
    r4 = r2*r2
    three_r = 3*r
    r_minus_2 = r - 2
    sin_incl2 = (math.sin(incl))**2

    X = chi2*(chi2*(3*chi2+4*r*(2*r-3))+r2*(15*r*(r-4)+28))-6*r4*(r2-4)
    Y = chi4 * (chi4+r2*(7 * r * (three_r-4)+36))+6 * r * r_minus_2 * \
        (chi4*chi2+2*r2*r*(chi2*(three_r+2)+3*r2*r_minus_2))
    Z = ISCO_eq(r, chi)

    return r4*r4*Z+chi2*sin_incl2*(chi2*sin_incl2*Y-2*r4*X)

# ISSO radius solver
def PG_ISSO_solver(chi,incl):
    """Function that determines the radius of the innermost stable
    spherical orbit (ISSO) for a Kerr BH and a generic inclination
    angle between the BH spin and the orbital angular momentum.
    This function finds the appropriat root of PG_ISSO_eq.

    Parameters
    ----------
    chi: float
        the BH dimensionless spin parameter
    incl: float
        the inclination angle between the BH spin and the orbital
        angular momentum in radians

    Returns
    -------
    solution: float
        the radius of the orbit in BH mass units
    """
    # Auxiliary variables
    cos_incl=math.cos(incl)
    sgnchi = np.sign(cos_incl)*chi

    # ISCO radius for the given spin magnitude
    if sgnchi > 0.99:
        initial_guess = 2 # Works better than 6 for chi=1
    elif sgnchi < 0:
        initial_guess = 9
    else:
        initial_guess = 5 # Works better than 6 for chi=0.5
    rISCO_limit = scipy.optimize.fsolve(ISCO_eq, initial_guess, args=sgnchi)

    # ISSO radius for an inclination of pi/2
    if chi < 0:
        initial_guess = 9
    else:
        initial_guess = 6
    rISSO_at_pole_limit = scipy.optimize.fsolve(ISSO_eq_at_pole, initial_guess,
                                                args=chi)

    # If the inclination is 0 or pi, just output the ISCO radius
    if incl in [0, math.pi]:
        solution = rISCO_limit
    # If the inclination is pi/2, just output the ISSO radius at the pole(s)
    elif incl == math.pi/2:
        solution = rISSO_at_pole_limit
    # Otherwise, find the ISSO radius for a generic inclination
    else:
        initial_guess = max(rISCO_limit,rISSO_at_pole_limit)
        solution = scipy.optimize.fsolve(PG_ISSO_eq, initial_guess,
                                         args=(chi, incl))
        if solution < 1 or solution > 9:
            initial_guess = min(rISCO_limit,rISSO_at_pole_limit)
            solution = scipy.optimize.fsolve(PG_ISSO_eq, initial_guess,
                                             args=(chi, incl))

    return solution

##############################################################################
# 2H 2-piecewise polytropic EOS, NS non-rotating equilibrium sequence        #
# File format is: grav mass (Msun)   baryonic mass (Msun)    compactness     #
#                                                                            #
# Eventually, the function should call an NS sequence generator within LAL   #
# the EOS prescribed by the user and store it.                               #
##############################################################################
def load_ns_sequence(eos_name):
    """
    Load the data of an NS non-rotating equilibrium sequence
    generated using the equation of state (EOS) chosen by the
    user.  [Only the 2H 2-piecewise polytropic EOS is currently
    supported.  This yields NSs with large radiss (15-16km).]

    Parameters
    -----------
    eos_name: string
        NS equation of state label ('2H' is the only supported
        choice at the moment)

    Returns
    ----------
    ns_sequence: 3D-array
        contains the sequence data in the form NS gravitational
         mass (in solar masses), NS baryonic mass (in solar
         masses), NS compactness (dimensionless)
    max_ns_g_mass: float
        the maximum NS gravitational mass (in solar masses) in
        the sequence (this is the mass of the most massive stable
        NS)
    """
    ns_sequence = []

    if eos_name == '2H':
        ns_sequence_path = os.path.join(pycbc.tmpltbank.NS_SEQUENCE_FILE_DIRECTORY, 'equil_2H.dat')
        ns_sequence = np.loadtxt(ns_sequence_path)
    else:
        err_msg = "Only the 2H EOS is currently supported."
        err_msg += "If you plan to use a different NS EOS,"
        err_msg += "be sure not to filter too many templates!\n"
        logging.error(err_msg)
        sys.exit(1)
        raise Exception('Unsupported EOS!')

    max_ns_g_mass = max(ns_sequence[:,0])

    return (ns_sequence, max_ns_g_mass)

##############################################################################
# Given an NS equilibrium sequence and gravitational mass (in Msun), return  #
# the NS baryonic mass (in Msun).                                            #
##############################################################################
def ns_g_mass_to_ns_b_mass(ns_g_mass, ns_sequence):
    """
    Determines the baryonic mass of an NS given its gravitational
    mass and an NS equilibrium sequence.

    Parameters
    -----------
    ns_g_mass: float
        NS gravitational mass (in solar masses)
    ns_sequence: 3D-array
        contains the sequence data in the form NS gravitational
         mass (in solar masses), NS baryonic mass (in solar
         masses), NS compactness (dimensionless)

    Returns
    ----------
    float
        The NS baryonic mass (in solar massesr**3*(r**2*(r-6)+chi**2*(3*r+4))+
        chi**4*(3*r*(r-2)+chi**2))
    """
    x = ns_sequence[:,0]
    y = ns_sequence[:,1]
    f = scipy.interpolate.interp1d(x, y)

    return f(ns_g_mass)

##############################################################################
# Given an NS equilibrium sequence and gravitational mass (in Msun), return  #
# the NS compactness.                                                        #
##############################################################################
def ns_g_mass_to_ns_compactness(ns_g_mass, ns_sequence):
    """
    Determines the compactness of an NS given its
    gravitational mass and an NS equilibrium sequence.

    Parameters
    -----------
    ns_g_mass: float
        NS gravitational mass (in solar masses)
    ns_sequence: 3D-array
        contains the sequence data in the form NS gravitational
         mass (in solar masses), NS baryonic mass (in solar
         masses), NS compactness (dimensionless)

    Returns
    ----------
    float
        The NS compactness (dimensionless)
    """
    x = ns_sequence[:,0]
    y = ns_sequence[:,2]
    f = scipy.interpolate.interp1d(x, y)

    return f(ns_g_mass)

##############################################################################
# NS-BH merger remnant mass [Foucart, Hinderer, Nissanke PRD 98, 081501(R)   #
# (2018)].                                                                   #
#                                                                            #
# * THIS ASSUMES THE NS SPIN IS 0 (the user is warned about this).           #
# * Physical parameters passed to remnant_mass (see sanity checks below)     #
#   must not be unphysical.                                                  #
# * Tilted BH spin cases are handled by using rISSO(chi,incl) where rISCO    #
#   appears in the formula. This approximation was suggested in Stone, Loeb, #
#   Berger, PRD 87, 084053 (2013) for the previous NS-BH remnant mass fit of #
#   Foucart, PRD 86, 124007 (2012).                                          #
##############################################################################
def remnant_mass(eta, ns_g_mass, ns_sequence, chi, incl):
    """
    Function that determines the remnant disk mass of
    an NS-BH system using the fit to numerical-relativity
    results discussed in Foucart, Hinderer, Nissanke PRD 98, 081501(R) (2018).

    Parameters
    -----------
    eta: float
        the symmetric mass ratio of the binary
    ns_g_mass: float
        NS gravitational mass (in solar masses)
    ns_sequence: 3D-array
        contains the sequence data in the form NS gravitational
         mass (in solar masses), NS baryonic mass (in solar
         masses), NS compactness (dimensionless)
    chi: float
        the BH dimensionless spin parameter
    incl: float
        the inclination angle between the BH spin and the orbital
        angular momentum in radians

    Returns
    ----------
    remnant_mass: float
        The remnant mass in solar masses
    """

    # Sanity checks on eta and chi
    if not (eta>0. and eta<=0.25 and abs(chi)<=1 and ns_g_mass>0):
        err_msg = "The BH spin magnitude must be <= 1,"
        err_msg += "eta must be between 0 and 0.25,"
        err_msg += "and the NS mass must be positive."
        logging.error(err_msg)
        info_msg = "The function remnant_mass was launched with"
        info_msg += "ns_mass={0}, eta={1}, chi={2},"
        info_msg += "inclination={3}\n'.format(ns_g_mass, eta, chi, incl)"
        logging.info(info_msg)
        sys.exit(1)
        raise Exception('Unphysical parameters!')

    # NS compactness and rest mass
    ns_compactness = ns_g_mass_to_ns_compactness(ns_g_mass, ns_sequence)
    ns_b_mass = ns_g_mass_to_ns_b_mass(ns_g_mass, ns_sequence)

    # Fit parameters and tidal correction
    alpha = 0.406
    beta  = 0.139
    gamma = 0.255
    delta = 1.761
    # The remnant mass over the NS rest mass
    remnant_mass = (max(alpha/eta**(1./3.)*(1 - 2*ns_compactness)
                        - beta*ns_compactness / eta*PG_ISSO_solver(chi, incl)
                        + gamma, 0.)) ** delta
    # Convert to solar masses
    remnant_mass = remnant_mass*ns_b_mass

    return remnant_mass


#############################################################################
# Vectorized version of remnant_mass. The third argument (NS equilibrium    #
# sequence) is excluded from vectorisation.                                 #
#############################################################################
remnant_masses = np.vectorize(remnant_mass)
remnant_masses.excluded.add(2)
