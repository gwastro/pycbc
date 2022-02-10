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
import scipy.optimize

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
