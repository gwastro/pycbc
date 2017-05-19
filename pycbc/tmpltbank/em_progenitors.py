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

from __future__ import division
import math
import numpy as np
import pycbc
import pycbc.tmpltbank.em_progenitors
#from pycbc.tmpltbank import NS_SEQUENCE_FILE_DIRECTORY
import os.path
import scipy.optimize
import scipy.interpolate
#from lal import LAL_PI, LAL_MTSUN_SI

#############################################################################
# Innermost Stable Spherical Orbit (ISSO) solver in the Perez-Giz formalism #
# [see Stone, Loeb, Berger, PRD 87, 084053 (2013)].                         #
#############################################################################

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

# Equation that determines the ISCO radius (in BH mass units):
# arguments and sign are inverted wrt to ISCO_eq
def ISCO_eq_chi_first(chi,r):
    """
    Polynomial that enables the calculation of the Kerr
    inntermost stable circular orbit (ISCO) radius via its
    roots.  The arguments of the function and the sign of
    the polynomial are inverted with respect to ISCO_eq:
    this facilitates the job of the root-finder that calls
    this function.

    Parameters
    -----------
    chi: float
        the BH dimensionless spin parameter
    r: float
        the radial coordinate in BH mass units

    Returns
    ----------
    float
        -((r*(r-6))**2-chi**2*(2*r*(3*r+14)-9*chi**2))
    """
    return -ISCO_eq(r, chi)

# Equation that determines the ISSO radius (in BH mass units) at one of the
# poles
def ISSO_eq_at_pole(r, chi):
    """
    Polynomial that enables the calculation of the Kerr polar
    (inclination = +/- pi/2) innermost stable spherical orbit
    (ISSO) radius via its roots.  Physical solutions are
    between 6 and 1+sqrt[3]+sqrt[3+2sqrt[3]].

    Parameters
    -----------
    r: float
        the radial coordinate in BH mass units
    chi: float
        the BH dimensionless spin parameter

    Returns
    ----------
    float
        r**3*(r**2*(r-6)+chi**2*(3*r+4))+chi**4*(3*r*(r-2)+chi**2)
    """
    return r**3*(r**2*(r-6)+chi**2*(3*r+4))+chi**4*(3*r*(r-2)+chi**2)

# Equation that determines the ISSO radius (in BH mass units) for a generic
# orbital inclination
def PG_ISSO_eq(r, chi, ci):
    """
    Polynomial that enables the calculation of a generic innermost
    stable spherical orbit (ISSO) radius via its roots.  Physical
    solutions are between the equatorial ISSO (aka the ISCO) radius
    and the polar ISSO radius.
    [See Stone, Loeb, Berger, PRD 87, 084053 (2013)]

    Parameters
    -----------
    r: float
        the radial coordinate in BH mass units
    chi: float
        the BH dimensionless spin parameter
    ci: float
        cosine of the inclination angle between the BH spin
        and the orbital angular momentum

    Returns
    ----------
    float
        r**8*Z+chi**2*(1-ci**2)*(chi**2*(1-ci**2)*Y-2*r**4*X)
        where
        X=chi**2*(chi**2*(3*chi**2+4*r*(2*r-3))+r**2*(15*r*(r-4)+28))-6*r**4*(r**2-4)
        Y=chi**4*(chi**4+r**2*(7*r*(3*r-4)+36))+6*r*(r-2)*(chi**6+2*r**3*(chi**2*(3*r+2)+3*r**2*(r-2)))
        Z=ISCO_eq(r,chi)
    """
    X=chi**2*(chi**2*(3*chi**2+4*r*(2*r-3))+r**2*(15*r*(r-4)+28))-6*r**4*(r**2-4)
    Y=chi**4*(chi**4+r**2*(7*r*(3*r-4)+36))+6*r*(r-2)*(chi**6+2*r**3*(chi**2*(3*r+2)+3*r**2*(r-2)))
    Z=ISCO_eq(r, chi)

    return r**8*Z+chi**2*(1-ci**2)*(chi**2*(1-ci**2)*Y-2*r**4*X)

# ISSO radius solver
def PG_ISSO_solver(chi,incl):
    """
    Function that determines the radius of the innermost stable
    spherical orbit (ISSO) for a Kerr BH and a generic inclination
    angle between the BH spin and the orbital angular momentum.
    This function finds the appropriat root of PG_ISSO_eq.

    Parameters
    -----------
    chi: float
        the BH dimensionless spin parameter
    incl: float
        the inclination angle between the BH spin and the orbital
        angular momentum in radians

    Returns
    ----------
    solution: float
        the radius of the orbit in BH mass units
    """
    # Auxiliary variables
    ci=math.cos(incl)
    sgnchi = np.sign(ci)*chi

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
    rISSO_at_pole_limit = scipy.optimize.fsolve(ISSO_eq_at_pole, initial_guess, args=chi)

    # If the inclination is 0 or pi, just output the ISCO radius
    if incl in [0, math.pi]:
        solution = rISCO_limit
    # If the inclination is pi/2, just output the ISSO radius at the pole(s)
    elif incl == math.pi/2:
        solution = rISSO_at_pole_limit
    # Otherwise, find the ISSO radius for a generic inclination
    else:
        initial_guess = max(rISCO_limit,rISSO_at_pole_limit)
        solution = scipy.optimize.fsolve(PG_ISSO_eq, initial_guess, args=(chi, ci))
        if solution < 1 or solution > 9:
            initial_guess = min(rISCO_limit,rISSO_at_pole_limit)
            solution = scipy.optimize.fsolve(PG_ISSO_eq, initial_guess, args=(chi, ci))

    return solution

############################################################################
# Effective aligned spin of a NS-BH binary with tilted BH spin, as defined #
# in Stone, Loeb, Berger, PRD 87, 084053 (2013)].                          #
############################################################################
# Branch of positive chi solutions to rISCO(chi_eff) = rISSO(chi,incl)
def pos_branch(incl, chi):
    """
    Determines the effective [as defined in Stone, Loeb,
    Berger, PRD 87, 084053 (2013)] aligned dimensionless
    spin parameter of a NS-BH binary with tilted BH spin.
    This means finding the root chi_eff of
    ISCO_eq_chi_first(chi_eff, PG_ISSO_solver(chi,incl)).
    The result returned by this function belongs to the
    branch of the greater solutions, i.e. the greater of
    the two possible solutions is returned.

    Parameters
    -----------
    incl: float
        the inclination angle between the BH spin and the orbital
        angular momentum in radians
    chi: float
        the BH dimensionless spin parameter

    Returns
    ----------
    chi_eff: float
        the (greater) effective dimensionless spin parameter solution
    """
    if incl == 0:
        chi_eff = chi
    else:
        rISSO = PG_ISSO_solver(chi,incl)
        chi_eff = scipy.optimize.fsolve(ISCO_eq_chi_first, 1.0, args=(rISSO))

    return chi_eff

# Solution to rISCO(chi_eff) = rISSO(chi,incl) with the correct sign
def bh_effective_spin(chi,incl):
    """
    Determines the effective [as defined in Stone, Loeb,
    Berger, PRD 87, 084053 (2013)] aligned dimensionless
    spin parameter of a NS-BH binary with tilted BH spin.
    This means finding the root chi_eff of
    ISCO_eq_chi_first(chi_eff, PG_ISSO_solver(chi,incl))
    with the correct sign.

    Parameters
    -----------
    chi: float
        the BH dimensionless spin parameter
    incl: float
        the inclination angle between the BH spin and the orbital
        angular momentum in radians

    Returns
    ----------
    chi_eff: float
        the effective dimensionless spin parameter solution
    """
    if incl == 0:
        chi_eff = chi
    else:
        # ISSO radius for the given BH spin magnitude and inclination
        rISSO = PG_ISSO_solver(chi,incl)
        # Angle at which the branch of positive solutions has its minumum
        incl_flip = scipy.optimize.fmin(pos_branch, math.pi/4, args=tuple([chi]), full_output=False, disp=False)[-1]
        # Use incl_flip to determine the initial guess: the sign difference
        # in the initial_guess ensures that chi_eff has the correct sign
        if incl>incl_flip:
            initial_guess = -1.1
        else:
            initial_guess = 1.0
        chi_eff = scipy.optimize.fsolve(ISCO_eq_chi_first, initial_guess, args=(rISSO))

    return chi_eff

############################################################################
# 2H 2-piecewise polytropic EOS, NS non-rotating equilibrium sequence      #
# File format is: grav mass (Msun)   baryonic mass (Msun)    compactness   #
#                                                                          #
# Eventually, the function should call an NS sequence generator within LAL #
# the EOS prescribed by the user and store it.                             #
############################################################################
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
        #ns_sequence_path = os.path.join(NS_SEQUENCE_FILE_DIRECTORY, 'equil_2H.dat')
        ns_sequence = np.loadtxt(ns_sequence_path)
    else:
        print('Only the 2H EOS is currently supported!')
        print('If you plan to use a different NS EOS, be sure not to filter')
        print('too many templates!\n')
        raise Exception('Unsupported EOS!')

    max_ns_g_mass = max(ns_sequence[:,0])

    return (ns_sequence, max_ns_g_mass)

####################################################################################
# Given an NS equilibrium sequence and gravitational mass (in Msun), return the NS #
# baryonic mass (in Msun).                                                         #
####################################################################################
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
        The NS baryonic mass (in solar massesr**3*(r**2*(r-6)+chi**2*(3*r+4))+chi**4*(3*r*(r-2)+chi**2))
    """
    x = ns_sequence[:,0]
    y = ns_sequence[:,1]
    f = scipy.interpolate.interp1d(x, y)

    return f(ns_g_mass)

####################################################################################
# Given an NS equilibrium sequence and gravitational mass (in Msun), return the NS #
# compactness.                                                                     #
####################################################################################
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

########################################################################################
# Physical parameters passed to remnant_mass (see later) must not be unphysical        #
########################################################################################

########################################################################################
# Remnant mass for a NS-BH merger [Foucart PRD 86, 124007 (2012)]                      #
# The result is shifted by a quantity (shift, in solar masses) passed as an argument   #
# of remnant_mass: this can effectively used as a remnant mass threshold when solving  #
# the constraint remnant_mass(...)=0.                                                  #
# Allowing for negative remanant mass to be able to solve remnant mass == 0 if neeeded #
# THIS ASSUMES THE NS SPIN IS 0 (the user is warned about this)                        #
########################################################################################
def xi_eq(x, kappa, chi_eff, q):
    """
    The roots of this equation determine the orbital radius
    at the onset of NS tidal disruption in a nonprecessing
    NS-BH binary [(7) in Foucart PRD 86, 124007 (2012)]

    Parameters
    -----------
    x: float
        orbital separation in units of the NS radius
    kappa: float
        the BH mass divided by the NS radius
    chi_eff: float
        the BH dimensionless spin parameter
    q: float
        the binary mass ratio (BH mass / NS mass)

    Returns
    ----------
    float
        x**3*(x**2-3*kappa*x+2*chi_eff*kappa*sqrt[kappa*x)
        -3*q*(x**2-2*kappa*x+(chi_eff*kappa)**2)
    """
    return x**3*(x**2-3*kappa*x+2*chi_eff*kappa*math.sqrt(kappa*x))-3*q*(x**2-2*kappa*x+(chi_eff*kappa)**2)

def remnant_mass(eta, ns_g_mass, ns_sequence, chi, incl, shift):
    """
    Function that determines the remnant disk mass of
    an NS-BH system using the fit to numerical-relativity
    results discussed in Foucart PRD 86, 124007 (2012).

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
    shift: float
        an amount to be subtracted to the remnant mass predicted
        by the model (in solar masses)

    Returns
    ----------
    remnant_mass: float
        The remnant mass in solar masses
    """
    # Sanity checks
    if not (eta>0. and eta<=0.25 and abs(chi)<=1):
        print('The BH spin magnitude must be <=1 and eta must be between 0 and 0.25')
        print('This script was launched with ns_mass={0}, eta={1}, chi={2}, inclination={3}\n'.format(ns_g_mass, eta, chi, incl))
        raise Exception('Unphysical parameters!')

    # Binary mass ratio define to be > 1
    q = (1+math.sqrt(1-4*eta)-2*eta)/eta*0.5

    # NS compactness and rest mass
    ns_compactness = ns_g_mass_to_ns_compactness(ns_g_mass, ns_sequence)
    ns_b_mass = ns_g_mass_to_ns_b_mass(ns_g_mass, ns_sequence)

    # Sanity checks
    if not (ns_compactness>0 and q>=1):
        print('A positive NS compactness and a mass ratio that is >1 must be obtained.')
        print('This script was launched with ns_mass={0}, eta={1}, chi={2}, inclination={3}'.format(ns_b_mass, eta, chi, incl))
        print('and obtained ns_compactness={0} and q={1}.'.format(ns_compactness, q))
        print('SOMETHING WENT WRONG!!\n')
        raise Exception('Unphysical parameters!')

    # Calculate the dimensionless parameter kappa
    kappa = q*ns_compactness

    # Effective equatorial spin parameter needed to determine the torus mass*)
    chi_eff = bh_effective_spin(chi, incl)

    #Sanity checks
    if not abs(chi_eff)<=1:
        print('The effective BH spin magnitude must be <=1')
        print('This script was launched with ns_mass={0}, eta={1}, chi={2}, inclination={3}'.format(ns_b_mass, eta, chi, incl))
        print('and obtained chi_eff={0}.'.format(chi_eff))
        print('SOMETHING WENT WRONG!!\n')
        raise Exception('Unphysical parameters!')

    # Taking the 1st element with full_output=1 avoids some annoying messages on stdout
    xi = scipy.optimize.fsolve(xi_eq, 100., args=(kappa,chi_eff,q), full_output=1)[0]

    # Fit parameters and tidal correction
    alpha = 0.296 # +/- 0.011
    beta  = 0.171 # +/- 0.008
    # The remnant mass over the NS rest mass
    remnant_mass = alpha*xi*(1-2*ns_compactness)-beta*kappa*PG_ISSO_solver(chi_eff,0)

    # The remnant mass in the same units as the NS rest mass (presumably solar masses)
    remnant_mass = remnant_mass*ns_b_mass - shift

    return remnant_mass

######################################################################################
# Calculate an upper limit to the remnant mass by setting the BH spin magnitude to 1 #
# (its maximum possible value) and allowing the BH spin vector to be tilted so that  #
# chi_z*cos(tilt) = 1.  An unreasonably large remnant disk mass is returned if the   #
# maximum possible NS mass is exceeded.  Works with list and single numbers.         #
######################################################################################
def remnant_mass_ulim(eta, ns_g_mass, bh_spin_z, ns_sequence, max_ns_g_mass, shift):
    """
    Function that determines the maximum remnant disk mass
    for an NS-BH system with given symmetric mass ratio,
    NS mass, and BH spin parameter component along the
    orbital angular momentum.  This is a wrapper to
    the function remnant_mass.  Maximization is achieved
    by setting the BH dimensionless spin magntitude to unity.
    An unreasonably large remnant disk mass (100 solar masses)
    is returned if the maximum possible NS mass is exceeded
    in applying the model of Foucart PRD 86, 124007 (2012).

    Parameters
    -----------
    eta: float
        the symmetric mass ratio of the binary
    ns_g_mass: float
        NS gravitational mass (in solar masses)
    bh_spin_z: float
        the BH dimensionless spin parameter for the spin projection
        along the orbital angular momentum
    ns_sequence: 3D-array
        contains the sequence data in the form NS gravitational
         mass (in solar masses), NS baryonic mass (in solar
         masses), NS compactness (dimensionless)
    shift: float
        an amount to be subtracted to the remnant mass upper limit
        predicted by the model (in solar masses)

    Returns
    ----------
    remnant_mass_upper_limit: float
        The remnant mass upper limit in solar masses
    """
    # Sanity checks
    if not (eta > 0. and eta <=0.25 and abs(bh_spin_z)<=1):
        raise Exception("""The absolute value of the BH spin z-component must be <=1.
           Eta must be between 0 and 0.25.
           The function remnant_mass_ulim was launched with eta={0} and chi_z={1}.
           Unphysical parameters!""".format(eta, bh_spin_z))
    # To maximise the remnant mass, allow for the BH spin magnitude to be maximum
    bh_spin_magnitude = 1.
    # Unreasonably large remnant disk mass
    default_remnant_mass = 100.
    if not ns_g_mass > max_ns_g_mass:
        bh_spin_inclination = np.arccos(bh_spin_z/bh_spin_magnitude)
        remnant_mass_upper_limit = pycbc.tmpltbank.em_progenitors.remnant_mass(eta, ns_g_mass, ns_sequence, bh_spin_magnitude, bh_spin_inclination, shift)
    else:
        remnant_mass_upper_limit = default_remnant_mass

    return remnant_mass_upper_limit

################################################################################
# Given a NS mass, a BH spin z-component, and and EOS, find the minimum value  #
# of the symmetric mass ratio (eta) required to produce and EM counterpart.    #
# The user must specify the remnant disk mass threshold (thershold) and a      #
# default value to be assigned to eta if the NS gravitational mass exceeds the #
# maximum NS mass allowed by the EOS (eta_default).                            #
################################################################################
def find_em_constraint_data_point(mNS, sBH, eos_name, threshold, eta_default):
    """
    Function that determines the minimum symmetric mass ratio
    for an NS-BH system with given NS mass, BH spin parameter
    component along the orbital angular momentum, and NS equation
    of state (EOS), required for the remnant disk mass to exceed
    a certain threshold value specified by the user.  A default
    value specified by the user is returned if the NS gravitational
    mass exceeds the maximum NS mass allowed by the EOS.

    Parameters
    -----------
    mNS: float
        the NS mass
    sBH: float
        BH dimensionless spin parameter for the BH spin
        component along the orbital angular momentum
    eos_name: string
        NS equation of state label ('2H' is the only supported
        choice at the moment)
    eta_default: float
        the value to be returned if the input NS mass is too high
    threshold: float
        an amount to be subtracted to the remnant mass upper limit
        predicted by the model (in solar masses)

    Returns
    ----------
    eta_sol: float
        The miniumum symmetric mass ratio for which the remnant
        disk mass exceeds the threshold
    """
    ns_sequence, max_ns_g_mass = load_ns_sequence(eos_name)
    if mNS > max_ns_g_mass:
        eta_sol = eta_default
    else:
        eta_min = 0.01 #mBH_max*mNS/(mBH_max+mNS)**2
        disk_mass_down = remnant_mass_ulim(eta_min, mNS, sBH, ns_sequence, max_ns_g_mass, threshold)
        eta_max = 0.25 #mBH_min*mNS/(mBH_min+mNS)**2
        disk_mass_up = remnant_mass_ulim(eta_max, mNS, sBH, ns_sequence, max_ns_g_mass, threshold)
        if disk_mass_down*disk_mass_up < 0:
            # Methods that work are (in order of performance speed): brentq, brenth, ridder, bisect
            eta_sol = scipy.optimize.brentq(remnant_mass_ulim, eta_min, eta_max, args=(mNS, sBH, ns_sequence, max_ns_g_mass, threshold))
        elif disk_mass_down > 0:
            eta_sol = 0.        # EM counterpart requires eta<eta_min: penalize this point
        elif disk_mass_up < 0:
            eta_sol = 0.2500001 # EM counterpart is impossible
        elif disk_mass_down == 0:
            eta_sol = eta_min
        else:
            eta_sol = eta_max

    return eta_sol

################################################################################
# Vectorized version of find_em_constraint_data_point.  Numpy v1.7 and above   #
# allows one to list the arguments to exclude from vectorization.  An older    #
# version of numpy is available on several clusters, so for now we will have   #
# to work around this and make sure we pass all arguments as vectors of the    #
# same length.                                                                 #
################################################################################
find_em_constraint_data_points = np.vectorize(find_em_constraint_data_point)

################################################################################
# Generate the (bh_spin_z, ns_g_mass, eta) surface above which EM counterparts #
# are possible.  The user must specify the remnant disk mass threshold (shift) #
# and the default eta to be assigned to points for which ns_g_mass exceeds the #
# maximum NS mass allowed by the EOS (eta_default).                            #
################################################################################
def generate_em_constraint_data(mNS_min, mNS_max, delta_mNS, sBH_min, sBH_max, delta_sBH, eos_name, threshold, eta_default):
    """
    Wrapper that calls find_em_constraint_data_point over a grid
    of points to generate the bh_spin_z x ns_g_mass x eta surface
    above which NS-BH binaries yield a remnant disk mass that
    exceeds the threshold required by the user.  The user must also
    specify the default symmetric mass ratio value to be assigned
    to points for which the NS mass exceeds the maximum NS mass
    allowed by the chosend NS equation of state.
    The 2D surface that is generated is saved to file in two formats:
    constraint_em_bright.npz and constraint_em_bright.npz.

    Parameters
    -----------
    mNS_min: float
        lower boundary of the grid in the NS mass direction
    mNS_max: float
        upper boundary of the grid in the NS mass direction
    delta_mNS: float
        grid spacing in the NS mass direction
    sBH_min: float
        lower boundary of the grid in the direction of the
        BH dimensionless spin component along the orbital
        angular momentum
    sBH_max: float
        upper boundary of the grid in the direction of the
        BH dimensionless spin component along the orbital
        angular momentum
    delta_sBH: float
        grid spacing in the direction of the BH dimensionless
        spin component along the orbital angular momentum
    eos_name: string
        NS equation of state label ('2H' is the only supported
        choice at the moment)
    threshold: float
        an amount to be subtracted to the remnant mass upper limit
        predicted by the model (in solar masses)
    eta_default: float
        the value to be returned for points in the grids in which
        the NS mass is too high
    """
    # Build a grid of points in the mNS x sBHz space,
    # making sure maxima and minima are included
    mNS_nsamples = complex(0,int(np.ceil((mNS_max-mNS_min)/delta_mNS)+1))
    sBH_nsamples = complex(0,int(np.ceil((sBH_max-sBH_min)/delta_sBH)+1))
    mNS_vec, sBH_vec = np.mgrid[mNS_min:mNS_max:mNS_nsamples, sBH_min:sBH_max:sBH_nsamples] # pylint:disable=invalid-slice-index
    mNS_locations = np.array(mNS_vec[:,0])
    sBH_locations = np.array(sBH_vec[0])
    mNS_sBH_grid = zip(mNS_vec.ravel(), sBH_vec.ravel())
    mNS_sBH_grid = np.array(mNS_sBH_grid)
    mNS_vec = np.array(mNS_sBH_grid[:,0])
    sBH_vec = np.array(mNS_sBH_grid[:,1])

    # Until a numpy v>=1.7 is available everywhere, we have to use a silly
    # vectorization of find_em_constraint_data_point and pass to it a bunch of
    # constant arguments as vectors with one entry repeated several times
    eos_name_vec=[eos_name for _ in range(len(mNS_vec))]
    eos_name_vec=np.array(eos_name_vec)
    threshold_vec=np.empty(len(mNS_vec))
    threshold_vec.fill(threshold)
    eta_default_vec=np.empty(len(mNS_vec))
    eta_default_vec.fill(eta_default)

    # Compute the minimum etas at all point in the mNS x sBHz grid
    eta_sol = find_em_constraint_data_points(mNS_vec, sBH_vec, eos_name_vec, threshold_vec, eta_default_vec)
    eta_sol = eta_sol.reshape(-1,len(sBH_locations))
    # Save the results
    np.savez('constraint_em_bright', mNS_pts=mNS_locations, sBH_pts=sBH_locations, eta_mins=eta_sol)
    # Cast the results in a format that is quick to plot from textfile
    constraint_data = zip(mNS_vec.ravel(), sBH_vec.ravel(), eta_sol.ravel())
    np.savetxt('constraint_em_bright.dat', constraint_data)

################################################################################
# Given a BH spin z-component and an NS mass, return the minumum symmetric     #
# mass ratio (eta) required to produce an EM counterpart.  The function uses   #
# data of a sampled constraint surface eta(bh_spin_z, ns_g_mass).  The remant  #
# mass threshold to generate a counterpart must be set when generating such    #
# constraint data (see generate_em_constraint_data).                           #
################################################################################
def min_eta_for_em_bright(bh_spin_z, ns_g_mass, mNS_pts, sBH_pts, eta_mins):
    """
    Function that uses the end product of generate_em_constraint_data
    to swipe over a set of NS-BH binaries and determine the minimum
    symmetric mass ratio required by each binary to yield a remnant
    disk mass that exceeds a certain threshold.  Each binary passed
    to this function consists of a NS mass and a BH spin parameter
    component along the orbital angular momentum.  Unlike
    find_em_constraint_data_point, which solves the problem at
    a given point in the paremter space and is more generic, this
    function interpolates the results produced by
    generate_em_constraint_data at the desired locations:
    generate_em_constraint_data must be run once prior to calling
    min_eta_for_em_bright.

    Parameters
    -----------
    bh_spin_z: array
        desired values of the BH dimensionless spin parameter for the
        spin projection along the orbital angular momentum
    ns_g_mass: array
        desired values of the NS gravitational mass (in solar masses)
    mNS_pts: array
        NS mass values (in solar masses) from the output of
        generate_em_constraint_data
    sBH_pts: array
        BH dimensionless spin parameter values along the orbital
        angular momentum from the output of generate_em_constraint_data
    eta_mins: array
        minimum symmetric mass ratio values to exceed a given remnant
        disk mass threshold from the output of generate_em_constraint_data

    Returns
    ----------
    eta_min: array
        the minimum symmetric mass ratio required by each binary in the
        input to yield a remnant disk mass that exceeds a certain
        threshold
    """
    f = scipy.interpolate.RectBivariateSpline(mNS_pts, sBH_pts, eta_mins, kx=1, ky=1)
    # If bh_spin_z is a numpy array (assuming ns_g_mass has the same size)
    if isinstance(bh_spin_z, np.ndarray):
        eta_min = np.empty(len(bh_spin_z))
        for i in range(len(bh_spin_z)):
            eta_min[i] = f(ns_g_mass[i], bh_spin_z[i])
    # Else (assuming ns_g_mass and bh_spin_z are single numbers)
    else:
        eta_min = f(ns_g_mass, bh_spin_z)

    return eta_min
