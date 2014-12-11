# Copyright (C) 2014 Francesco Pannarale
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
import os.path
import scipy.optimize
import scipy.interpolate
#from lal import LAL_PI, LAL_MTSUN_SI

# Setup the directory with the NS equilibrium sequence(s)
NS_SEQUENCE_FILE_DIRECTORY = os.path.join(os.path.dirname(__file__), 'ns_sequences')

#############################################################################
# Innermost Stable Spherical Orbit (ISSO) solver in the Perez-Giz formalism #
# [see Stone, Loeb, Berger, PRD 87, 084053 (2013)].                         #
#############################################################################

# Equation that determines the ISCO radius (in BH mass units)
def ISCO_eq(r, chi):
    return (r*(r-6))**2-chi**2*(2*r*(3*r+14)-9*chi**2)

# Equation that determines the ISCO radius (in BH mass units):
# arguments and sign are inverted wrt to ISCO_eq
def ISCO_eq_chi_first(chi,r):
    return -((r*(r-6))**2-chi**2*(2*r*(3*r+14)-9*chi**2))

# Equation that determines the ISSO radius (in BH mass units) at one of the
# poles 
def ISSO_eq_at_pole(r, chi):
    return r**3*(r**2*(r-6)+chi**2*(3*r+4))+chi**4*(3*r*(r-2)+chi**2)

# Equation that determines the ISSO radius (in BH mass units) for a generic
# orbital inclination 
def PG_ISSO_eq(r, chi, ci):
    Z=(r*(r-6))**2-chi**2*(2*r*(3*r+14)-9*chi**2)
    X=chi**2*(chi**2*(3*chi**2+4*r*(2*r-3))+r**2*(15*r*(r-4)+28))-6*r**4*(r**2-4)
    Y=chi**4*(chi**4+r**2*(7*r*(3*r-4)+36))+6*r*(r-2)*(chi**6+2*r**3*(chi**2*(3*r+2)+3*r**2*(r-2)))

    return r**8*Z+chi**2*(1-ci**2)*(chi**2*(1-ci**2)*Y-2*r**4*X)

# ISSO radius solver
def PG_ISSO_solver(chi,incl):
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
    #rISCO_limit = scipy.optimize.root(ISCO_eq, initial_guess, args=sgnchi, method='hybr').x[0]
    rISCO_limit = scipy.optimize.fsolve(ISCO_eq, initial_guess, args=sgnchi)

    # ISSO radius for an inclination of pi/2
    if chi < 0:
        initial_guess = 9
    else:
        initial_guess = 6
    #rISSO_at_pole_limit = scipy.optimize.root(ISSO_eq_at_pole, initial_guess, args=chi, method='hybr').x[0]
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
        #solution = scipy.optimize.root(PG_ISSO_eq, initial_guess, args=(chi, ci), method='hybr').x[0]
        solution = scipy.optimize.fsolve(PG_ISSO_eq, initial_guess, args=(chi, ci))
        if solution < 1 or solution > 9:
            initial_guess = min(rISCO_limit,rISSO_at_pole_limit)
            #solution = scipy.optimize.root(PG_ISSO_eq, initial_guess, args=(chi, ci), method='hybr').x[0]
            solution = scipy.optimize.fsolve(PG_ISSO_eq, initial_guess, args=(chi, ci))

    return solution

############################################################################
# Effective aligned spin of a NS-BH binary with tilted BH spin, as defined #
# in Stone, Loeb, Berger, PRD 87, 084053 (2013)].                          #
############################################################################
# Branch of positive chi solutions to rISCO(chi_eff) = rISSO(chi,incl) 
def pos_branch(incl, chi):
    if incl == 0:
        chi_eff = chi
    else:
        rISSO = PG_ISSO_solver(chi,incl) 
        #chi_eff = scipy.optimize.root(ISCO_eq_chi_first, 1.0, args=(rISSO),method='hybr').x[0]
        chi_eff = scipy.optimize.fsolve(ISCO_eq_chi_first, 1.0, args=(rISSO))

    return chi_eff

# Solution to rISCO(chi_eff) = rISSO(chi,incl) with the correct sign
def bh_effective_spin(chi,incl):
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
        #chi_eff = scipy.optimize.root(ISCO_eq_chi_first, initial_guess, args=(rISSO),method='hybr').x[0]
        chi_eff = scipy.optimize.fsolve(ISCO_eq_chi_first, initial_guess, args=(rISSO))

    return chi_eff

########################################################################################
# 2H 2-piecewise polytropic EOS, NS non-rotating equilibrium sequence                  #
# File format is: grav mass (Msun)   baryonic mass (Msun)    compactness               #
#                                                                                      # 
# Eventually, the script should call an NS sequence generator within LAL with the EOS  #
# prescribed by the user and store it.                                                 #
########################################################################################
def load_ns_sequence(eos_name):
    ns_sequence = []

    if eos_name == '2H':
        ns_sequence_path = os.path.join(NS_SEQUENCE_FILE_DIRECTORY, 'equil_2H.dat')
        ns_sequence = np.loadtxt(ns_sequence_path)
    else:
        print 'Only the 2H EOS is currently supported!'
        print 'If you plan to use a different NS EOS, be sure not to filter'
        print 'too many templates!\n'
        raise Exception('Unsupported EOS!')

    max_ns_g_mass = max(ns_sequence[:,0])

    return (ns_sequence, max_ns_g_mass)

####################################################################################
# Given an NS equilibrium sequence and gravitational mass (in Msun), return the NS #
# baryonic mass (in Msun).                                                         #
####################################################################################
def ns_g_mass_to_ns_b_mass(ns_g_mass, ns_sequence):
    x = ns_sequence[:,0] 
    y = ns_sequence[:,1]
    f = scipy.interpolate.interp1d(x, y)

    return f(ns_g_mass)

####################################################################################
# Given an NS equilibrium sequence and gravitational mass (in Msun), return the NS #
# compactness.                                                                     #
####################################################################################
def ns_g_mass_to_ns_compactness(ns_g_mass, ns_sequence):
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
# THIS ASSUMES THE NS SPIN IS 0. The user should be warned about this!                 #
########################################################################################
def xi_eq(x, kappa, chi_eff, q):
   return x**3*(x**2-3*kappa*x+2*chi_eff*kappa*math.sqrt(kappa*x))-3*q*(x**2-2*kappa*x+(chi_eff*kappa)**2)

def remnant_mass(eta, ns_g_mass, ns_sequence, chi, incl, shift):
    #Sanity checks
    if not (eta>0. and eta<=0.25 and abs(chi)<=1):
        print 'The BH spin magnitude must be <=1 and eta must be between 0 and 0.25'
        print 'This script was launched with ns_mass={0}, eta={1}, chi={2}, inclination={3}\n'.format(ns_b_mass, eta, chi, incl) 
        raise Exception('Unphysical parameters!')

    # Binary mass ratio define to be > 1
    q = (1+math.sqrt(1-4*eta)-2*eta)/eta*0.5

    # NS compactness and rest mass
    ns_compactness = ns_g_mass_to_ns_compactness(ns_g_mass, ns_sequence) 
    ns_b_mass   = ns_g_mass_to_ns_b_mass(ns_g_mass, ns_sequence)

    # Sanity checks
    if not (ns_compactness>0 and q>=1):
        print 'A positive NS compactness and a mass ratio that is >1 must be obtained.'
        print 'This script was launched with ns_mass={0}, eta={1}, chi={2}, inclination={3}'.format(ns_b_mass, eta, chi, incl) 
        print 'and obtained ns_compactness={0} and q={1}.'.format(ns_compactness, q)
        print 'SOMETHING WENT WRONG!!\n'
        raise Exception('Unphysical parameters!')

    # Calculate the dimensionless parameter kappa
    kappa = q*ns_compactness

    # Effective equatorial spin parameter needed to determine the torus mass*)
    chi_eff = bh_effective_spin(chi, incl)

    #Sanity checks
    if not abs(chi_eff)<=1:
        print 'The effective BH spin magnitude must be <=1'
        print 'This script was launched with ns_mass={0}, eta={1}, chi={2}, inclination={3}'.format(ns_b_mass, eta, chi, incl) 
        print 'and obtained chi_eff={0}.'.format(chi_eff)
        print 'SOMETHING WENT WRONG!!\n'
        raise Exception('Unphysical parameters!')

    # Fit parameters and tidal correction
    alpha = 0.296 # +/- 0.011
    beta  = 0.171 # +/- 0.008
    #xi = scipy.optimize.root(xi_eq, 100., args=(kappa,chi_eff,q), method='hybr').x[0]
    xi = scipy.optimize.fsolve(xi_eq, 100., args=(kappa,chi_eff,q))

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
    # Sanity checks
    if not (eta > 0. and eta <=0.25 and abs(bh_spin_z)<=1):
        print 'The absolute value of the BH spin z-component must be <=1.'
        print 'Eta must be between 0 and 0.25.'
        print 'This script was launched with eta={0} and chi_z={1}\n'.format(eta, bh_spin_z) 
        raise Exception('Unphysical parameters!')
    # To maximise the remnant mass, allow for the BH spin magnitude to be maximum
    bh_spin_magnitude = 1.
    # Unreasonably large remnant disk mass
    default_remnant_mass = 100.
    if not ns_g_mass > max_ns_g_mass:
        bh_spin_inclination = np.arccos(bh_spin_z/bh_spin_magnitude) 
        remnant_mass = pycbc.tmpltbank.em_progenitors.remnant_mass(eta, ns_g_mass, ns_sequence, bh_spin_magnitude, bh_spin_inclination, shift)
    else:
        remnant_mass = default_remnant_mass 

    return remnant_mass

################################################################################
# Given a NS mass, a BH spin z-component, and and EOS, find the minimum value  #
# of the symmetric mass ratio (eta) required to produce and EM counterpart.    #
# The user must specify the remnant disk mass threshold (shift) and a default  #
# value to be assigned to eta if the NS gravitational mass exceeds the maximum #
# NS mass allowed by the EOS (eta_default).i                                   #
################################################################################
def find_em_constraint_data_point(mNS, sBH, mBH_min, mBH_max, eos_name, shift, eta_default):
    ns_sequence, max_ns_g_mass = load_ns_sequence(eos_name)
    if mNS > max_ns_g_mass:
        eta_sol = eta_default
    else:
        eta_min = 0.01#mBH_max*mNS/(mBH_max+mNS)**2
        disk_mass_down = remnant_mass_ulim(eta_min, mNS, sBH, ns_sequence, max_ns_g_mass, shift) 
        eta_max = 0.25 #mBH_min*mNS/(mBH_min+mNS)**2
        disk_mass_up = remnant_mass_ulim(eta_max, mNS, sBH, ns_sequence, max_ns_g_mass, shift) 
        #print eta_min, disk_mass_down, eta_max, disk_mass_up
        if disk_mass_down*disk_mass_up < 0:
            try:
                # Initial guess in the middle of the range
                eta_sol = scipy.optimize.fsolve(remnant_mass_ulim, 0.5*(eta_max-eta_min), args=(mNS, sBH, ns_sequence, max_ns_g_mass, shift))
            except: 
                try:
                    # Initial guess is low
                    print 'Problems with mNS={0}, sBH,z={1}.: trying a different initial guess in the root-finder'.format(mNS, sBH)
                    eta_sol = scipy.optimize.fsolve(remnant_mass_ulim, eta_min+0.01, args=(mNS, sBH, ns_sequence, max_ns_g_mass, shift))
                except: 
                    # Initial guess is high
                    print 'Problems with mNS={0}, sBH,z={1}.: trying a different initial guess in the root-finder'.format(mNS, sBH)
                    eta_sol = scipy.optimize.fsolve(remnant_mass_ulim, eta_max-0.00001, args=(mNS, sBH, ns_sequence, max_ns_g_mass, shift))
        elif disk_mass_down > 0:
            eta_sol = 0.        # EM counterpart requires eta<eta_min: penalize this point 
        elif disk_mass_up < 0:
            eta_sol = 0.2500001 # EM counterpart is impossible
        elif disk_mass_down == 0:
            eta_sol = eta_min
        else:
            eta_sol = eta_max
        print mNS, sBH, eta_sol

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
def generate_em_constraint_data(mNS_min, mNS_max, delta_mNS, sBH_min, sBH_max, delta_sBH, mBH_min, mBH_max, eos_name, shift, eta_default): 
    import warnings
    warnings.filterwarnings('ignore', 'The iteration is not making good progress')

    # Build a grid of points in the mNS x sBHz space,
    # making sure maxima and minima are included
    mNS_nsamples = complex(0,int(np.ceil((mNS_max-mNS_min)/delta_mNS)+1))
    sBH_nsamples = complex(0,int(np.ceil((sBH_max-sBH_min)/delta_sBH)+1))
    mNS_vec, sBH_vec = np.mgrid[mNS_min:mNS_max:mNS_nsamples, sBH_min:sBH_max:sBH_nsamples]
    mNS_locations = np.array(mNS_vec[:,0])
    sBH_locations = np.array(sBH_vec[0])
    mNS_sBH_grid = zip(mNS_vec.ravel(), sBH_vec.ravel())
    mNS_sBH_grid = np.array(mNS_sBH_grid)
    mNS_vec = np.array(mNS_sBH_grid[:,0])
    sBH_vec = np.array(mNS_sBH_grid[:,1])

    # Until a numpy v>=1.7 is available everywhere, we have to used a silly
    # vectorization of Find_em_constraint_data_point and pass to it a bunch of
    # constant arguments as vectors with one entry repeated several times
    mBH_min_vec=np.empty(len(mNS_vec))
    mBH_min_vec.fill(mBH_min)
    mBH_max_vec=np.empty(len(mNS_vec))
    mBH_max_vec.fill(mBH_max)
    eos_name_vec=[eos_name for i in range(len(mNS_vec))]
    eos_name_vec=np.array(eos_name_vec)
    shift_vec=np.empty(len(mNS_vec))
    shift_vec.fill(shift)
    eta_default_vec=np.empty(len(mNS_vec))
    eta_default_vec.fill(eta_default)

    # Compute the minimum etas at all point in the mNS x sBHz grid
    eta_sol = find_em_constraint_data_points(mNS_vec, sBH_vec, mBH_min_vec, mBH_max_vec, eos_name_vec, shift_vec, eta_default_vec)
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
    f = scipy.interpolate.RectBivariateSpline(mNS_pts, sBH_pts, eta_mins, kx=1, ky=1)
    # If bh_spin_z a numpy array (assuming ns_g_mass has the same size)
    if isinstance(bh_spin_z, np.ndarray):
        eta_min = np.empty(len(bh_spin_z))
        for i in range(len(bh_spin_z)):
            eta_min[i] = f(ns_g_mass[i], bh_spin_z[i]) 
    # Else (assuming ns_g_mass and bh_spin_z are single numbers)
    else:
        eta_min = f(ns_g_mass, bh_spin_z)

    return eta_min

######################
# A few simple tests #
######################
# # TEST 1
# for i in np.arange(-1, 1.01, 0.1):
#     print i, PG_ISSO_solver(i,0), PG_ISSO_solver(i,math.pi/2)
# 
# # TEST 2: compare with Figure 4 of ApJ 668, 417 (2007)
# testfile = open('test.dat', 'w')
# testchis = [1, 0.998, 0.9, 0.5, 0, -0.5, -1]
# for chi in testchis:
#     for incl in np.arange(0, math.pi+0.001, math.pi/61):
#         testfile.write(str(PG_ISSO_solver(chi, incl))+' '+str(incl*180/math.pi)+'\n')
#     testfile.write('\n')
# 
# # TEST 3: compare with page 16 of 1212.4810 v1
# print bh_effective_spin(0.9,20*math.pi/180)
# print bh_effective_spin(0.9,40*math.pi/180)
# print bh_effective_spin(0.9,60*math.pi/180)
# 
# # TEST 4: Compute a remnant mass for a binary and EOS
# # 1) Set the properties of the binary 
# ns_g_mass = 1.4 # Msun
# eta = 0.16 # Dimensionless 
# bh_spin_magnitude = 0 # Dimensionless 
# bh_spin_inclination = 0 # Radians
# eos = '2H' # String
# # 2) Load the desired NS sequence table
# ns_sequence = load_ns_sequence(eos)
# # 3) Compute the remnant mass in solar mass units
# print 'The remnant mass for mNS={0}Msun, eta={1}, chi = {2}, i={3}, and the {4} EOS is {5}Msun'.format(ns_g_mass, eta, bh_spin_magnitude, bh_spin_inclination, eos, remnant_mass(ns_g_mass, eta, ns_sequence, bh_spin_magnitude, bh_spin_inclination, shift))


#ns_g_mass = 2.31890951636 # Msun
#bh_g_mass = 13.6596155569 # Msun
#eta = ns_g_mass*bh_g_mass/(bh_g_mass+ns_g_mass)**2 # Dimensionless 
#bh_spin_magnitude = 0.919209381906 # Dimensionless 
#bh_spin_inclination = 0 # Radians
#eos = '2H' # String
#ns_sequence, max_ns_g_mass = load_ns_sequence(eos)
#print 'The maximum NS mass is {0}Msun'.format(max_ns_g_mass)
#print 'The remnant mass for mNS={0}Msun, eta={1}, chi = {2}, i={3}, and the {4} EOS is:'.format(ns_g_mass, eta, bh_spin_magnitude, bh_spin_inclination, eos)
#print '     {0}Msun'.format(remnant_mass(eta, ns_g_mass, ns_sequence, bh_spin_magnitude, bh_spin_inclination, shift))
