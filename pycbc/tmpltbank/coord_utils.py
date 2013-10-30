#!/usr/bin/env python

# Copyright (C) 2013 Ian W. Harry
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
import numpy
from lal import LAL_PI, LAL_MTSUN_SI
from pycbc.tmpltbank.lambda_mapping import get_chirp_params
from pycbc.tmpltbank.lambda_mapping import get_beta_sigma_from_aligned_spins

def estimate_mass_range(numPoints, order, evals, evecs, maxmass1, minmass1, \
                        maxmass2, minmass2, maxspin, f0, covary=True, \
                        evecsCV=None, maxBHspin=None, \
                        minTotalMass=None, maxTotalMass=None):
    """
    This function will generate a large set of points with random masses and
    spins (using pycbc.tmpltbank.get_random_mass) and translate these points
    into the xi_i coordinate system.

    Parameters
    ----------
    numPoints : int
        Number of systems to simulate
    order : string
        The Post-Newtonian order that is used to translate from masses and
        spins to the lambda_i coordinate system.
    evals : numpy.array
        Array of the eigenvalues of the metric in lambda_i coordinates. Used
        to rotate to a Cartesian coordinate system.
    evecs : numpy.matrix
        Matrix of the eigenvectors of the metric in lambda_i coordinates. Used
        to rotate to a Cartesian coordinate system.
    maxmass1 : float
        Maximum allowed mass of the first body.
    minmass1 : float
        Minimum allowed mass of the first body.
    maxmass2 : float
        Maximum allowed mass of the second body.
    minmass2 : float
        Minimum allowed mass of the second body.
    maxspin : float
        Maximum spin allowed on both bodies if maxBHspin is not given. If
        maxBHspin is given, this will be the maximum allowed spin for neutron
        stars only (mass < 3 solar masses).
    f0 : float
        This is an arbitrary scaling factor introduced to avoid the potential
        for numerical overflow when calculating this. Generally the default
        value (70) is safe here. **IMPORTANT, if you want to calculate the
        ethinca metric components later this MUST be set equal to f_low.**
        This value must also be used consistently (ie. don't change its value
        when calling different functions!).
    covary : boolean, optional (default = True)
        If this is given then evecsCV will be used to rotate from the Cartesian
        coordinate system into the principal coordinate direction (xi_i). If
        not given then points in the original Caretesian coordinates are 
        returned.
    evecsCV : numpy.matrix, optional (default = None)
        If covary=True this matrix is used to perform the rotation to the xi_i
        coordinate system.
    maxBHspin : float, optional (default = None)
        If this is given it will be the maximum allowed spin for black holes
        (mass > 3 solar masses). If not given maxspin will be used for both
        black hole and neutron star spin.
    minTotalMass : float, optional (default = None)
        If this float is given the total mass of the generated systems must be
        larger than this number.
    maxTotalMass : float, optional (default = None)
        If this float is given the total mass of the generated systems must be
        smaller than this number.

    Return
    -------
    xis : numpy.array
        A list of the positions of each point in the xi_i coordinate system.
    """
    valsF = get_random_mass(numPoints, minmass1, maxmass1, minmass2, \
          maxmass2, maxspin, maxBHspin=maxBHspin, minTotalMass=minTotalMass, \
          maxTotalMass=maxTotalMass)
    valsF = numpy.array(valsF)
    mass = valsF[0]
    eta = valsF[1]
    beta = valsF[2]
    sigma = valsF[3]
    gamma = valsF[4]
    chis = 0.5*(valsF[5] + valsF[6])
    if covary:
        lambdas = get_cov_params(mass, eta, beta, sigma, gamma, chis, f0, \
                                 evecs, evals, evecsCV, order)
    else:
        lambdas = get_conv_params(mass, eta, beta, sigma, gamma, chis, f0, \
                                  evecs, evals, order)

    return numpy.array(lambdas)

def get_random_mass(numPoints, minmass1, maxmass1, minmass2, \
                    maxmass2, maxspin, maxBHspin=None, \
                    minTotalMass=None, maxTotalMass=None, qm_scalar_fac=1):
    """
    This function will generate a large set of points within the chosen mass
    and spin space. It will also return the corresponding PN spin coefficients
    for ease of use later (though these may be removed at some future point).

    Parameters
    ----------
    numPoints : int
        Number of systems to simulate
    minmass1 : float
        Maximum allowed mass of the first body.
    maxmass1 : float
        Minimum allowed mass of the first body.
    minmass2 : float
        Maximum allowed mass of the second body.
    maxmass2 : float
        Minimum allowed mass of the second body.
    maxspin : float
        Maximum spin allowed on both bodies if maxBHspin is not given. If
        maxBHspin is given, this will be the maximum allowed spin for neutron
        stars only (mass < 3 solar masses).
    maxBHspin : float, optional (default = None)
        If this is given it will be the maximum allowed spin for black holes
        (mass > 3 solar masses). If not given maxspin will be used for both
        black hole and neutron star spin.
    minTotalMass : float, optional (default = None)
        If this float is given the total mass of the generated systems must be
        larger than this number.
    maxTotalMass : float, optional (default = None)
        If this float is given the total mass of the generated systems must be
        smaller than this number.
    qm_scalar_fac : float, optional (default = 1)
        The 2PN sigma spin coefficient has an equation-of-state dependent term.
        For black holes this value is equal to 1, for NS this value is larger,
        will depend on masses and is not known. If this value is given then
        it will be used when computing the sigma coefficient to model NS EOS
        effects on this 2PN coefficient.
        NOTE: This is currently a placeholder and does nothing. It needs to
        be implemented in pycbc.tmpltbank.get_beta_sigma_from_aligned_spins().

    Returns
    --------
    mass : numpy.array
        List of the total masses.
    eta : numpy.array
        List of the symmetric mass ratios
    beta : numpy.array
        List of the 1.5PN beta spin coefficients
    sigma : numpy.array
        List of the 2PN sigma spin coefficients
    gamma : numpy.array
        List of the 2.5PN gamma spin coefficients
    spin1z : numpy.array
        List of the spin on the heavier body. NOTE: Body 1 is **always** the
        heavier body to remove mass,eta -> m1,m2 degeneracy
    spin2z : numpy.array
        List of the spin on the smaller body. NOTE: Body 2 is **always** the
        smaller body to remove mass,eta -> m1,m2 degeneracy
    """

    # WARNING: We expect mass1 > mass2 ALWAYS
    # Check inputs, and set the minimum/maximum total masses
    minmass = minmass1 + minmass2
    maxmass = maxmass1 + maxmass2
    if minTotalMass and (minTotalMass > minmass):
        minmass = minTotalMass
    if maxTotalMass and (maxTotalMass < maxmass):
        maxmass = maxTotalMass

    if (minmass2 > minmass1) or (maxmass2 > maxmass1):
        errMsg = "Mass1 must be larger than mass2. Check input options."
        raise ValueError(errMsg)

    if (minmass2 > maxmass2) or (minmass1 > maxmass1):
        errMsg = "Minimum masses cannot be larger than maximum masses."
        errMsg += "Check input options."
        raise ValueError(errMsd)

    mincompmass = minmass2
    maxcompmass = maxmass1

    # First we choose the total masses from a unifrom distribution in mass
    # to the -5/3. power.
    mass = numpy.random.random(numPoints) * \
          (minmass**(-5./3.)-maxmass**(-5./3.)) + maxmass**(-5./3.)
    mass = mass**(-3./5.)

    # Next we choose the mass ratios, this will take different limits based on
    # the value of total mass
    maxmass2 = numpy.minimum(mass/2.,maxmass2)
    minmass1 = numpy.maximum(minmass1,mass/2.)
    mineta = numpy.maximum(mincompmass * (mass-mincompmass)/(mass*mass),\
                           maxcompmass*(mass-maxcompmass)/(mass*mass))
    maxeta = numpy.minimum(0.25,maxmass2 * (mass - maxmass2) / (mass*mass))
    maxeta = numpy.minimum(maxeta,minmass1 * (mass - minmass1) / (mass*mass))
    if (maxeta < mineta).any():
        errMsg = "WARNING: Max eta is smaller than min eta!!"
        raise ValueError(errMsg)     
    eta = numpy.random.random(numPoints) * (maxeta - mineta) + mineta

    # Also calculate the component masses; mass1 > mass2
    diff = (mass*mass * (1-4*eta))**0.5
    mass1 = (mass + diff)/2.
    mass2 = (mass - diff)/2.
    # Check the masses are where we want them to be (allowing some floating
    # point rounding error).
    if (mass1 > maxmass1*1.001).any() or (mass1 < minmass1*0.999).any():
        errMsg = "Mass1 is not within the specified mass range."
        raise ValueError(errMsg)
    if (mass2 > maxmass2*1.001).any() or (mass2 < minmass2*0.999).any():
        errMsg = "Mass2 is not within the specified mass range."
        raise ValueError(errMsg)

    # Next up is the spins. First check if we have non-zero spins
    if maxspin == 0 and not maxBHspin:
        spinspin = numpy.zeros(numPoints,dtype=float)
        spin1z = numpy.zeros(numPoints,dtype=float)
        spin2z = numpy.zeros(numPoints,dtype=float)
    else:
        # Spin 1 first
        mspin = numpy.zeros(len(mass1))
        mspin += maxspin
        if maxBHspin:
            mspin[mass1 > 3] = maxBHspin
        spin1z = (2*numpy.random.random(numPoints) - 1) * mspin
        # Then spin 2
        mspin = numpy.zeros(len(mass2))
        mspin += maxspin
        if maxBHspin:
            mspin[mass2 > 3] = maxBHspin
        spin2z = (2*numpy.random.random(numPoints) - 1) * mspin
        spinspin = spin1z*spin2z

    # And compute the PN components that come out of this
    beta, sigma, gamma, chiS = get_beta_sigma_from_aligned_spins(\
                                 mass, eta, spin1z, spin2z)

    return mass,eta,beta,sigma,gamma,spin1z,spin2z

def get_cov_params(totmass, eta, beta, sigma, gamma, chis, f0, evecs, evals, \
                   evecsCV, order):
    """
    Function to convert between masses and spins and locations in the xi
    parameter space. Xi = Cartesian metric and rotated to principal components.

    Parameters
    -----------
    totmass : float or numpy.array
        Total mass(es) of the system(s)
    eta : float or numpy.array
        Symmetric mass ratio(s) of the system(s)
    beta : float or numpy.array
        1.5PN spin coefficienct(s) of the system(s)
    sigma: float or numpy.array
        2PN spin coefficienct(s) of the system(s)
    gamma : float or numpy.array
        2.5PN spin coefficienct(s) of the system(s)
    chis : float or numpy.array
        0.5 * (spin1z + spin2z) for the system(s)
    f0 : float
        This is an arbitrary scaling factor introduced to avoid the potential
        for numerical overflow when calculating this. Generally the default
        value (70) is safe here. **IMPORTANT, if you want to calculate the
        ethinca metric components later this MUST be set equal to f_low.**
        This value must also be used consistently (ie. don't change its value
        when calling different functions!).
    order : string
        The Post-Newtonian order that is used to translate from masses and
        spins to the lambda_i coordinate system.
    evals : numpy.array
        Array of the eigenvalues of the metric in lambda_i coordinates. Used
        to rotate to a Cartesian coordinate system.
    evecs : numpy.matrix
        Matrix of the eigenvectors of the metric in lambda_i coordinates. Used
        to rotate to a Cartesian coordinate system.
    evecsCV : numpy.matrix, optional (default = None)
        If covary=True this matrix is used to perform the rotation to the xi_i
        coordinate system.

    Returns
    --------
    xis : list of floats or numpy.arrays
        Position of the system(s) in the xi coordinate system
    """

    # Do this by doing masses - > lambdas -> mus
    mus = get_conv_params(totmass, eta, beta, sigma, gamma, chis, f0, evecs, \
                          evals, order)
    # and then mus -> xis
    xis = get_covaried_params(mus, evecsCV)
    return xis

def get_conv_params(totmass, eta, beta, sigma, gamma,chis, f0, evecs, evals, \
                    order):
    """
    Function to convert between masses and spins and locations in the mu
    parameter space. Mu = Cartesian metric, but not principal components.

    Parameters
    -----------
    totmass : float or numpy.array
        Total mass(es) of the system(s)
    eta : float or numpy.array
        Symmetric mass ratio(s) of the system(s)
    beta : float or numpy.array
        1.5PN spin coefficienct(s) of the system(s)
    sigma: float or numpy.array
        2PN spin coefficienct(s) of the system(s)
    gamma : float or numpy.array
        2.5PN spin coefficienct(s) of the system(s)
    chis : float or numpy.array
        0.5 * (spin1z + spin2z) for the system(s)
    f0 : float
        This is an arbitrary scaling factor introduced to avoid the potential
        for numerical overflow when calculating this. Generally the default
        value (70) is safe here. **IMPORTANT, if you want to calculate the
        ethinca metric components later this MUST be set equal to f_low.**
        This value must also be used consistently (ie. don't change its value
        when calling different functions!).
    order : string
        The Post-Newtonian order that is used to translate from masses and
        spins to the lambda_i coordinate system.
    evals : numpy.array
        Array of the eigenvalues of the metric in lambda_i coordinates. Used
        to rotate to a Cartesian coordinate system.
    evecs : numpy.matrix
        Matrix of the eigenvectors of the metric in lambda_i coordinates. Used
        to rotate to a Cartesian coordinate system.
    evecsCV : numpy.matrix
        This matrix is used to perform the rotation to the xi_i
        coordinate system.

    Returns
    --------
    mus : list of floats or numpy.arrays
        Position of the system(s) in the mu coordinate system
    """

    # Do this by masses -> lambdas
    lambdas = get_chirp_params(totmass,eta,beta,sigma,gamma,chis,f0,order)
    # and lambdas -> mus
    mus = get_mu_params(lambdas, f0, evecs, evals)
    return mus

def get_mu_params(lambdas, f0, evecs, evals):
    """
    Function to rotate from the lambda coefficients into position in the mu
    coordinate system. Mu = Cartesian metric, but not principal components.

    Parameters
    -----------
    lambdas : list of floats or numpy.arrays
        Position of the system(s) in the lambda coefficients
    f0 : float
        This is an arbitrary scaling factor introduced to avoid the potential
        for numerical overflow when calculating this. Generally the default
        value (70) is safe here. **IMPORTANT, if you want to calculate the
        ethinca metric components later this MUST be set equal to f_low.**
        This value must also be used consistently (ie. don't change its value
        when calling different functions!).
    evals : numpy.array
        Array of the eigenvalues of the metric in lambda_i coordinates. Used
        to rotate to a Cartesian coordinate system.
    evecs : numpy.matrix
        Matrix of the eigenvectors of the metric in lambda_i coordinates. Used
        to rotate to a Cartesian coordinate system.

    Returns
    --------
    mus : list of floats or numpy.arrays
        Position of the system(s) in the mu coordinate system
    """

    mus = []
    for i in range(len(evals)):
        mus.append(rotate_vector(evecs,lambdas,numpy.sqrt(evals[i]),i))
    return mus

def get_covaried_params(mus,evecsCV):
    """
    Function to rotate from position(s) in the mu_i coordinate system into the
    position(s) in the xi_i coordinate system

    Parameters
    -----------
    mus : list of floats or numpy.arrays
        Position of the system(s) in the mu coordinate system
    evecsCV : numpy.matrix
        This matrix is used to perform the rotation to the xi_i
        coordinate system.

    Returns
    --------
    xis : list of floats or numpy.arrays
        Position of the system(s) in the xi coordinate system
    """
    xis = []
    for i in range(len(evecsCV)):
        xis.append(rotate_vector(evecsCV,mus,1.,i))
    return xis

def rotate_vector(evecs,old_vector,rescale_factor,index):
    """
    Function to find the position of the system(s) in one of the xi_i or mu_i
    directions. 

    Parameters
    -----------
    evecs : numpy.matrix
        Matrix of the eigenvectors of the metric in lambda_i coordinates. Used
        to rotate to a Cartesian coordinate system.
    old_vector : list of floats or numpy.arrays
        The position of the system(s) in the original coordinates
    rescale_factor : float
        Scaling factor to apply to resulting position(s)
    index : int
        The index of the final coordinate system that is being computed. Ie.
        if we are going from mu_i -> xi_j, this will give j.
    
    Returns
    --------
    positions : float or numpy.array
        Position of the point(s) in the resulting coordinate.
    """
    temp = 0
    for i in range(len(evecs)):
        temp += evecs[i,index] * old_vector[i]
    temp *= rescale_factor
    return temp

def get_point_distance(point1, point2, evals, evecs, evecsCV, order, f0):
    """
    Function to calculate the mismatch between two points, supplied in terms
    of the masses and spins, using the xi_i parameter space metric to
    approximate the mismatch of the two points. Can also take one of the points
    as an array of points and return an array of mismatches (but only one can
    be an array!)

    point1 : List of floats or numpy.arrays
        point1[0] contains the mass(es) of the heaviest body(ies).
        point1[1] contains the mass(es) of the smallest body(ies).
        point1[2] contains the spin(es) of the heaviest body(ies).
        point1[3] contains the spin(es) of the smallest body(ies).
    point2 : List of floats
        point2[0] contains the mass of the heaviest body.
        point2[1] contains the mass of the smallest body.
        point2[2] contains the spin of the heaviest body.
        point2[3] contains the spin of the smallest body.
    evals : numpy.array
        Array of the eigenvalues of the metric in lambda_i coordinates. Used
        to rotate to a Cartesian coordinate system.
    evecs : numpy.matrix
        Matrix of the eigenvectors of the metric in lambda_i coordinates. Used
        to rotate to a Cartesian coordinate system.
    evecsCV : numpy.matrix, optional (default = None)
        If covary=True this matrix is used to perform the rotation to the xi_i
        coordinate system.
    order : string
        The Post-Newtonian order that is used to translate from masses and
        spins to the lambda_i coordinate system.
    f0 : float
        This is an arbitrary scaling factor introduced to avoid the potential
        for numerical overflow when calculating this. Generally the default
        value (70) is safe here. **IMPORTANT, if you want to calculate the
        ethinca metric components later this MUST be set equal to f_low.**
        This value must also be used consistently (ie. don't change its value
        when calling different functions!).

    Returns
    --------
    dist : float or numpy.array
        Distance between the point2 and all points in point1
    xis1 : List of floats or numpy.arrays
        Position of the input point1(s) in the xi_i parameter space
    xis2 : List of floats
        Position of the input point2 in the xi_i parameter space
    """
    aMass1 = point1[0]
    aMass2 = point1[1]
    aSpin1 = point1[2]
    aSpin2 = point1[3]
    try:
        leng = len(aMass1)
        aArray = True
    except:
        aArray = False

    bMass1 = point2[0]
    bMass2 = point2[1]
    bSpin1 = point2[2]
    bSpin2 = point2[3]
    bArray = False

    aTotMass = aMass1 + aMass2
    aEta = (aMass1 * aMass2) / (aTotMass * aTotMass)
    aCM = aTotMass * aEta**(3./5.)

    bTotMass = bMass1 + bMass2
    bEta = (bMass1 * bMass2) / (bTotMass * bTotMass)
    bCM = bTotMass * bEta**(3./5.)

    abeta, asigma, agamma, achis = get_beta_sigma_from_aligned_spins(aTotMass,\
                                     aEta, aSpin1, aSpin2)
    bbeta, bsigma, bgamma, bchis = get_beta_sigma_from_aligned_spins(bTotMass,\
                                     bEta, bSpin1, bSpin2)

    aXis = get_cov_params(aTotMass, aEta, abeta, asigma, agamma, achis, f0,\
                          evecs, evals, evecsCV, order)

    bXis = get_cov_params(bTotMass, bEta, bbeta, bsigma, bgamma, bchis, f0,\
                          evecs, evals, evecsCV, order)

    dist = (aXis[0] - bXis[0])**2
    for i in range(1,len(aXis)):
        dist += (aXis[i] - bXis[i])**2

    return dist, aXis, bXis

def calc_point_dist(vsA, entryA, MMdistA):
    """
    This function is used to determine if the distance between two points is
    less than that stated by the minimal match.

    Parameters
    ----------
    vsA : list or numpy.array or similar
        An array of point 1's position in the \chi_i coordinate system
    entryA : list or numpy.array or similar
        An array of point 2's position in the \chi_i coordinate system
    MMdistA : float
        The minimal mismatch allowed between the points

    Returns
    --------
    Boolean
        True if the points have a mismatch < MMdistA
        False if the points have a mismatch > MMdistA
    """
    val = (vsA[0] - entryA[0])**2
    for i in range(1,len(vsA)):
        val += (vsA[i] - entryA[i])**2
    return (val < MMdistA)

def calc_point_dist_vary(mus1, fUpper1, mus2, fUpper2, fMap, MMdistA):
    """
    Function to determine if two points, with differing upper frequency cutoffs
    have a mismatch < MMdistA for *both* upper frequency cutoffs.

    Parameters
    ----------
    mus1 : List of numpy arrays
        mus1[i] will give the array of point 1's position in the \chi_j
        coordinate system. The i element corresponds to varying values of the
        upper frequency cutoff. fMap is used to map between i and actual
        frequencies
    fUpper1 : float
        The upper frequency cutoff (ISCO) of point 1.
    mus2 : List of numpy arrays
        mus2[i] will give the array of point 2's position in the \chi_j
        coordinate system. The i element corresponds to varying values of the
        upper frequency cutoff. fMap is used to map between i and actual
        frequencies
    fUpper2 : float
        The upper frequency cutoff (ISCO) of point 2.
    fMap : dictionary
        fMap[fUpper] will give the index needed to get the \chi_j coordinates
        in the two sets of mus
    MMdistA
        The minimal mismatch allowed between the points

    Returns
    --------
    Boolean
        True if the points have a mismatch < MMdistA
        False if the points have a mismatch > MMdistA
    """
    idx1 = fMap[fUpper1]
    vecs1 = mus1[idx1]
    vecs2 = mus2[idx1]
    val = (vecs1[0] - vecs2[0])**2
    for i in range(1,len(vecs1)):
        val += (vecs1[i] - vecs2[i])**2
    if (val > MMdistA):
        return False
    idx2 = fMap[fUpper2]
    vecs1 = mus1[idx2]
    vecs2 = mus2[idx2]
    val = (vecs1[0] - vecs2[0])**2
    for i in range(1,len(vecs1)):
        val += (vecs1[i] - vecs2[i])**2
    return (val < MMdistA)

def return_nearest_isco(totmass, freqs):
    """
    Given a value for total mass and a list of discrete frequencies, this will
    return the frequency in the list closest to the ISCO.

    Parameters
    ----------
    totmass : numpy.array
        The total mass of the input systems
    freqs : list of floats
        A list of frequencies

    Returns
    -------
    numpy.array
        The frequencies closest to the ISCO frequency corresponding to each
        value of totmass. 
    """
    # FIXME: I'm not entirely sure how this works! Documentation may be wrong.
    fISCO = (1/6.)**(3./2.) / (LAL_PI * totmass * LAL_MTSUN_SI)
    refEv = numpy.zeros(len(fISCO),dtype=float)
    for i in range(len(freqs)):
        if (i == 0):
            logicArr = fISCO < ((freqs[0] + freqs[1])/2.)
        elif (i == (len(freqs)-1)):
            logicArr = fISCO > ((freqs[-2] + freqs[-1])/2.)
        else:
            logicArrA = fISCO > ((freqs[i-1] + freqs[i])/2.)
            logicArrB = fISCO < ((freqs[i] + freqs[i+1])/2.)
            logicArr = numpy.logical_and(logicArrA,logicArrB)
        if logicArr.any():
            refEv[logicArr] = freqs[i]
    return refEv

