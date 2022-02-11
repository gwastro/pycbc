from __future__ import division
import math
import numpy as np
import pycbc
import os.path
import scipy.interpolate
import logging
import sys

##############################################################################
# 2H 2-piecewise polytropic EOS, NS non-rotating equilibrium sequence        #
# File format is: grav mass (Msun)   baryonic mass (Msun)    compactness     #
#                                                                            #       THIS FUNCTION IS BEING CALLED IN TMPLTBANK, does it need to be moved back?
# Eventually, the function should call an NS sequence generator within LAL   #       or should I go and modify those files to reflect new location?
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
