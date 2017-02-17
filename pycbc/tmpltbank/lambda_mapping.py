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
import re
import numpy
from lal import MTSUN_SI, GAMMA, PI, CreateDict
import lalsimulation
from pycbc import pnutils

# PLEASE ENSURE THESE ARE KEPT UP TO DATE WITH THE REST OF THIS FILE
pycbcValidTmpltbankOrders = ['zeroPN','onePN','onePointFivePN','twoPN',\
      'twoPointFivePN','threePN','threePointFivePN']

pycbcValidOrdersHelpDescriptions="""
     * zeroPN: Will only include the dominant term (proportional to chirp mass)
     * onePN: Will only the leading orbit term and first correction at 1PN
     * onePointFivePN: Will include orbit and spin terms to 1.5PN.
     * twoPN: Will include orbit and spin terms to 2PN.
     * twoPointFivePN: Will include orbit and spin terms to 2.5PN.
     * threePN: Will include orbit terms to 3PN and spin terms to 2.5PN.
     * threePointFivePN: Include orbit terms to 3.5PN and spin terms to 2.5PN
"""


def generate_mapping(order):
    """
    This function will take an order string and return a mapping between
    components in the metric and the various Lambda components. This must be
    used (and consistently used) when generating the metric *and* when
    transforming to/from the xi_i coordinates to the lambda_i coordinates.

    NOTE: This is not a great way of doing this. It would be nice to clean
    this up. Hence pulling this function out. The valid PN orders are
    {}

    Parameters
    ----------
    order : string
        A string containing a PN order. Valid values are given above.

    Returns
    --------
    mapping : dictionary
        A mapping between the active Lambda terms and index in the metric
    """
    mapping = {}
    mapping['Lambda0'] = 0
    if order == 'zeroPN':
        return mapping
    mapping['Lambda2'] = 1
    if order == 'onePN':
        return mapping
    mapping['Lambda3'] = 2
    if order == 'onePointFivePN':
        return mapping
    mapping['Lambda4'] = 3
    if order == 'twoPN':
        return mapping
    mapping['LogLambda5'] = 4
    if order == 'twoPointFivePN':
        return mapping
    mapping['Lambda6'] = 5
    mapping['LogLambda6'] = 6
    if order == 'threePN':
        return mapping
    mapping['Lambda7'] = 7
    if order == 'threePointFivePN':
        return mapping
    raise ValueError("Order %s is not understood." %(order))

# Override doc so the PN orders are added automatically to online docs
generate_mapping.__doc__ = \
    generate_mapping.__doc__.format(pycbcValidOrdersHelpDescriptions)

def generate_inverse_mapping(order):
    """Genereate a lambda entry -> PN order map.

    This function will generate the opposite of generate mapping. So where
    generate_mapping gives dict[key] = item this will give
    dict[item] = key. Valid PN orders are:
    {}
    
    Parameters
    ----------
    order : string
        A string containing a PN order. Valid values are given above.

    Returns
    --------
    mapping : dictionary
        An inverse mapping between the active Lambda terms and index in the
        metric
    """
    mapping = generate_mapping(order)
    inv_mapping = {}
    for key,value in mapping.items():
        inv_mapping[value] = key

    return inv_mapping

generate_inverse_mapping.__doc__ = \
    generate_inverse_mapping.__doc__.format(pycbcValidOrdersHelpDescriptions)

def get_ethinca_orders():
    """
    Returns the dictionary mapping TaylorF2 PN order names to twice-PN 
    orders (powers of v/c)
    """
    ethinca_orders = {"zeroPN"           : 0,
                      "onePN"            : 2,
                      "onePointFivePN"   : 3,
                      "twoPN"            : 4,
                      "twoPointFivePN"   : 5,
                      "threePN"          : 6,
                      "threePointFivePN" : 7
                     }
    return ethinca_orders

def ethinca_order_from_string(order):
    """
    Returns the integer giving twice the post-Newtonian order 
    used by the ethinca calculation. Currently valid only for TaylorF2 metric

    Parameters
    ----------
    order : string
    
    Returns
    -------
    int
    """
    if order in get_ethinca_orders().keys():
      return get_ethinca_orders()[order]
    else: raise ValueError("Order "+str(order)+" is not valid for ethinca"
                           "calculation! Valid orders: "+
                           str(get_ethinca_orders().keys()))

def get_chirp_params_new(mass1, mass2, spin1z, spin2z, f0, order):
    """
    Take a set of masses and spins and convert to the various lambda
    coordinates that describe the orbital phase. Accepted PN orders are:
    {}
 
    Parameters
    ----------
    mass1 : float or array
        Mass1 of input(s).
    mass2 : float or array
        Mass2 of input(s).
    spin1z : float or array
        Parallel spin component(s) of body 1.
    spin2z : float or array
        Parallel spin component(s) of body 2.
    f0 : float
        This is an arbitrary scaling factor introduced to avoid the potential
        for numerical overflow when calculating this. Generally the default
        value (70) is safe here. **IMPORTANT, if you want to calculate the
        ethinca metric components later this MUST be set equal to f_low.**
        This value must also be used consistently (ie. don't change its value
        when calling different functions!).
    order : string
        The Post-Newtonian order that is used to translate from masses and
        spins to the lambda_i coordinate system. Valid orders given above.

    Returns
    --------
    lambdas : list of floats or numpy.arrays
        The lambda coordinates for the input system(s)
    """

    # Determine whether array or single value input
    sngl_inp = False
    try:
        num_points = len(mass1)
    except TypeError:
        sngl_inp = True
        # If you care about speed, you aren't calling this function one entry
        # at a time.
        mass1 = numpy.array([mass1])
        mass2 = numpy.array([mass2])
        spin1z = numpy.array([spin1z])
        spin2z = numpy.array([spin2z])
        num_points = 1
    lal_pars = CreateDict()
    phasing_vs = numpy.zeros([num_points, 13])
    phasing_vlogvs = numpy.zeros([num_points, 13])
    phasing_vlogvsqs = numpy.zeros([num_points, 13])
    for i in xrange(num_points):
        phasing = lalsimulation.SimInspiralTaylorF2AlignedPhasing(
                            mass1[i], mass2[i], spin1z[i], spin2z[i], lal_pars)
        phasing_vs[i] = phasing.v
        phasing_vlogvs[i] = phasing.vlogv
        phasing_vlogvsqs[i] = phasing.vlogvsq

    pmf = PI * (mass1 + mass2)*MTSUN_SI * f0
    pmf13 = pmf**(1./3.)

    mapping = generate_inverse_mapping(order)
    lambdas = []
    lambda_str = '^Lambda([0-9]+)'
    loglambda_str = '^LogLambda([0-9]+)'
    logloglambda_str = '^LogLogLambda([0-9]+'
    for idx in xrange(len(mapping.keys())):
        # RE magic engage!
        rematch = re.match(lambda_str, mapping[idx])
        if rematch:
            pn_order = int(rematch.groups()[0])
            lambdas.append(phasing_vs[:,pn_order] * pmf13**(-5+pn_order))
            continue
        rematch = re.match(loglambda_str, mapping[idx])
        if rematch:
            pn_order = int(rematch.groups()[0])
            lambdas.append(phasing_vlogvs[:,pn_order] * pmf13**(-5+pn_order))
            continue
        rematch = re.match(logloglambda_str, mapping[idx])
        if rematch:
            pn_order = int(rematch.groups()[0])
            lambdas.append(phasing_vlogvsqs[:,pn_order] * pmf13**(-5+pn_order))
            continue
        err_msg = "Failed to parse " +  mapping[idx]
        raise ValueError(err_msg)

    if sngl_inp:
        return [l[0] for l in lambdas]
    else:
        return lambdas

get_chirp_params_new.__doc__ = \
    get_chirp_params_new.__doc__.format(pycbcValidOrdersHelpDescriptions)

def get_chirp_params_old(mass1, mass2, spin1z, spin2z, f0, order):
    """
    Take a set of masses and spins and convert to the various lambda
    coordinates that describe the orbital phase. Accepted PN orders are:
    {}
 
    Parameters
    ----------
    mass1 : float or array
        Mass1 of input(s).
    mass2 : float or array
        Mass2 of input(s).
    spin1z : float or array
        Parallel spin component(s) of body 1.
    spin2z : float or array
        Parallel spin component(s) of body 2.
    f0 : float
        This is an arbitrary scaling factor introduced to avoid the potential
        for numerical overflow when calculating this. Generally the default
        value (70) is safe here. **IMPORTANT, if you want to calculate the
        ethinca metric components later this MUST be set equal to f_low.**
        This value must also be used consistently (ie. don't change its value
        when calling different functions!).
    order : string
        The Post-Newtonian order that is used to translate from masses and
        spins to the lambda_i coordinate system. Valid orders given above.

    Returns
    --------
    lambdas : list of floats or numpy.arrays
        The lambda coordinates for the input system(s)
    """

    totmass, eta = pnutils.mass1_mass2_to_mtotal_eta(mass1, mass2)
    beta, sigma, gamma = pnutils.get_beta_sigma_from_aligned_spins(\
               eta, spin1z, spin2z)

    # Convert mass to seconds
    totmass = totmass * MTSUN_SI
    pi = numpy.pi
    mapping = generate_inverse_mapping(order)
    lambdas = []

    for idx in xrange(len(mapping.keys())):
        if mapping[idx] == 'Lambda0':
            lambda0 = 3. / (128. * eta * (pi * totmass * f0)**(5./3.))
            lambdas.append(lambda0)
        elif mapping[idx] == 'Lambda2':
            lambda2 = 5. / (96. * pi * eta * totmass * f0) \
                      * (743./336. + 11./4. * eta)
            lambdas.append(lambda2)
        elif mapping[idx] == 'Lambda3':
            lambda3 = (-3. * pi**(1./3.))/(8. * eta * (totmass*f0)**(2./3.)) \
                      * (1. - beta/ (4. * pi))
            lambdas.append(lambda3)
        elif mapping[idx] == 'Lambda4':
            lambda4 = 15. / (64. * eta * (pi * totmass * f0)**(1./3.)) * \
                  (3058673./1016064. + 5429./1008. * eta + 617./144. * \
                   eta**2 - sigma)
            lambdas.append(lambda4)
        # No Lambda5 term is present as that would be equivalent to a constant
        # phase offset, and thus completely degenerate with the initial orbital
        # phase.
        elif mapping[idx] == 'LogLambda5':
            loglambda5 = 3. * (38645.*pi/756. - 65.*pi*eta/9. - gamma)
            loglambda5 = loglambda5 * (3./(128.*eta))
            lambdas.append(loglambda5)
        elif mapping[idx] == 'Lambda6':
            lambda6 = 11583231236531./4694215680. - (640.*pi*pi)/3.\
                      - (6848.*GAMMA)/21.
            lambda6 -= (6848./21.) * numpy.log(4 * (pi*totmass*f0)**(1./3.))
            lambda6 += (-15737765635/3048192. + 2255.*pi*pi/12.)*eta
            lambda6 += (76055.*eta*eta)/1728. - (127825.*eta*eta*eta)/1296.
            lambda6 = lambda6 * 3./(128.*eta) * (pi * totmass * f0)**(1/3.)
            lambdas.append(lambda6)
        elif mapping[idx] == 'LogLambda6':
            loglambda6 =  -( 6848./21) 
            loglambda6 = loglambda6 * 3./(128.*eta)\
                         * (pi * totmass * f0)**(1/3.)
            lambdas.append(loglambda6)
        elif mapping[idx] == 'Lambda7':
            lambda7 = (77096675.*pi)/254016. + (378515.*pi*eta)/1512. \
                      - (74045.*pi*eta*eta)/756.
            lambda7 = lambda7 * 3./(128.*eta) * (pi * totmass * f0)**(2/3.)
            lambdas.append(lambda7)
        else:
            err_msg = "Do not understand term {}.".format(mapping[idx])
            raise ValueError(err_msg)
                 
    return lambdas

get_chirp_params_old.__doc__ = \
    get_chirp_params_old.__doc__.format(pycbcValidOrdersHelpDescriptions)

get_chirp_params = get_chirp_params_old
