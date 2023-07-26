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
import re
import numpy
import pycbc.libutils
from lal import MTSUN_SI, PI, CreateREAL8Vector

lalsimulation = pycbc.libutils.import_optional('lalsimulation')

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
    # For some as-of-yet unknown reason, the tidal terms are not giving correct
    # match estimates when enabled. So, for now, this order is commented out.
    #if order == 'tidalTesting':
    #    mapping['Lambda10'] = 8
    #    mapping['Lambda12'] = 9
    #    return mapping
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

def get_chirp_params(mass1, mass2, spin1z, spin2z, f0, order,
                     quadparam1=None, quadparam2=None, lambda1=None,
                     lambda2=None):
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
        if quadparam1 is not None:
            quadparam1 = numpy.array([quadparam1])
        if quadparam2 is not None:
            quadparam2 = numpy.array([quadparam2])
        if lambda1 is not None:
            lambda1 = numpy.array([lambda1])
        if lambda2 is not None:
            lambda2 = numpy.array([lambda2])
        num_points = 1

    if quadparam1 is None:
        quadparam1 = numpy.ones(len(mass1), dtype=float)
    if quadparam2 is None:
        quadparam2 = numpy.ones(len(mass1), dtype=float)
    if lambda1 is None:
        lambda1 = numpy.zeros(len(mass1), dtype=float)
    if lambda2 is None:
        lambda2 = numpy.zeros(len(mass1), dtype=float)

    mass1_v = CreateREAL8Vector(len(mass1))
    mass1_v.data[:] = mass1[:]
    mass2_v = CreateREAL8Vector(len(mass1))
    mass2_v.data[:] = mass2[:]
    spin1z_v = CreateREAL8Vector(len(mass1))
    spin1z_v.data[:] = spin1z[:]
    spin2z_v = CreateREAL8Vector(len(mass1))
    spin2z_v.data[:] = spin2z[:]
    lambda1_v = CreateREAL8Vector(len(mass1))
    lambda1_v.data[:] = lambda1[:]
    lambda2_v = CreateREAL8Vector(len(mass1))
    lambda2_v.data[:] = lambda2[:]
    dquadparam1_v = CreateREAL8Vector(len(mass1))
    dquadparam1_v.data[:] = quadparam1[:] - 1.
    dquadparam2_v = CreateREAL8Vector(len(mass1))
    dquadparam2_v.data[:] = quadparam2[:] - 1.

    phasing_arr = lalsimulation.SimInspiralTaylorF2AlignedPhasingArray\
        (mass1_v, mass2_v, spin1z_v, spin2z_v, lambda1_v, lambda2_v,
         dquadparam1_v, dquadparam2_v)

    vec_len = lalsimulation.PN_PHASING_SERIES_MAX_ORDER + 1;
    phasing_vs = numpy.zeros([num_points, vec_len])
    phasing_vlogvs = numpy.zeros([num_points, vec_len])
    phasing_vlogvsqs = numpy.zeros([num_points, vec_len])

    lng = len(mass1)
    jmp = lng * vec_len
    for idx in range(vec_len):
        phasing_vs[:,idx] = phasing_arr.data[lng*idx : lng*(idx+1)]
        phasing_vlogvs[:,idx] = \
            phasing_arr.data[jmp + lng*idx : jmp + lng*(idx+1)]
        phasing_vlogvsqs[:,idx] = \
            phasing_arr.data[2*jmp + lng*idx : 2*jmp + lng*(idx+1)]

    pim = PI * (mass1 + mass2)*MTSUN_SI
    pmf = pim * f0
    pmf13 = pmf**(1./3.)
    logpim13 = numpy.log((pim)**(1./3.))

    mapping = generate_inverse_mapping(order)
    lambdas = []
    lambda_str = '^Lambda([0-9]+)'
    loglambda_str = '^LogLambda([0-9]+)'
    logloglambda_str = '^LogLogLambda([0-9]+)'
    for idx in range(len(mapping.keys())):
        # RE magic engage!
        rematch = re.match(lambda_str, mapping[idx])
        if rematch:
            pn_order = int(rematch.groups()[0])
            term = phasing_vs[:,pn_order]
            term = term + logpim13 * phasing_vlogvs[:,pn_order]
            lambdas.append(term * pmf13**(-5+pn_order))
            continue
        rematch = re.match(loglambda_str, mapping[idx])
        if rematch:
            pn_order = int(rematch.groups()[0])
            lambdas.append((phasing_vlogvs[:,pn_order]) * pmf13**(-5+pn_order))
            continue
        rematch = re.match(logloglambda_str, mapping[idx])
        if rematch:
            raise ValueError("LOGLOG terms are not implemented")
            #pn_order = int(rematch.groups()[0])
            #lambdas.append(phasing_vlogvsqs[:,pn_order] * pmf13**(-5+pn_order))
            #continue
        err_msg = "Failed to parse " +  mapping[idx]
        raise ValueError(err_msg)

    if sngl_inp:
        return [l[0] for l in lambdas]
    else:
        return lambdas

get_chirp_params.__doc__ = \
    get_chirp_params.__doc__.format(pycbcValidOrdersHelpDescriptions)
