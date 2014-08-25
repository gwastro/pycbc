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
from lal import MTSUN_SI, GAMMA

# PLEASE ENSURE THESE ARE KEPT UP TO DATE WITH THE REST OF THIS FILE
pycbcValidTmpltbankOrders = ['zeroPN','onePN','onePointFivePN','twoPN',\
      'twoPointFivePN','threePN','threePointFivePN','taylorF4_45PN']

pycbcValidOrdersHelpDescriptions="""
* zeroPN: Will only include the dominant term (proportional to chirp mass)
* onePN: Will only the leading orbit term and first correction at 1PN
* onePointFivePN: Will include orbit and spin terms to 1.5PN.
* twoPN: Will include orbit and spin terms to 2PN.
* twoPointFivePN: Will include orbit and spin terms to 2.5PN.
* threePN: Will include orbit terms to 3PN and spin terms to 2.5PN.
* threePointFivePN: Include orbit terms to 3.5PN and spin terms to 2.5PN
* taylorF4_45PN: Will use the R2F4 metric to 4.5PN. This includes partial
spin terms from 3 to 4.5PN and partial orbit terms from 4 to 4.5PN.
"""


def generate_mapping(order):
    """
    This function will take an order string and return a mapping between
    components in the metric and the various Lambda components. This must be
    used (and consistently used) when generating the metric *and* when
    transforming to/from the xi_i coordinates to the lambda_i coordinates.

    NOTE: This is not a great way of doing this. It would be nice to clean
    this up. Hence pulling this function out. The valid PN orders are: %s

    Parameters
    ----------
    order : string
        A string containing a PN order. Valid values are given above.
    Returns
    --------
    mapping: dictionary
        A mapping between the active Lambda terms and index in the metric
    """ %(pycbcValidOrdersHelpDescriptions)
    mapping = {}
    if order == 'taylorF4_45PN':
        mapping['Lambda0'] = 0
        mapping['Lambda2'] = 1
        mapping['Lambda3'] = 2
        mapping['Lambda4'] = 3
        mapping['LogLambda5'] = 4
        mapping['Lambda6'] = 5
        mapping['LogLambda6'] = 6
        mapping['Lambda7'] = 7
        mapping['LogLambda8'] = 8
        mapping['LogLogLambda8'] = 9
        mapping['Lambda9'] = 10
        mapping['LogLambda9'] = 11
        return mapping
    else:
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

def generate_inverse_mapping(order):
    """
    This function will generate the opposite of generate mapping. So where
    generate_mapping gives:
    dict[key] = item
    This will give
    dict[item] = key
    Valid PN orders are: %s
    
    Parameters
    ----------
    order : string
        A string containing a PN order. Valid values are given above.
    Returns
    --------
    mapping: dictionary
        An inverse mapping between the active Lambda terms and index in the
        metric
    """ %(pycbcValidOrdersHelpDescriptions)
    mapping = generate_mapping(order)
    inv_mapping = {}
    for key,value in mapping.items():
        inv_mapping[value] = key

    return inv_mapping

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

def get_chirp_params(totmass, eta, beta, sigma, gamma, chis, f0, order):
    """
    Take a set of masses and spins and convert to the various lambda
    coordinates that describe the orbital phase. Accepted PN orders are: %s
 
    Parameters
    ----------
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
        spins to the lambda_i coordinate system. Valid orders given above.

    Returns
    --------
    lambdas : list of floats or numpy.arrays
        The lambda coordinates for the input system(s)
    """ %(pycbcValidOrdersHelpDescriptions)

    # Convert mass to seconds
    totmass = totmass * MTSUN_SI
    pi = numpy.pi
    mapping = generate_inverse_mapping(order)
    if order[0:8] == 'taylorF4':
        r2f4Terms = True
    else:
        r2f4Terms = False
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
            if not r2f4Terms:
                lambdas.append(lambda6)
                continue
            lambda6spin = 502.6548245743669 * beta + 88.45238095238095 * sigma\
                          + (110. * eta * sigma) - 20. * beta * beta
            lambda6spin = lambda6spin * 3./(128.*eta) \
                          * (pi * totmass * f0)**(1/3.)
            lambda6 += lambda6spin
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
            if not r2f4Terms:
                lambdas.append(lambda7)
                continue
            lambda7spin = -510.0603994773098*beta \
                          - 368.01525846326734*beta*eta \
                          + 1944.363555525455*chis*eta \
                          - 502.6548245743669*sigma \
                          + 40.*beta*sigma + 241.47615535889872*beta*eta*eta \
                          + 2961.654024441635*chis*eta*eta \
                          + 676.0619469026549*chis*eta*eta*eta
            lambda7spin = lambda7spin * 3./(128.*eta) * \
                          (pi * totmass * f0)**(2/3.)
            lambdas.append(lambda7 + lambda7spin)
        elif mapping[idx] == 'Lambda8':
            # This term should not be present as it is equivalent to a
            # constant time offset
            lambda8orbit = 342.6916926002232 + 2869.024558661873*eta - \
                3773.659169914512*eta*eta + 172.0609149438239*eta*eta*eta - \
                24.284336419753085*eta*eta*eta*eta
            lambda8spin = 936.7471880419974*beta \
                          - 311.03929625364435*beta*eta \
                          - 2455.4171568883194*chis*eta \
                          + 195.39588893571195*beta*chis*eta \
                          + 48.491201779065534*sigma \
                          + 101.92901234567901*eta*sigma \
                          - 58.81315259633844*beta*beta \
                          + 8.918387413962636*eta*beta*beta \
                          - 686.5167663065837*chis*eta*eta \
                          + 54.631268436578175*beta*chis*eta*eta \
                          + 71.69753086419753*sigma*eta*eta \
                          - 4.444444444444445*sigma*sigma
            lambda8 = (lambda8 + lambda8spin) * 3./(128.*eta) \
                      * (pi * totmass * f0)**(3./3.)
            lambdas.append(lambda8)
        elif mapping[idx] == 'LogLambda8':
            loglambda8orbit = -1028.0750778006693 - 8607.073675985623*eta \
                         + 11320.977509743536*eta*eta \
                         - 516.1827448314717*eta*eta*eta \
                         + 72.85300925925927*eta*eta*eta*eta
            loglambda8spin = -2810.241564125992*beta \
                             + 933.117888760933*beta*eta \
                             + 7366.251470664957*chis*eta \
                             - 586.1876668071359*beta*chis*eta \
                             - 145.4736053371966*sigma \
                             - 305.78703703703707*eta*sigma \
                             + 176.4394577890153*beta*beta \
                             - 26.755162241887906*eta*beta*beta \
                             + 2059.5502989197507*chis*eta*eta \
                             - 163.89380530973452*beta*chis*eta*eta \
                             - 215.0925925925926*sigma*eta*eta \
                             + 13.333333333333334*sigma*sigma
            loglambda8 = (loglambda8orbit + loglambda8spin) * 3./(128.*eta) \
                         * (pi * totmass * f0)**(3./3.)  
            lambdas.append(loglambda8)
        elif mapping[idx] == 'LogLogLambda8':
            logloglambda8 = 480.7316704459562 + 597.8412698412699*eta
            logloglambda8 = logloglambda8 * 3./(128.*eta) \
                            * (pi * totmass * f0)**(3./3.)
            lambdas.append(logloglambda8)
        elif mapping[idx] == 'Lambda9':
            lambda9orbit = 20021.24571514093 - 42141.161261993766*eta \
                           - 4047.211701119762*eta*eta \
                           - 2683.4848475303043*eta*eta*eta
            lambda9spin = 2105.9471080635244*beta \
                          + 3909.271818583914*beta*eta \
                          - 2398.354686411564*chis*eta \
                          + 1278.4225104920606*sigma \
                          - 198.6688790560472*beta*sigma \
                          + 589.0486225480862*eta*sigma \
                          - 62.43362831858406*beta*eta*sigma \
                          + 439.6407501053519*chis*eta*sigma \
                          - 376.99111843077515*beta*beta \
                          + 10.*beta*beta*beta \
                          - 202.1451909795383*beta*eta*eta \
                          - 5711.929102446965*chis*eta*eta \
                          + 122.9203539823009*chis*sigma*eta*eta \
                          - 493.00738145963066*beta*eta*eta*eta \
                          - 4955.659484448894*chis*eta*eta*eta \
                          - 991.4721607669617*chis*eta*eta*eta*eta
            lambda9 = (lambda9orbit + lambda9spin) * 3./(128.*eta) \
                      * (pi * totmass * f0)**(4./3.)
            lambdas.append(lambda9)
        elif mapping[idx] == 'LogLambda9':
            loglambda9orbit = -4097.833617482457
            loglambda9spin = 326.0952380952381*beta
            loglambda9 = (loglambda9orbit + loglambda9spin) * 3./(128.*eta) \
                         * (pi * totmass * f0)**(4./3.) 
            lambdas.append(loglambda9)
        else:
            errMsg = "Do not understand term %s." %(mapping[idx])
            raise ValueError(errMsg)
                 
    return lambdas
