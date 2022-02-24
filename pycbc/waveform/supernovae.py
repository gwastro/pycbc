"""Generate core-collapse supernovae waveform for core bounce and
subsequent postbounce oscillations.
"""

import numpy
import h5py
from pycbc.types import TimeSeries

_pc_dict = {}


def get_corecollapse_bounce(**kwargs):
    """ Generates core bounce and postbounce waveform by using principal
    component basis vectors from a .hdf file. The waveform parameters are the
    coefficients of the principal components and the distance. The number of
    principal components used can also be varied.
    """

    try:
        principal_components = _pc_dict['principal_components']
    except KeyError:
        with h5py.File(kwargs['principal_components_file'], 'r') as pc_file:
            principal_components = numpy.array(pc_file['principal_components'])
            _pc_dict['principal_components'] = principal_components

    if 'coefficients_array' in kwargs:
        coefficients_array = kwargs['coefficients_array']
    else:
        coeffs_keys = [x for x in kwargs if x.startswith('coeff_')]
        coeffs_keys = numpy.sort(numpy.array(coeffs_keys))
        coefficients_array = numpy.array([kwargs[x] for x in coeffs_keys])

    no_of_pcs = int(kwargs['no_of_pcs'])
    coefficients_array = coefficients_array[:no_of_pcs]
    principal_components = principal_components[:no_of_pcs]

    pc_len = len(principal_components)
    assert len(coefficients_array) == pc_len

    distance = kwargs['distance']
    mpc_conversion = 3.08567758128e+22
    distance *= mpc_conversion

    strain = numpy.dot(coefficients_array, principal_components) / distance
    delta_t = kwargs['delta_t']
    outhp = TimeSeries(strain, delta_t=delta_t)
    outhc = TimeSeries(numpy.zeros(len(strain)), delta_t=delta_t)
    return outhp, outhc


# Approximant names ###########################################################
supernovae_td_approximants = {'CoreCollapseBounce': get_corecollapse_bounce}
