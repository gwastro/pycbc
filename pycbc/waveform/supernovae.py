"""Generate core-collapse supernovae waveforms
"""

import numpy
import math
import h5py
from pycbc.types import TimeSeries, FrequencySeries, float64, complex128, zeros
from pycbc.waveform.waveform import get_obj_attrs


def get_td_corecollapse_bounce_signal(template=None, **kwargs):
    """generates CCSNe waveform
    """
    
    # eos_dict = dict({0:'SFHo', 1:'SFHx', 2:'LS180', 3:'HSIUF', 
    #                  4:'LS220', 5:'GSHenFSU2.1', 6:'GShenFSU1.7', 
    #                  7:'LS375', 8:'HSTMA', 9:'HSFSG', 10:'HSDD2', 
    #                  11:'BHBL', 12:'BHBLP', 13:'all_eos'})

    # eos_index = math.floor(numpy.float(kwargs['eos_index']))

    # check if a hdf file with principal components is provided as an arg:
    if 'pc_hdf_file' in kwargs:
        pc_file = h5py.File(kwargs['pc_hdf_file'], 'r')
        # eos_name = pc_file.get(eos_dict[eos_index])
        # eos_data = pc_file.get(eos_name)
        # principal_components = numpy.array(eos_data.get('principal_components'))
        principal_components = numpy.array(pc_file.get('principal_components'))

    if 'principal_components' in kwargs:
        principal_components = kwargs['principal_components']

    
    
    if 'coefficients_array' in kwargs:
        coefficients_array = kwargs['coefficients_array']
    else:
        coeffs_keys = [x for x in kwargs if x.startswith('coeff_')]
        coeffs_keys = numpy.sort(numpy.array(coeffs_keys))
        coefficients_array = numpy.array([kwargs[x] for x in coeffs_keys])

    if 'tc' in kwargs:
        print('tc is being passed here')
        t_bounce = kwargs['tc']
        
    no_of_pcs = int(kwargs['no_of_pcs'])
    
    coefficients_array = coefficients_array[:no_of_pcs]
    principal_components = principal_components[:no_of_pcs]

    pc_len = len(principal_components)
    assert len(coefficients_array) == pc_len

    distance = kwargs['distance']
    mpc_conversion = 3.08567758128e+22
    distance *=  mpc_conversion

    wf = numpy.dot(coefficients_array, principal_components) / distance
    
    delta_t = kwargs['delta_t']
    outhp = TimeSeries(wf, delta_t=delta_t)
    outhc = TimeSeries(numpy.zeros(len(wf)), delta_t=delta_t)
    
#    outhp.start_time = t_bounce - 0.5
#    print "tstart from the supernovae.py code: ", outhp.start_time
    # returning the same output for hp, hc as 2D waveforms don't have 
    # polarization info
    
    return outhp, outhc


# Approximant names ###########################################################
supernovae_td_approximants = {'CoreCollapseBounce' : get_td_corecollapse_bounce_signal}
