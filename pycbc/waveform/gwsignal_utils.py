'''
Code to covert pycbc waveform gen params to gwsignal params
'''
from astropy import units as u
from lalsimulation.gwsignal.core.parameter_conventions import (default_dict,
        common_units_dictionary, full_parameter_list)

gws_units = common_units_dictionary.copy()
gws_units.update({k: v.unit for k, v in default_dict.items()})


pycbc_to_gws = {
        'delta_t': 'deltaT',
        'delta_f': 'deltaF',
        'f_lower': 'f22_start',
        'f_ref': 'f22_ref',
        'coa_phase': 'phi_ref',
        'long_asc_nodes': 'longAscNodes',
        'mean_per_ano': 'meanPerAno',
        'f_final': 'f_max',
        }


def to_gwsignal_dict(par):
    '''convert param dict to gws dict
    '''
    params = {pycbc_to_gws.get(k, k): v for k, v in par.items() if (v is not
        None and pycbc_to_gws.get(k, k) in gws_units)}

    for key in params:
        params[key] *= gws_units.get(key, u.dimensionless_unscaled)

    _ = params.setdefault('condition', 0)

    return params
