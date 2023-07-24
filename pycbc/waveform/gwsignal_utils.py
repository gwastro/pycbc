'''
Code to covert pycbc waveform gen params to gwsignal params
'''
# import lalsimulation as lalsim
from lalsimulation.gwsignal.core.parameter_conventions import (default_dict,
        common_units_dictionary)

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
## There is additional f_ref/min (other than f22_ref/min) in gws parameters. Don't know
## what that is?

def to_gwsignal_dict(par):
    '''convert param dict to gws dict
    '''
    params = par.copy()
    for key in par:
        # if par[key]:
        knew = pycbc_to_gws.get(key, key)
        params[knew] = params.pop(key)
        params[knew] = (params[knew]*gws_units.get(knew)) if params[knew] \
                else params[knew]

    _ = params.setdefault('condition', 1)

    return params
