#!/usr/bin/env python
# -*- coding: utf-8 -*-
import astropy.units as u
from lalsimulation import gwsignal
from lalsimulation.gwsignal.core.parameter_conventions import (mass_params, mass_params_, default_dict, common_units_dictionary)


mass_dimensionless = {mass_dimensionless:u.dimensionless_unscaled for mass_dimensionless in mass_params if mass_dimensionless not in mass_params_}
gws_units = {k: v.unit for k,v in default_dict.items()}

pycbc_to_gws = {
        'delta_t': 'deltaT',
        'delta_f': 'deltaF',
        'f_lower': 'f22_start',
        'f_ref': 'f22_ref',
        'coa_phase': 'phi_ref',
        'long_asc_nodes': 'longAscNodes',
        'mean_per_ano': 'meanPerAno'
        }


def to_gwsignal_dict(p):
    params = p.copy()
    for key in pycbc_to_gws:
        params[pycbc_to_gws.get(key)] = params.pop(key)

    for key in p:
        if key in pycbc_to_gws:
            params[pycbc_to_gws.get(key)] = params.pop(key)

        params[key] *= gws_units.get(key)

    return params
