#!/usr/bin/env python
# -*- coding: utf-8 -*-
import astropy.units as u

pycbc_to_gws = {
        'delta_t': 'deltaT',
        'delta_f': 'deltaF',
        'f_lower': 'f22_start',
        'f_ref': 'f22_ref',
        'coa_phase': 'phi_ref',
        'long_asc_nodes': 'longAscNodes',
        'mean_per_ano': 'meanPerAno'
        }

gws_units = {'mass1': u.solMass,
             'mass2': u.solMass,
             'spin1x': u.dimensionless_unscaled,
             'spin1y': u.dimensionless_unscaled,
             'spin1z': u.dimensionless_unscaled,
             'spin2x': u.dimensionless_unscaled,
             'spin2y': u.dimensionless_unscaled,
             'spin2z': u.dimensionless_unscaled,
             'deltaT': u.s,
             'deltaF': u.Hz,
             'f22_stat': u.Hz,
             'f22_ref': u.Hz,
             'phi_ref': u.rad,
             'distance': u.Mpc,
             'inclination': u.rad,
             'eccentricity': u.dimensionless_unscaled,
             'longAscNodes': u.rad,
             'meanPerAno': u.rad}

def to_gwsignal_dict(p):
    params = {}
    for k in p:
        knew = pycbc_to_gws.get(k)
        if not knew:
            knew = k

        if knew in gws_units and p[k]:
            params[knew] = p[k] * gws_units[knew]

    return params
