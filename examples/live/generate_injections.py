#!/usr/bin/env python

import sys
from pycbc.io import FieldArray
from pycbc.inject import InjectionSet


dtype = [('mass1', float), ('mass2', float),
         ('spin1z', float), ('spin2z', float),
         ('tc', float), ('distance', float),
         ('approximant', 'S32')]

static_params = {'f_lower': 18.0, 'f_ref': 18.0,
                 'taper': 'start', 'ra': 45.0, 'dec': 45.0,
                 'inclination': 0.0, 'coa_phase': 0.0, 'polarization': 0.0}

samples = FieldArray(2, dtype=dtype)

# The following 'magic numbers' are intended to match the highest
# and lowest mass templates in the template bank
samples['mass1'] = [290.929321, 1.1331687]
samples['mass2'] = [3.6755455, 1.010624]
samples['spin1z'] = [0.9934847, 0.029544285]
samples['spin2z'] = [0.92713535, 0.020993788]
samples['tc'] = [1272790100.1, 1272790260.1]
samples['distance'] = [301.5, 36.0]
samples['approximant'] = ['SEOBNRv4', 'SpinTaylorT4']

InjectionSet.write('injections.hdf', samples, static_args=static_params,
                   injtype='cbc', cmd=" ".join(sys.argv))
