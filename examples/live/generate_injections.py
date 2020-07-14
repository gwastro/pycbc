#!/usr/bin/env python

import os, sys
from pycbc.io import record
import h5py
from pycbc.inject import InjectionSet

if os.path.exists('./test_inj.hdf'):
    raise OSError("output-file already exists")

# 2 injections 

static_params = { 'f_lower': 18.0, 'f_ref': 18.0, \
               'taper': 'start', 'ra': 45.0, 'dec': 45.0, 'inclination': 0.0, \
               'coa_phase': 0.0, 'polarization': 0.0}

samples = record.FieldArray(2, dtype=[('approximant', h5py.string_dtype(encoding='utf-8')), \
                                      ('mass1', float), ('mass2', float), \
                                      ('spin1z', float), ('spin2z', float), \
                                      ('tc', float), ('distance', float)])

    
# The following 'magic numbers' are intended to match the highest and lowest
# mass injections in the template bank
samples['approximant'] = [u'SEOBNRv4', u'SpinTaylorT4']
samples['mass1'] = [290.929321, 1.1331687]
samples['mass2'] = [3.6755455, 1.010624]
samples['spin1z'] = [0.9934847, 0.029544285]
samples['spin2z'] = [0.92713535, 0.020993788]
samples['tc'] = [1272790100.1, 1272790260.1]
samples['distance'] = [292.0, 35.0]
     
InjectionSet.write('test_inj.hdf', samples, static_args = static_params,
                   injtype = 'cbc', cmd=" ".join(sys.argv))