#!/usr/bin/env python

import os, sys
from pycbc.io import record
from pycbc.inject import InjectionSet

if os.path.exists('./test_inj1.hdf'):
    raise OSError("output-file 1 already exists")
    
if os.path.exists('./test_inj2.hdf'):
    raise OSError("output-file 2 already exists")

    
# injection 1    
static_params = { 'f_lower': 18.0, 'f_ref': 18.0, 'approximant': 'SEOBNRv4', 
               'taper': 'start', 'ra': 45.0, 'dec': 45.0, 'inclination': 0.0, 
               'coa_phase': 0.0, 'polarization': 0.0}

samples = record.FieldArray(2, dtype=[('mass1', float), ('mass2', float), 
                                      ('spin1z', float), ('spin2z', float), 
                                      ('tc', float), ('distance', float)])

    
# The following 'magic numbers' are intended to match the highest
# mass injection in the template bank
samples['mass1'] = [290.929321]
samples['mass2'] = [3.6755455]
samples['spin1z'] = [0.9934847]
samples['spin2z'] = [0.92713535]
samples['tc'] = [1272790100.1]
samples['distance'] = [603.0]
     
InjectionSet.write('test_inj1.hdf', samples, static_args = static_params,
                   injtype = 'cbc', cmd=" ".join(sys.argv))



#injection 2    
static_params = { 'f_lower': 18.0, 'f_ref': 18.0, 'approximant': 'SpinTaylorT4', 
               'taper': 'start', 'ra': 45.0, 'dec': 45.0, 'inclination': 0.0, 
               'coa_phase': 0.0, 'polarization': 0.0}

samples = record.FieldArray(2, dtype=[('mass1', float), ('mass2', float), 
                                      ('spin1z', float), ('spin2z', float), 
                                      ('tc', float), ('distance', float)])

    
# The following 'magic numbers' are intended to match the lowest
# mass injection in the template bank
samples['mass1'] = [1.1331687]
samples['mass2'] = [1.010624]
samples['spin1z'] = [0.029544285]
samples['spin2z'] = [0.020993788]
samples['tc'] = [1272790260.1]
samples['distance'] = [72.0]
     
InjectionSet.write('test_inj2.hdf', samples, static_args = static_params,
                   injtype = 'cbc', cmd=" ".join(sys.argv))