#!/usr/bin/env python

import sys
import numpy as np
from pycbc.io import FieldArray
from pycbc.inject import InjectionSet


dtype = [('mass1', float), ('mass2', float),
         ('spin1z', float), ('spin2z', float),
         ('tc', float), ('distance', float),
         ('ra', float), ('dec', float),
         ('approximant', 'S32')]

static_params = {'f_lower': 17.,
                 'f_ref': 17.,
                 'taper': 'start',
                 'inclination': 0.,
                 'coa_phase': 0.,
                 'polarization': 0.}

samples = FieldArray(3, dtype=dtype)

# masses and spins are intended to match the highest
# and lowest mass templates in the template bank
# Last injection is designed to be found as an EM-bright single
samples['mass1'] = [290.929321, 1.1331687, 2.2756491]
samples['mass2'] = [3.6755455, 1.010624, 1.1077247]
samples['spin1z'] = [0.9934847, 0.029544285, -0.59105825]
samples['spin2z'] = [0.92713535, 0.020993788, 0.047548451]

# distance and sky locations for coincs to have network SNRs ~15
# and for single to pass SNR cuts
samples['tc'] = [1272790100.1, 1272790260.1, 1272790490.2]
samples['distance'] = [178., 79., 47.]
samples['ra'] = [np.deg2rad(45), np.deg2rad(10), np.deg2rad(10)]
samples['dec'] = [np.deg2rad(45), np.deg2rad(-45), np.deg2rad(-45)]

samples['approximant'] = ['SEOBNRv4_opt', 'SpinTaylorT4', 'SpinTaylorT4']

InjectionSet.write('injections.hdf', samples, static_args=static_params,
                   injtype='cbc', cmd=" ".join(sys.argv))
