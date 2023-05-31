#!/usr/bin/env python

import sys
import numpy as np
import h5py

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

samples = np.array([None] * 3, dtype=dtype)

# masses and spins are intended to match the highest
# and lowest mass templates in the template bank
# Last injection is designed to be found as an EM-bright single
samples['mass1'] = [290.929321, 41.1331687, 2.2756491]
samples['mass2'] = [3.6755455, 31.010624, 1.1077247]
samples['spin1z'] = [0.9934847, 0.029544285, -0.59105825]
samples['spin2z'] = [0.92713535, 0.020993788, 0.047548451]

# Injections must be between 1186740069 and 1186743653
samples['tc'] = [1186741000.1294786, 1186742000.91204, 1186743000.1298376]
samples['distance'] = [178., 130., 47.]
samples['ra'] = [np.deg2rad(45), np.deg2rad(10), np.deg2rad(-10)]
samples['dec'] = [np.deg2rad(45), np.deg2rad(-45), np.deg2rad(45)]

samples['approximant'] = ['SEOBNRv4_opt', 'SEOBNRv4_opt', 'SpinTaylorT4']


with h5py.File('injections.hdf','w') as injout:
    for k, dt in dtype:
        injout[k] = samples[k]
    for p, v in static_params.items():
        injout.attrs[p] = v
    injout.attrs['static_args'] = list(static_params.keys())
    injout.attrs['cmd'] = " ".join(sys.argv)
    injout.attrs['injtype'] = "cbc"
