#!/usr/bin/env python
# Read a pycbc_multi_inspiral HDF5 trigger file and check that it contains
# triggers compatible with mock GW170817-like injections
# 2022 Andrew Williamson, Tito Dal Canton

import sys
import logging
import h5py
import numpy as np
from pycbc import init_logging

init_logging(True)
gw170817_time = 1187008882
end_times = (np.arange(3) - 1) * 300 + gw170817_time
pols = ['standard', 'left', 'right', 'left+right']
refs = {
    'standard': np.array([38.8, 18.4, 39.4]),
    'left': np.array([23.5, 17.0, 38.9]),
    'right': np.array([38.1, 17.1, 24.3]),
    'left+right': np.array([38.1, 17.1, 38.9])
    }
status = 0
for pol in pols:
    with h5py.File(pol + '.hdf', 'r') as f:
        snrs = [f['network/end_time_gc'][:], f['network/coherent_snr'][:]]
    # search for compatible trigs
    mask = np.logical_and(
        abs(end_times - snrs[0]) < 0.1,
        snrs[1] > 0.9 * refs[pol],
        snrs[1] < 1.1 * refs[pol]
        )
    n = mask.sum()
    result = 'PASS' if n == 3 else 'FAIL'
    if result == 'FAIL':
        status = 1
    logging.info('"%s" polarization: %s (%d/3 triggers found)', pol, result, n)
sys.exit(status)
