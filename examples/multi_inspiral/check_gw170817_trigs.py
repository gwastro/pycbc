#!/usr/bin/env python

# Read a pycbc_multi_inspiral HDF5 trigger file and check that it contains
# triggers compatible with GW170817
# 2022 Andrew Williamson, Tito Dal Canton

import sys
import logging
import h5py


expected_time = 1187008882.43

with h5py.File('GW170817_test_output.hdf', 'r') as f:
    end_time = f['network/end_time_gc'][:]
    coh_snr = f['network/coherent_snr'][:]
    rw_snr = f['network/reweighted_snr'][:]
    slide_id = f['network/slide_id'][:]

# search for compatible trigs
mask = (
    (abs(expected_time - end_time) < 0.1)
    & (coh_snr > 25)
    & (rw_snr > 25)
    & (slide_id == 0)
)
n = mask.sum()
if n > 0:
    print(
        'PASS: GW170817 found with coherent SNR = %.2f; reweighted SNR %.2f',
        coh_snr[mask],
        rw_snr[mask]
    )
    status = 0
else:
    print('FAIL: GW170817 not found')
    status = 1

sys.exit(status)
