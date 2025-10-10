#!/usr/bin/env python

# Read a pycbc_multi_inspiral HDF5 trigger file and check that it contains
# triggers compatible with GW170817
# 2022 Andrew Williamson, Tito Dal Canton

import sys
import logging
import h5py


with h5py.File('GW170817_test_output.hdf', 'r') as f:
    end_time = f['network/end_time_gc'][:]
    coh_snr = f['network/coherent_snr'][:]
    rw_snr = f['network/reweighted_snr'][:]
    slide_id = f['network/slide_id'][:]
    present_detectors = ''.join(
        d for d in 'HLV' if f'{d}1' in f.keys()
    )

# basic sanity checks
assert (end_time > 0).all()
assert (coh_snr > 0).all()
assert (rw_snr > 0).all()

# search for trigs compatible with GW170817
expected_time = 1187008882.43
required_net_snr = {'H': 17, 'HL': 28, 'HLV': 28}
mask = (
    (abs(expected_time - end_time) < 0.1)
    & (coh_snr > required_net_snr[present_detectors])
    & (rw_snr > required_net_snr[present_detectors])
    & (slide_id == 0)
)
n = mask.sum()
if n > 0:
    print(
        f'PASS: GW170817 found with max coherent SNR {max(coh_snr[mask]):.2f}, '
        f'max reweighted SNR {max(rw_snr[mask]):.2f}'
    )
    status = 0
else:
    print('FAIL: GW170817 not found')
    status = 1

sys.exit(status)
