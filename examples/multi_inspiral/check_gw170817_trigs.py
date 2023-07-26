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
gw170817_time = 1187008882.43
status = 0
with h5py.File('GW170817_test_output.hdf', 'r') as f:
    snrs = [
        f['network/end_time_gc'][:],
        f['network/coherent_snr'][:],
        f['network/reweighted_snr'][:],
        f['network/slide_id'][:]]
# search for compatible trigs
mask = (
    (abs(gw170817_time - snrs[0]) < 0.1)
    & (snrs[1] > 25)
    & (snrs[2] > 25)
    & (snrs[3] == 0)
    )
n = mask.sum()
if n > 0:
    result = 'PASS'
    status = 0
else:
    result = 'FAIL'
    status = 1
logging.info(
    '%s: GW170817 found with coherent SNR = %.2f; reweighted SNR %.2f', result,
    snrs[1][mask], snrs[2][mask])
sys.exit(status)
