#!/usr/bin/env python
# Read a pycbc_inspiral HDF5 trigger file and check that it contains triggers
# compatible with GW150914
# 2016 Tito Dal Canton

import sys
import h5py
import numpy as np

gw150914_time = 1126259462.4
gw150914_snr = {'H1': 19.71, 'L1': 13.28}
gw150914_chi2r = {'H1': 1.05, 'L1': 0.45}

f = h5py.File(sys.argv[1], 'r')
detector = f.keys()[0]
end_times = f[detector]['end_time'][:]
snrs = f[detector]['snr'][:]
chi2rs = f[detector]['chisq'][:] / (2 * f[detector]['chisq_dof'][:] - 2)

# search for trigs compatible with GW150914
mask = np.logical_and.reduce([abs(end_times - gw150914_time) < 0.1,
                              snrs > 0.8 * gw150914_snr[detector],
                              snrs < 1.2 * gw150914_snr[detector],
                              chi2rs > 0.8 * gw150914_chi2r[detector],
                              chi2rs < 1.2 * gw150914_chi2r[detector]])

if mask.any():
    print('Pass: %d GW150914-like triggers' % sum(mask))
    print('end_time snr reduced_chi2')
    for t, s, c in zip(end_times[mask], snrs[mask], chi2rs[mask]):
        print('%.3f %.3f %.3f' % (t, s, c))
    sys.exit(0)
else:
    print('Fail: no GW150914-like triggers')
    sys.exit(1)
