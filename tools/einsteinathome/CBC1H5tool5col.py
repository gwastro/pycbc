#!/usr/bin/env python

from __future__ import print_function
import sys
import h5py

f = h5py.File(sys.argv[1], 'r')
detector = f.keys()[0]

template_hashs = f[detector]['template_hash'][:]
chisqs         = f[detector]['chisq'][:]
coa_phases     = f[detector]['coa_phase'][:]
end_times      = f[detector]['end_time'][:]
snrs           = f[detector]['snr'][:]

for trigger in zip(template_hashs, chisqs, coa_phases, end_times, snrs):
    for value in trigger:
        print (" ", value, end="")
    print ("")
