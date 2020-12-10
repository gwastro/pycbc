#!/bin/bash
set -e

pycbc_dtphase \
--ifos H1 L1 \
--relative-sensitivities .7 1 \
--sample-size 200000 \
--snr-ratio 2.0 \
--seed 10 \
--output-file statHL.hdf \
--smoothing-sigma 1 \
--verbose

pycbc_dtphase \
--ifos L1 V1 \
--relative-sensitivities 1 0.3 \
--sample-size 200000 \
--snr-ratio 2.0 \
--seed 10 \
--output-file statLV.hdf \
--smoothing-sigma 1 \
--verbose

pycbc_dtphase \
--ifos H1 V1 \
--relative-sensitivities .7 .3 \
--sample-size 200000 \
--snr-ratio 2.0 \
--seed 10 \
--output-file statHV.hdf \
--smoothing-sigma 1 \
--verbose


pycbc_dtphase \
--ifos H1 L1 V1 \
--relative-sensitivities .7 1 .3 \
--sample-size 50000 \
--timing-uncertainty .01 \
--snr-ratio 2.0 \
--seed 10 \
--output-file statHLV.hdf \
--smoothing-sigma 1 \
--verbose
