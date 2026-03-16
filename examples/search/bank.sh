#!/bin/bash
set -e
export OMP_NUM_THREADS=1
pycbc_brute_bank \
--verbose \
--output-file bank.hdf \
--minimal-match 0.95 \
--tolerance .005 \
--buffer-length 2 \
--sample-rate 2048 \
--tau0-threshold 0.5 \
--tau0-crawl 5 \
--tau0-start 0 \
--tau0-end 15 \
--input-config ${BASH_SOURCE%/*}/bank.ini \
--psd-model aLIGOZeroDetLowPower \
--seed 1 \
--low-frequency-cutoff 20.0
