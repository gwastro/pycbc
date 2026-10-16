#!/bin/bash
set -e

# Step 1: build the fine template bank (high minimal match) that we want SNR
# time series for, without filtering every template directly -- their SNR is
# reconstructed from automatically selected reference templates via short FIR
# filters (see make_fir_bank.sh).

export OMP_NUM_THREADS=1
LDIR=$(dirname -- "${BASH_SOURCE[0]}")

pycbc_brute_bank --verbose \
    --input-config $LDIR/bank.ini \
    --output-file bank_full.hdf \
    --psd-model aLIGOZeroDetHighPower \
    --low-frequency-cutoff 20.0 \
    --tau0-crawl 100 \
    --tau0-start 1000 \
    --tau0-end 1050 \
    --tau0-threshold 0.5 \
    --minimal-match 0.97 \
    --sample-rate 2048 \
    --tolerance 0.001 \
    --buffer-length 2 \
    --tau0-cutoff-frequency 10 \
    --nprocesses 1 \
    --seed 1
