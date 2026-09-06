#!/bin/bash
set -e

# Step 1: build the two template banks the FIR/ratio method needs:
#  - a loose "coarse" bank (low minimal match): these are the templates that
#    actually get matched-filtered against the data.
#  - a tight "fine" bank (high minimal match): the templates we want SNR
#    time series for, but never filter directly -- their SNR is instead
#    reconstructed from a coarse template's SNR via a short FIR filter
#    (see make_fir_bank.sh).

export OMP_NUM_THREADS=1
LDIR=$(dirname -- "${BASH_SOURCE[0]}")

pycbc_brute_bank --verbose \
    --input-config $LDIR/bank.ini \
    --output-file bank_coarse.hdf \
    --psd-model aLIGOZeroDetHighPower \
    --low-frequency-cutoff 20.0 \
    --tau0-crawl 100 \
    --tau0-start 1000 \
    --tau0-end 1050 \
    --tau0-threshold 50 \
    --minimal-match 0.3 \
    --sample-rate 2048 \
    --tolerance 0.01 \
    --buffer-length 2 \
    --tau0-cutoff-frequency 10 \
    --nprocesses 1 \
    --seed 1

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
