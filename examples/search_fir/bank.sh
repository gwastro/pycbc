#!/bin/bash
set -e

LDIR=$(dirname -- "${BASH_SOURCE[0]}")
BANK_INI="bank.ini"
if [ ! -f "${BANK_INI}" ]; then
    BANK_INI="${LDIR}/bank.ini"
fi

export OMP_NUM_THREADS=1

pycbc_brute_bank \
    --verbose \
    --output-file bank_raw.hdf \
    --minimal-match 0.97 \
    --tolerance 0.01 \
    --buffer-length 1 \
    --use-trimmed-buffer \
    --sample-rate 2048 \
    --tau0-threshold 0.2 \
    --tau0-crawl 1000 \
    --tau0-start 10 \
    --tau0-end 500 \
    --tau0-cutoff-frequency 20.0 \
    --psd-model aLIGOZeroDetHighPower \
    --seed 1 \
    --low-frequency-cutoff 20.0 \
    --input-config "${BANK_INI}"

mkdir -p bank_full
mkdir -p bank_inj

pycbc_fir_bank \
        --fine-bank bank_raw.hdf \
        --output-file bank_full/bank.hdf \
        --psd-model aLIGOZeroDetHighPower \
        --segment-length 512 \
        --sample-rate 2048 \
        --max-taps 501 6001 \
        --min-match 0.999 0.9999 \
        --n-processes 4 \
        --force \
        --verbose --verbose 
cp bank_full/bank.hdf bank_inj/bank.hdf
