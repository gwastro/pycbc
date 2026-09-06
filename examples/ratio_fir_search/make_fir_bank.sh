#!/bin/bash
set -e

# Step 2: for every "fine" template, find its best-matching "coarse"
# reference template and least-squares fit a short FIR filter that
# approximates their frequency-domain ratio (fine / coarse). At search
# time, applying this filter to the coarse template's SNR time series
# reconstructs an approximation of the fine template's own SNR, without
# ever matched-filtering the fine template directly.
#
# --min-match is a target on the fit's own quality (a template overlap
# check between the fine template and the ratio-filter reconstruction),
# not a general property of the method -- it's high here mainly to make
# comparing against a plain pycbc_inspiral run easier to interpret.

export OMP_NUM_THREADS=1

pycbc_fir_bank \
    --coarse-bank bank_coarse.hdf \
    --fine-bank bank_full.hdf \
    --output-file fir_full_bank.hdf \
    --psd-model aLIGOZeroDetHighPower \
    --sample-rate 2048 \
    --f-low 20.0 \
    --n-taps 251 \
    --max-taps 1501 \
    --tap-step 100 \
    --min-match 0.99999 \
    --delta-f 0.5 \
    --regularization-type ridge \
    --regularization-magnitude .0001 \
    --ratio-cap 1000 \
    --max-tap-value 100000 \
    --search-depth 1 \
    --n-processes 1 \
    --verbose
