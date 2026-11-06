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
#
# Bank storage note: pycbc_fir_bank automatically writes the fine and coarse
# template bank waveforms with compressed representations enabled, and stores
# the FIR filter datasets (taps, actual_tap_count, sigmas) with HDF5 gzip
# compression (level 4).

export OMP_NUM_THREADS=1

pycbc_fir_bank \
    --fine-bank bank_full.hdf \
    --output-file fir_full_bank.hdf \
    --psd-model aLIGOZeroDetHighPower \
    --segment-length ${SEGMENT_LENGTH:-512} \
    --sample-rate 2048 \
    --max-taps ${MAX_TAPS:-1501} \
    --min-match ${MIN_MATCH:-0.99999} \
    --n-processes 1 \
    --verbose "$@"
