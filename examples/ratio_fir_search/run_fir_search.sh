#!/bin/bash
set -e

# The FIR/ratio search itself: matched-filters only the coarse templates,
# then reconstructs each fine template's SNR by applying its FIR filter
# (built in make_fir_bank.sh) to the matching coarse template's SNR time
# series. Same data and bank as run_reference_search.sh, so the two
# outputs (fir.hdf here, reference.hdf there) can be compared directly --
# see compare_triggers.py.

export OMP_NUM_THREADS=1

pycbc_inspiral_fir \
    --pad-data 8 \
    --strain-high-pass 15 \
    --sample-rate 2048 \
    --segment-length 512 \
    --segment-start-pad 120 \
    --segment-end-pad 16 \
    --allow-zero-padding \
    --taper-data 1 \
    --psd-estimation median \
    --psd-segment-length 16 \
    --psd-segment-stride 8 \
    --psd-inverse-length 16 \
    --psd-num-segments 63 \
    --low-frequency-cutoff 20 \
    --high-frequency-cutoff 900 \
    --snr-threshold 5.0 \
    --cluster-window 1 \
    --chisq-bins 64 \
    --autochi-stride 2 \
    --autochi-number-points 20 \
    --sgchisq-locations "mtotal>30:20-15,20-30,20-45,20-60,20-75,20-90,20-105,20-120" \
    --signal-consistency-chain snr:6:lower:skip power traditional:10:upper:remove newsnr:6:lower_inc:skip auto sg \
    --channel-name L1:FAKE-STRAIN \
    --gps-start-time 1000000000 \
    --gps-end-time 1000005000 \
    --injection-filter-rejector-chirp-time-window 20000 \
    --output fir.hdf \
    --bank-file fir_full_bank.hdf \
    --fake-strain-seed 0 \
    --fake-strain aLIGOZeroDetHighPower \
    --fake-strain-sample-rate 2048 \
    --fake-strain-flow 12 \
    --verbose \
    --autogating-threshold 15 \
    --autogating-max-iterations 50 \
    --autogating-cluster 0.5 \
    --autogating-width 0.25 \
    --autogating-taper 0.25 \
    --autogating-pad 16
