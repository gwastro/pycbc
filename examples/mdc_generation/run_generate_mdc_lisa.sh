#!/bin/bash

# Generate injection file
./make_injections_lisa.sh

# Generate noise without signals frame files
pycbc_generate_mock_data \
    --ifo-list LISA_A LISA_E LISA_T \
    --low-frequency-cutoff LISA_A:1e-4 LISA_E:1e-4 LISA_T:1e-4 \
    --psd-model LISA_A:analytical_psd_lisa_tdi_AE \
                LISA_E:analytical_psd_lisa_tdi_AE \
                LISA_T:analytical_psd_lisa_tdi_T \
    --fake-strain-seed LISA_A:100 LISA_E:150 LISA_T:200 \
    --gps-start-time 0 \
    --gps-end-time 31536000 \
    --sample-rate 0.2 \
    --len-arm 2.5e9 \
    --acc-noise-level 2.4e-15 \
    --oms-noise-level 7.9e-12 \
    --tdi '2.0' \
    --channel-name LISA_A:SIMULATED_STRAIN \
                   LISA_E:SIMULATED_STRAIN \
                   LISA_T:SIMULATED_STRAIN \
    --tag LISA_NOISE \
    --fake-strain-filter-duration 31536000

# Generate noise + signals frame files
    pycbc_generate_mock_data \
        --ifo-list LISA_A LISA_E LISA_T \
        --low-frequency-cutoff LISA_A:1e-4 LISA_E:1e-4 LISA_T:1e-4 \
        --psd-model LISA_A:analytical_psd_lisa_tdi_AE \
                    LISA_E:analytical_psd_lisa_tdi_AE \
                    LISA_T:analytical_psd_lisa_tdi_T \
        --fake-strain-seed LISA_A:100 LISA_E:150 LISA_T:200 \
        --gps-start-time 0 \
        --gps-end-time 31536000 \
        --sample-rate 0.2 \
        --len-arm 2.5e9 \
        --acc-noise-level 2.4e-15 \
        --oms-noise-level 7.9e-12 \
        --tdi '2.0' \
        --channel-name LISA_A:SIMULATED_STRAIN \
                       LISA_E:SIMULATED_STRAIN \
                       LISA_T:SIMULATED_STRAIN \
        --injection-file injection_lisa_1_smbhb.hdf \
        --tag LISA_NOISE_SIGNAL \
        --fake-strain-filter-duration 31536000

# Generate Signal (with Zeronoise) frame files 
    pycbc_generate_mock_data \
        --ifo-list LISA_A LISA_E LISA_T \
        --low-frequency-cutoff LISA_A:1e-4 LISA_E:1e-4 LISA_T:1e-4 \
        --psd-model LISA_A:zeroNoise LISA_E:zeroNoise LISA_T:zeroNoise \
        --fake-strain-seed LISA_A:100 LISA_E:150 LISA_T:200 \
        --gps-start-time 0 \
        --gps-end-time 31536000 \
        --sample-rate 0.2 \
        --len-arm 2.5e9 \
        --acc-noise-level 2.4e-15 \
        --oms-noise-level 7.9e-12 \
        --tdi '2.0' \
        --channel-name LISA_A:SIMULATED_STRAIN \
                       LISA_E:SIMULATED_STRAIN \
                       LISA_T:SIMULATED_STRAIN \
        --injection-file injection_lisa_1_smbhb.hdf \
        --tag LISA_ZERONOISE_SIGNAL \
        --fake-strain-filter-duration 31536000
