#!/bin/bash

# Generate injection file
./make_injections_hlv.sh

# Generate Noise without signals frame files
pycbc_generate_mock_data --ifo-list H1 L1 V1 --low-frequency-cutoff H1:20 L1:20 V1:20 --psd-model H1:aLIGOZeroDetHighPower L1:aLIGOZeroDetHighPower V1:aLIGOZeroDetHighPower --fake-strain-seed H1:1234 L1:2345 V1:3456 --gps-start-time 1257294808 --gps-end-time 1257294908 --sample-rate 2048 --channel-name H1:SIMULATED_STRAIN L1:SIMULATED_STRAIN V1:SIMULATED_STRAIN --tag HLV_NOISE 

# Generate Noise + Signals frame files
pycbc_generate_mock_data --ifo-list H1 L1 V1 --low-frequency-cutoff H1:20 L1:20 V1:20 --psd-model H1:aLIGOZeroDetHighPower L1:aLIGOZeroDetHighPower V1:aLIGOZeroDetHighPower --fake-strain-seed H1:1234 L1:2345 V1:3456 --gps-start-time 1257294808 --gps-end-time 1257294908 --sample-rate 2048 --channel-name H1:SIMULATED_STRAIN L1:SIMULATED_STRAIN V1:SIMULATED_STRAIN --tag HLV_NOISE_SIGNAL --injection-file injection_hlv_10_bbh.hdf

# Generate Signal (with Zeronoise) frame files 
pycbc_generate_mock_data --ifo-list H1 L1 V1 --low-frequency-cutoff H1:20 L1:20 V1:20 --psd-model H1:zeroNoise L1:zeroNoise V1:zeroNoise --fake-strain-seed H1:1234 L1:2345 V1:3456 --gps-start-time 1257294808 --gps-end-time 1257294908 --sample-rate 2048 --channel-name H1:SIMULATED_STRAIN L1:SIMULATED_STRAIN V1:SIMULATED_STRAIN --tag HLV_ZERONOISE_SIGNAL --injection-file injection_hlv_10_bbh.hdf
