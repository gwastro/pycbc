#!/bin/bash

# Generate injection file
./make_injections_lisa.sh

# Generate Noise without signals frame files
pycbc_generate_mock_data --ifo-list LISA_A LISA_E LISA_T --low-frequency-cutoff LISA_A:1e-4 LISA_E:1e-4 LISA_T:1e-4 --psd-model LISA_A:analytical_psd_lisa_tdi_AE LISA_E:analytical_psd_lisa_tdi_AE LISA_T:analytical_psd_lisa_tdi_T --fake-strain-seed LISA_A:1234 LISA_E:2345 LISA_T:3456 --gps-start-time 1257294808 --gps-end-time 1257294908 --sample-rate 5 --channel-name LISA_A:SIMULATED_STRAIN LISA_E:SIMULATED_STRAIN LISA_T:SIMULATED_STRAIN --tag HLV_NOISE 

# Generate Noise + Signals frame files
pycbc_generate_mock_data --ifo-list LISA_A LISA_E LISA_T --low-frequency-cutoff LISA_A:1e-4 LISA_E:1e-4 LISA_T:1e-4 --psd-model LISA_A:analytical_psd_lisa_tdi_AE LISA_E:analytical_psd_lisa_tdi_AE LISA_T:analytical_psd_lisa_tdi_T --fake-strain-seed LISA_A:1234 LISA_E:2345 LISA_T:3456 --gps-start-time 1257294808 --gps-end-time 1257294908 --sample-rate 5 --channel-name LISA_A:SIMULATED_STRAIN LISA_E:SIMULATED_STRAIN LISA_T:SIMULATED_STRAIN --tag HLV_NOISE_PLUS_SIGNAL --injection-file injections_lisa.hdf

# Generate Signal (with Zeronoise) frame files 
pycbc_generate_mock_data --ifo-list LISA_A LISA_E LISA_T --low-frequency-cutoff LISA_A:1e-4 LISA_E:1e-4 LISA_T:1e-4 --psd-model LISA_A:zeroNoise LISA_E:zeroNoise LISA_T:zeroNoise --fake-strain-seed LISA_A:1234 LISA_E:2345 LISA_T:3456 --gps-start-time 1257294808 --gps-end-time 1257294908 --sample-rate 5 --channel-name LISA_A:SIMULATED_STRAIN LISA_E:SIMULATED_STRAIN LISA_T:SIMULATED_STRAIN --tag HLV_ZERONOISE_SIGNAL --injection-file injections_lisa.hdf

