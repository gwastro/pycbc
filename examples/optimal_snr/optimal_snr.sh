#!/bin/bash

# This example shows how to draw a small population of simulated merger signals
# and calculate their optimal SNRs using the `pycbc_create_injections` and
# `pycbc_optimal_snr` commands. It also includes variants that employ the
# `lalapps_inspinj` command and/or the LIGOLW XML format.

set -e

echo Creating some injections in HDF format

pycbc_create_injections \
    --config-files hdf_injection_definition.ini \
    --variable-params-section variable_params \
    --static-params-section static_params \
    --dist-section prior \
    --ninjections 5 \
    --output-file injections.hdf \
    --force \
    --verbose

echo Creating some injections in LIGOLW XML format

lalapps_inspinj \
    --gps-start-time 1276524000 \
    --gps-end-time 1276525000 \
    --time-step 100 \
    --time-interval 0 \
    --l-distr random \
    --i-distr uniform \
    --d-distr uniform \
    --min-distance 10000 \
    --max-distance 500000 \
    --f-lower 15 \
    --m-distr componentMass \
    --min-mass1 10 \
    --max-mass1 40 \
    --min-mass2 10 \
    --max-mass2 40 \
    --waveform IMRPhenomDpseudoFourPN \
    --disable-spin

echo Calculating optimal SNR for HDF injections, saving to LIGOLW XML

pycbc_optimal_snr \
    --input-file injections.hdf \
    --output-file injections_optimal_snr_hdf2xml.xml.gz \
    --snr-columns \
        H1:alpha1 \
        L1:alpha2 \
        V1:alpha3 \
    --psd-model aLIGOZeroDetHighPower

echo Calculating optimal SNR for HDF injections, saving to HDF

pycbc_optimal_snr \
    --input-file injections.hdf \
    --output-file injections_optimal_snr_hdf2hdf.hdf \
    --snr-columns \
        H1:optimal_snr_lho \
        L1:optimal_snr_llo \
        V1:optimal_snr_virgo \
    --psd-model aLIGOZeroDetHighPower

echo Calculating optimal SNR for LIGOLW XML injections, saving to LIGOLW XML

pycbc_optimal_snr \
    --input-file HL-INJECTIONS_1-1276524000-1000.xml \
    --output-file injections_optimal_snr_xml2xml.xml.gz \
    --snr-columns \
        H1:alpha1 \
        L1:alpha2 \
        V1:alpha3 \
    --psd-model aLIGOZeroDetHighPower
