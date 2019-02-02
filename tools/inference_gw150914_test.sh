#!/bin/bash

set -e

echo -e "\\n>> [`date`] Testing pycbc_inference on GW150914 data"

cat > gw150914_inference_test.ini <<EOL
[model]
name = gaussian_noise
low-frequency-cutoff = 20

[sampler]
name = emcee_pt
nwalkers = 30
ntemps = 2
niterations = 20
checkpoint-interval = 10

[variable_params]
; waveform parameters that will vary in MCMC
tc =
mass1 =
mass2 =
spin1_a =
spin1_azimuthal =
spin1_polar =
spin2_a =
spin2_azimuthal =
spin2_polar =
distance =
inclination =
polarization =
ra =
dec =
coa_phase =

[static_params]
; waveform parameters that will not change in MCMC
approximant = IMRPhenomPv2
f_lower = 18
f_ref = 20

[prior-tc]
; coalescence time prior
name = uniform
min-tc = 1126259462.32
max-tc = 1126259462.52

[prior-mass1]
name = uniform
min-mass1 = 10.
max-mass1 = 80.

[prior-mass2]
name = uniform
min-mass2 = 10.
max-mass2 = 80.

[prior-spin1_a]
name = uniform
min-spin1_a = 0.0
max-spin1_a = 0.99

[prior-spin1_polar+spin1_azimuthal]
name = uniform_solidangle
polar-angle = spin1_polar
azimuthal-angle = spin1_azimuthal

[prior-spin2_a]
name = uniform
min-spin2_a = 0.0
max-spin2_a = 0.99

[prior-spin2_polar+spin2_azimuthal]
name = uniform_solidangle
polar-angle = spin2_polar
azimuthal-angle = spin2_azimuthal

[prior-distance]
; distance prior
name = uniform_radius
min-distance = 10
max-distance = 1000

[prior-inclination]
; inclination prior
name = sin_angle

[prior-ra+dec]
; sky position prior
name = uniform_sky

[prior-polarization]
; polarization prior
name = uniform_angle

[prior-coa_phase]
; coalescence phase prior
name = uniform_angle


;   Sampling transforms

[sampling_params]
; parameters on the left will be sampled in
; parametes on the right
mass1, mass2 : mchirp, q

[sampling_transforms-mchirp+q]
; inputs mass1, mass2
; outputs mchirp, q
name = mass1_mass2_to_mchirp_q
EOL


FRAMES="H1:../test/data/H-H1_GWOSC_4KHZ_R1-1126257415-4096.gwf L1:../test/data/L-L1_GWOSC_4KHZ_R1-1126257415-4096.gwf"
# trigger parameters
TRIGGER_TIME=1126259462.42

# data to use
# the longest waveform covered by the prior must fit in these times
SEARCH_BEFORE=6
SEARCH_AFTER=2
CHANNELS="H1:GWOSC-4KHZ_R1_STRAIN L1:GWOSC-4KHZ_R1_STRAIN"

# use an extra number of seconds of data in addition to the data specified
PAD_DATA=8

# PSD estimation options
PSD_ESTIMATION="H1:median L1:median"
PSD_INVLEN=4
PSD_SEG_LEN=8
PSD_STRIDE=4
PSD_DATA_LEN=1024

# sampler parameters
IFOS="H1 L1"
SAMPLE_RATE=2048
F_HIGHPASS=15
F_MIN=20

# get coalescence time as an integer
TRIGGER_TIME_INT=${TRIGGER_TIME%.*}

# start and end time of data to read in
GPS_START_TIME=$((TRIGGER_TIME_INT - SEARCH_BEFORE - PSD_INVLEN))
GPS_END_TIME=$((TRIGGER_TIME_INT + SEARCH_AFTER + PSD_INVLEN))

# start and end time of data to read in for PSD estimation
PSD_START_TIME=$((TRIGGER_TIME_INT - PSD_DATA_LEN/2))
PSD_END_TIME=$((TRIGGER_TIME_INT + PSD_DATA_LEN/2))

NPROCS=4
OUTPUT_FILE=gw150914_inference_test.hdf

pycbc_inference --verbose \
    --config-file gw150914_inference_test.ini \
    --output-file ${OUTPUT_FILE} \
    --seed 39392 \
    --instruments ${IFOS} \
    --gps-start-time ${GPS_START_TIME} \
    --gps-end-time ${GPS_END_TIME} \
    --channel-name ${CHANNELS} \
    --frame-files ${FRAMES} \
    --strain-high-pass ${F_HIGHPASS} \
    --pad-data ${PAD_DATA} \
    --psd-estimation ${PSD_ESTIMATION} \
    --psd-start-time ${PSD_START_TIME} \
    --psd-end-time ${PSD_END_TIME} \
    --psd-segment-length ${PSD_SEG_LEN} \
    --psd-segment-stride ${PSD_STRIDE} \
    --psd-inverse-length ${PSD_INVLEN} \
    --sample-rate ${SAMPLE_RATE} \
    --low-frequency-cutoff ${F_MIN} \
    --nprocesses ${NPROCS}

rm -f ${OUTPUT_FILE}

echo -e "\\n>> [`date`] Test GW150914 inference run complete"

exit 0
