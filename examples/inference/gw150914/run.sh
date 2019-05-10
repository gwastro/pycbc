#!/bin/sh
# trigger parameters
TRIGGER_TIME=1126259462.42

# data to use
# the longest waveform covered by the prior must fit in these times
SEARCH_BEFORE=6
SEARCH_AFTER=2

# use an extra number of seconds of data in addition to the data specified
PAD_DATA=8

# PSD estimation options
PSD_ESTIMATION="H1:median L1:median"
PSD_INVLEN=4
PSD_SEG_LEN=16
PSD_STRIDE=8
PSD_DATA_LEN=1024

# sampler parameters
CONFIG_PATH=inference.ini
OUTPUT_PATH=inference.hdf
IFOS="H1 L1"
SAMPLE_RATE=2048
F_HIGHPASS=15
F_MIN=20
PROCESSING_SCHEME=cpu

# the following sets the number of cores to use; adjust as needed to
# your computer's capabilities
NPROCS=10

# get coalescence time as an integer
TRIGGER_TIME_INT=${TRIGGER_TIME%.*}

# start and end time of data to read in
GPS_START_TIME=$((TRIGGER_TIME_INT - SEARCH_BEFORE - PSD_INVLEN))
GPS_END_TIME=$((TRIGGER_TIME_INT + SEARCH_AFTER + PSD_INVLEN))

# start and end time of data to read in for PSD estimation
PSD_START_TIME=$((TRIGGER_TIME_INT - PSD_DATA_LEN/2))
PSD_END_TIME=$((TRIGGER_TIME_INT + PSD_DATA_LEN/2))

# run sampler
# specifies the number of threads for OpenMP
# Running with OMP_NUM_THREADS=1 stops lalsimulation
# from spawning multiple jobs that would otherwise be used
# by inference and cause a reduced runtime.
OMP_NUM_THREADS=1 \
pycbc_inference --verbose \
    --seed 39392 \
    --instruments ${IFOS} \
    --gps-start-time ${GPS_START_TIME} \
    --gps-end-time ${GPS_END_TIME} \
    --channel-name ${CHANNELS} \
    ${FRAMES} \
    --strain-high-pass ${F_HIGHPASS} \
    --pad-data ${PAD_DATA} \
    --psd-estimation ${PSD_ESTIMATION} \
    --psd-start-time ${PSD_START_TIME} \
    --psd-end-time ${PSD_END_TIME} \
    --psd-segment-length ${PSD_SEG_LEN} \
    --psd-segment-stride ${PSD_STRIDE} \
    --psd-inverse-length ${PSD_INVLEN} \
    --sample-rate ${SAMPLE_RATE} \
    --data-conditioning-low-freq ${F_MIN} \
    --config-file ${CONFIG_PATH} \
    --output-file ${OUTPUT_PATH} \
    --processing-scheme ${PROCESSING_SCHEME} \
    --nprocesses ${NPROCS} \
    --force

