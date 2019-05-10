#!/bin/sh
TRIGGER_TIME_INT=1126259462

# sampler parameters
CONFIG_PATH=inference.ini
OUTPUT_PATH=inference.hdf
SEGLEN=8
PSD_INVERSE_LENGTH=4
IFOS="H1 L1"
STRAIN="H1:aLIGOZeroDetHighPower L1:aLIGOZeroDetHighPower"
SAMPLE_RATE=2048
F_MIN=20
PROCESSING_SCHEME=cpu

# the following sets the number of cores to use; adjust as needed to
# your computer's capabilities
NPROCS=10

# start and end time of data to read in
GPS_START_TIME=$((TRIGGER_TIME_INT - SEGLEN))
GPS_END_TIME=$((TRIGGER_TIME_INT + SEGLEN))

# run sampler
# Running with OMP_NUM_THREADS=1 stops lalsimulation
# from spawning multiple jobs that would otherwise be used
# by pycbc_inference and cause a reduced runtime.
OMP_NUM_THREADS=1 \
pycbc_inference --verbose \
    --seed 12 \
    --instruments ${IFOS} \
    --gps-start-time ${GPS_START_TIME} \
    --gps-end-time ${GPS_END_TIME} \
    --psd-model ${STRAIN} \
    --psd-inverse-length ${PSD_INVERSE_LENGTH} \
    --fake-strain ${STRAIN} \
    --fake-strain-seed 44 \
    --strain-high-pass ${F_MIN} \
    --sample-rate ${SAMPLE_RATE} \
    --data-conditioning-low-freq ${F_MIN} \
    --channel-name H1:FOOBAR L1:FOOBAR \
    --injection-file injection.hdf \
    --config-file ${CONFIG_PATH} \
    --output-file ${OUTPUT_PATH} \
    --processing-scheme ${PROCESSING_SCHEME} \
    --nprocesses ${NPROCS} \
    --force

