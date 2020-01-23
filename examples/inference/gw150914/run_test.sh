#!/bin/sh

# configuration files
PRIOR_CONFIG=../priors/gw150914_like.ini
DATA_CONFIG=data.ini
SAMPLER_CONFIG=../samplers/emcee_pt-gw150914_like.ini

OUTPUT_PATH=inference.hdf

# the following sets the number of cores to use; adjust as needed to
# your computer's capabilities
NPROCS=10

# run sampler
# Running with OMP_NUM_THREADS=1 stops lalsimulation
# from spawning multiple jobs that would otherwise be used
# by pycbc_inference and cause a reduced runtime.
OMP_NUM_THREADS=1 \
pycbc_inference --verbose \
    --seed 1897234 \
    --config-file ${PRIOR_CONFIG} ${DATA_CONFIG} ${SAMPLER_CONFIG} \
    --output-file ${OUTPUT_PATH} \
    --nprocesses ${NPROCS} \
    --config-delete "sampler:effective-nsamples" \
                    "sampler:max-samples-per-chain" \
                    "data:frame-files" \
    --config-overrides sampler:ntemps:2 \
                       sampler:nwalkers:30 \
                       sampler:niterations:20 \
                       sampler:checkpoint-interval:10 \
                       "data:frame-type:H1:LOSC L1:LOSC" \
                       "data:channel-name:H1:LOSC-STRAIN L1:LOSC-STRAIN" \
    --force
