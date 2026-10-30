#!/bin/sh
pycbc_inference --verbose \
        --config-files normal2d.ini \
        --output-file normal2d.hdf \
        --nprocesses 2 \
        --seed 10 \
        --force
