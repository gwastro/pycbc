#!/bin/sh
pycbc_inference --verbose \
        --config-files normal2d.ini \
        --output-file eryn_normal2d.hdf \
        --nprocesses 1 \
        --seed 10 \
        --force
