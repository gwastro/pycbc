#!/bin/bash

# Debugging: see what our conda environment looks like

echo "DEBUG: inspiral test"
echo "CONDA_PREFIX= $CONDA_PREFIX"
echo "Conda environment is:"
conda list
ls -lh $CONDA_PREFIX/lib/*gomp*
ldd $CONDA_PREFIX/lib/libgomp.so.1

pycbc_inference --verbose \
        --config-files normal2d.ini \
        --output-file normal2d.hdf \
        --nprocesses 2 \
        --seed 10 \
        --force
