#!/bin/bash

# Debugging: see what our conda environment looks like

echo "DEBUG: inspiral test"
echo "CONDA_PREFIX= $CONDA_PREFIX"
echo "Conda environment is:"
conda list
ls -lh $CONDA_PREFIX/lib/*gomp*
echo "ldd of libgomp:"
ldd $CONDA_PREFIX/lib/libgomp.so.1
echo "ldd of matchedfilter_cpu.cpython"
ldd $CONDA_PREFIX/lib/python3.11/site-packages/pycbc/filter/matchedfilter_cpu.cpython-311-x86_64-linux-gnu.so

pycbc_inference --verbose \
        --config-files normal2d.ini \
        --output-file normal2d.hdf \
        --nprocesses 2 \
        --seed 10 \
        --force
