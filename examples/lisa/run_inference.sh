#!/bin/sh
pycbc_inference --verbose \
   --config-files configs/0/model_debug-config.ini \
configs/0/data_debug-config.ini \
configs/0/sampler_debug-config.ini \
   --output-file inference_test.hdf \
   --force