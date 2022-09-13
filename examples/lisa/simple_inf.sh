#!/bin/sh
pycbc_inference --verbose \
   --config-files simple_run/model_debug-config.ini \
simple_run/data_debug-config.ini \
simple_run/sampler_debug-config.ini \
   --output-file simple_test.hdf \
   --force