#!/bin/sh
pycbc_create_injections --verbose \
        --config-files injection.ini \
        --ninjections 1 \
        --seed 10 \
        --output-file injection.hdf \
        --variable-params-section variable_params \
        --static-args-section static_params \
        --dist-section prior
