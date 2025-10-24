#!/bin/sh
pycbc_create_injections --verbose \
        --config-files injections_lisa.ini \
        --ninjections 1 \
        --seed 10 \
        --output-file injections_lisa.hdf \
        --variable-params-section variable_params \
        --static-params-section static_params \
        --dist-section prior \
        --force