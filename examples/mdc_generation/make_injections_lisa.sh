#!/bin/sh
pycbc_create_injections --verbose \
        --config-files injection_lisa_1_smbhb.ini \
        --ninjections 1 \
        --seed 10 \
        --output-file injection_lisa_1_smbhb.hdf \
        --variable-params-section variable_params \
        --static-params-section static_params \
        --dist-section prior \
        --force
