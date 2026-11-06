#!/bin/bash
set -e

# Make the injection: a 60 solar mass binary with a mass ratio of 1.5
pycbc_create_injections --verbose \
        --config-files injection.ini \
        --ninjections 1 \
        --seed 10 \
        --output-file injection.hdf \
        --variable-params-section variable_params \
        --static-params-section static_params \
        --dist-section prior \
        --force

OMP_NUM_THREADS=1 pycbc_inference \
--config-file `dirname "$0"`/mtotal_marg.ini \
--nprocesses 1 \
--processing-scheme cpu \
--output-file mtotal_marg.hdf \
--seed 0 \
--force \
--verbose

# This reconstructs the marginalized parameters, including the total mass
OMP_NUM_THREADS=1 pycbc_inference_model_stats \
--input-file mtotal_marg.hdf \
--output-file demarg_mtotal.hdf \
--nprocesses 1 \
--reconstruct-parameters \
--force \
--verbose

pycbc_inference_plot_posterior \
--input-file demarg_mtotal.hdf \
--output-file demarg_mtotal.png \
--parameters \
 mtotal q ra dec tc inclination coa_phase polarization distance \
--z-arg snr --plot-injection-parameters
