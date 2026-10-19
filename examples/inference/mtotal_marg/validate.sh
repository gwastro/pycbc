#!/bin/bash
# Runs the same analysis with the base model sampling the total mass directly,
# so the two can be compared. Run run.sh first.
set -e

OMP_NUM_THREADS=1 pycbc_inference \
--config-file `dirname "$0"`/validation.ini \
--nprocesses 1 \
--processing-scheme cpu \
--output-file validation.hdf \
--seed 0 \
--force \
--verbose

OMP_NUM_THREADS=1 pycbc_inference_model_stats \
--input-file validation.hdf \
--output-file demarg_validation.hdf \
--nprocesses 1 \
--reconstruct-parameters \
--force \
--verbose

# Overlay the total mass marginalized run on the one that sampled it directly
pycbc_inference_plot_posterior \
--input-file demarg_mtotal.hdf demarg_validation.hdf \
--input-file-labels "marginalized_mtotal" "marginalized_time" \
--output-file comparison.png \
--parameters mtotal q inclination distance \
--plot-injection-parameters
