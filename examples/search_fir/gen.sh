#!/bin/bash
set -e

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

pycbc_create_injections \
    --verbose \
    --config-files "${SCRIPT_DIR}/injections.ini" \
    --gps-start-time 1186725888 \
    --gps-end-time 1186754560 \
    --time-step 100 \
    --time-window 40 \
    --seed 123407 \
    --output-file injections.hdf \
    --force

pycbc_make_offline_search_workflow \
    --workflow-name gw \
    --output-dir output \
    --config-files analysis.ini plotting.ini executables.ini injections_minimal.ini \
    --config-overrides results_page:output-path:$(pwd)/html \
    --plan-now
