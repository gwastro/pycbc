#!/bin/bash
set -e

/home/gareth.cabourndavies/test_codes/pycbc/snr_timeseries_gdb_upload/env_snr_timeseries_gdb_upload/src/pycbc/bin/workflows/pycbc_make_offline_search_workflow \
--workflow-name gw \
--output-dir output_2 \
--config-files analysis.ini plotting.ini executables.ini injections_minimal.ini \
--config-overrides results_page:output-path:$(pwd)/html \
--cache-file reuse.cache \
--submit-now
