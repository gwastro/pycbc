#!/bin/bash
set -e

pycbc_make_offline_search_workflow \
--workflow-name gw \
--output-dir output_3 \
--config-files analysis.ini plotting.ini executables.ini injections_minimal.ini \
--config-overrides results_page:output-path:$(pwd)/html \
--submit-now \
--cache-file reuse.cache
