#!/bin/bash
set -e

pycbc_make_offline_search_workflow \
--workflow-name gw \
--output-dir output \
--config-files analysis.ini ../search/plotting.ini ../search/executables.ini ../search/injections_minimal.ini \
--config-overrides results_page:output-path:$(pwd)/html \
--submit-now
