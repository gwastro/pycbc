#!/bin/bash
set -e

pycbc_make_coinc_search_workflow \
--workflow-name gw \
--output-dir output \
--config-files analysis.ini plotting.ini executables.ini injections_minimal.ini \
--config-overrides results_page:output-path:$(pwd)/html
