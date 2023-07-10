#!/bin/bash
set -e

/home/gareth.cabourndavies/test_codes/pycbc/unknown_injections/env_unknown-injections/src/pycbc/bin/workflows/pycbc_make_offline_search_workflow \
--workflow-name gw \
--output-dir output \
--config-files analysis.ini plotting.ini executables.ini injections_minimal.ini injections_unknown.ini \
--config-overrides results_page:output-path:$(pwd)/html
