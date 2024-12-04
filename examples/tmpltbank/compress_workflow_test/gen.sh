#!/bin/bash
set -e

pycbc_make_bank_compression_workflow \
  --config-files compress.ini executables.ini \
  --workflow-name compress \
  --output-dir output \
  --config-overrides \
    results_page:output-path:$(pwd)/html \

