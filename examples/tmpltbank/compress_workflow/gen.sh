#!/bin/bash
set -e

pycbc_make_bank_compression_workflow \
--workflow-name compress \
--output-dir output \
--config-files \
  compress.ini \
  executables.ini \
--config-overrides \
  results_page:output-path:$(pwd)/html
