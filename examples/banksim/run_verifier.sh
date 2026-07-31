#!/bin/bash
WORKFLOW_NAME=example_bank_verifier

pycbc_make_bank_verifier_workflow \
    --workflow-name ${WORKFLOW_NAME} \
    --config-files bank_verifier_config.ini \
    --output-dir output \
    --submit-now
