#! /bin/bash

set -e

# name of the workflow
WORKFLOW_NAME="search_triggers"

# path to output dir
OUTPUT_DIR=output

# input configuration files
CONFIG_PATH=workflow.ini
INFERENCE_CONFIG_PATH=inference.ini

# directory that will be populated with HTML pages
HTML_DIR=${HOME}/public_html/inference_test

# run workflow generator on triggers from workflow
pycbc_make_inference_workflow --workflow-name ${WORKFLOW_NAME} \
    --config-files ${CONFIG_PATH} \
    --inference-config-file ${INFERENCE_CONFIG_PATH} \
    --output-dir ${OUTPUT_DIR} \
    --output-file ${WORKFLOW_NAME}.dax \
    --output-map ${OUTPUT_MAP_PATH} \
    --bank-file ${BANK_PATH} \
    --statmap-file ${STATMAP_PATH} \
    --single-detector-triggers ${SNGL_H1_PATHS} ${SNGL_L1_PATHS}
    --config-overrides workflow:start-time:${WORKFLOW_START_TIME} \
                       workflow:end-time:${WORKFLOW_END_TIME} \
                       workflow-inference:data-seconds-before-trigger:8 \
                       workflow-inference:data-seconds-after-trigger:8 \
                       results_page:output-path:${HTML_DIR} \
                       results_page:analysis-subtitle:${WORKFLOW_NAME}
