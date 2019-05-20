#! /bin/bash

set -e

# name of the workflow
WORKFLOW_NAME="event"

# path to output dir
OUTPUT_DIR=output

# input configuration files
CONFIG_PATH=workflow.ini
INFERENCE_CONFIG_PATH=inference.ini

# directory that will be populated with HTML pages
HTML_DIR=${HOME}/public_html/inference_test

# run workflow generator on specific GPS end time
pycbc_make_inference_workflow --workflow-name ${WORKFLOW_NAME} \
    --config-files ${CONFIG_PATH} \
    --inference-config-file ${INFERENCE_CONFIG_PATH} \
    --output-dir ${OUTPUT_DIR} \
    --output-file ${WORKFLOW_NAME}.dax \
    --output-map ${OUTPUT_MAP_PATH} \
    --gps-end-time ${GPS_END_TIME} \
    --config-overrides workflow:start-time:$((${GPS_END_TIME}-2)) \
                       workflow:end-time:$((${GPS_END_TIME}+2)) \
                       workflow-inference:data-seconds-before-trigger:2 \
                       workflow-inference:data-seconds-after-trigger:2 \
                       inference:psd-start-time:$((${GPS_END_TIME}-300)) \
                       inference:psd-end-time:$((${GPS_END_TIME}+1748)) \
                       results_page:output-path:${HTML_DIR} \
                       results_page:analysis-subtitle:${WORKFLOW_NAME}
