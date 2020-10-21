set -e

WORKFLOW_NAME=inference-dynesty-gw150914_gw170814
# Set the HTML_DIR to point to your public html page. This is where the results
# page will be written.
HTML_DIR=''
if [ "${HTML_DIR}" == '' ]; then
    echo "Please set an HTML_DIR"
    exit 1
fi
SEED=8827
pycbc_make_inference_workflow \
    --seed ${SEED} \
    --config-files workflow_config.ini events.ini \
    --workflow-name ${WORKFLOW_NAME} \
    --config-overrides results_page:output-path:${HTML_DIR}/${WORKFLOW_NAME}
