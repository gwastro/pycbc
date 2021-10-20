set -e

WORKFLOW_NAME=bbh_injections-dynesty
# Set the HTML_DIR to point to your public html page. This is where the results
# page will be written.
HTML_DIR=''
if [ "${HTML_DIR}" == '' ]; then
    echo "Please set an HTML_DIR"
    exit 1
fi
SEED=983124
# Set the number of injections to create. For a full PP test, we suggest using
# 100.
NINJ=10
pycbc_make_inference_inj_workflow \
    --seed ${SEED} \
    --num-injections 10 \
    --config-files workflow_config.ini injections_config.ini \
    --workflow-name ${WORKFLOW_NAME} \
    --config-overrides results_page:output-path:${HTML_DIR}/${WORKFLOW_NAME}
