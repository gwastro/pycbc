WORKFLOW_NAME=example_faithsim_workflow

pycbc_make_faithsim_workflow \
    --workflow-name ${WORKFLOW_NAME} \
    --config-files faithsim_workflow_config.ini \
    --output-dir output \
    --submit-now
