WORKFLOW_NAME=example_faithsim_workflow

pycbc_make_faithsim_workflow \
    --workflow-name ${WORKFLOW_NAME} \
    --config-files pegasus_workflow_conf.ini \
    --output-dir output
