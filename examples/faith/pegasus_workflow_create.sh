WORKFLOW_NAME=pycbc_make_faithsim_workflow

../pycbc_make_faithsim_workflow_finalize \
    --workflow-name ${WORKFLOW_NAME} \
    --config-files ../pycbc_faithsim_configuration.ini \
    --output-dir output
