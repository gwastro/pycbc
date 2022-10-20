WORKFLOW_NAME=pycbc_make_faithsim_workflow

./pycbc_make_faithsim_pegasus \
--workflow-name ${WORKFLOW_NAME} \
--config-files ./configuration_pycbc_faithsim.ini \
--output-dir output
