WORKFLOW_NAME=pycbc_make_faithsim_workflow

../pycbc_make_faithsim_workflow \
    --workflow-name ${WORKFLOW_NAME} \
    --config-files ../configuration_pycbc_faithsim.ini \
    --output-dir output2 \
    --config-overrides \
        workflow:file-retention-level:all_triggers \
        'pegasus_profile:condor|accounting_group:ligo.dev.o4.cbc.uber.pycbcoffline' \
        "pegasus_profile:condor|request_disk:100000" \
        "pegasus_profile:condor|request_memory:2048" \
        "pegasus_profile:condor|when_to_transfer_output:ON_SUCCESS" \
        "pegasus_profile:condor|success_exit_code:0" \
