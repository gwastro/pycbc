cd output

pycbc_submit_dax \
    --append-pegasus-property 'pegasus_profile:condor|accounting_group:ligo.dev.o4.cbc.uber.pycbcoffline' \
    --append-pegasus-property 'pegasus.transfer.bypass.input.staging=true' \
    --append-pegasus-property "pegasus_profile:condor|request_disk:10000000" \
    --append-pegasus-property "pegasus_profile:condor|request_gpus:1" \
    --append-pegasus-property "pegasus_profile:condor|request_memory:4096" \
    --append-pegasus-property "pegasus_profile:condor|when_to_transfer_output:ON_SUCCESS" \
    --append-pegasus-property "pegasus_profile:condor|success_exit_code:0" \
