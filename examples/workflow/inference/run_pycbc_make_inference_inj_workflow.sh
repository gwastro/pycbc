#! /bin/bash

set -e
workflow_name=emcee_pt_pptest
output_dir=run

export OMP_NUM_THREADS=1

pycbc_make_inference_inj_workflow \
    --output-dir ${output_dir} \
    --workflow-name ${workflow_name} \
    --data-type simulated_data \
    --output-file ${workflow_name}.dax \
    --inference-config-file inference.ini \
    --config-files workflow_config.ini \
    --config-overrides results_page:output-path:${PWD}/${output_dir}/results_html

echo "Now run:"
echo "cd ${output_dir}"
echo "pycbc_submit_dax --no-create-proxy --no-grid --dax ${workflow_name}.dax --accounting-group ${accounting_group} --enable-shared-filesystem"
