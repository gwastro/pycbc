H1=`realpath H*.gwf`
L1=`realpath L*.gwf`
V1=`realpath V*.gwf`

pycbc_make_inference_workflow \
    --seed 0 \
    --config-files workflow.ini \
    --workflow-name gw \
    --config-overrides \
    inference:config-overrides:data:"frame-files:'H1:${H1} L1:${L1} V1:${V1}'" \
    results_page:output-path:${PWD}/html

# The above adds an option to the workflow.ini to set a config override for the inference jobs
# the inference job then sets where it gets the frame files to our chosen location
