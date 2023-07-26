pycbc_inference \
    --config-files model.ini model-event1_relbin.ini model-event2_relbin.ini prior.ini data.ini dynesty.ini \
    --nprocesses 8 \
    --output-file hierarchical.hdf \
    --seed 10 \
    --force \
    --verbose
