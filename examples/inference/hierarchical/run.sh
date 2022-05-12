pycbc_inference \
    --config-files model.ini prior.ini data.ini sampler.ini \
    --nprocesses 4 \
    --output-file hierarchical.hdf \
    --seed 10 \
    --force \
    --verbose
