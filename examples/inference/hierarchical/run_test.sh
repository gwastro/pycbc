for model in  relbin gaussian margphase; do
    echo "========================"
    echo "With ${model}"
    pycbc_inference \
        --config-files model.ini model-event1_relbin.ini model-event2_${model}.ini prior.ini data.ini emcee.ini \
        --nprocesses 1 \
        --output-file hierarchical-e1_relbin-e2_${model}.hdf \
        --seed 10 \
        --force \
        --verbose
    echo "========================"
done
