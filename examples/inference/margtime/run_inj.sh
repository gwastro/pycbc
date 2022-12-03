pycbc_create_injections --verbose \
        --config-files injection.ini \
        --ninjections 1 \
        --seed 10 \
        --output-file injection.hdf \
        --variable-params-section variable_params \
        --static-params-section static_params \
        --dist-section prior \
        --force


OMP_NUM_THREADS=1 pycbc_inference \
--config-file `dirname "$0"`/margtime_inj.ini \
--nprocesses 1 \
--processing-scheme cpu \
--output-file marg_inj.hdf \
--seed 0 \
--force \
--verbose

# This reconstructs any marginalized parameters
OMP_NUM_THREADS=1 pycbc_inference_model_stats \
--input-file marg_inj.hdf \
--output-file demarg_inj.hdf \
--nprocesses 2 \
--reconstruct-parameters \
--force \
--verbose

pycbc_inference_plot_posterior \
--input-file demarg_inj.hdf \
--output-file demarg_inj.png \
--parameters \
 mass1 mass2 ra dec tc inclination coa_phase polarization distance \
--z-arg snr --plot-injection-parameters
