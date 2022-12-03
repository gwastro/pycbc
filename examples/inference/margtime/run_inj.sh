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
 "primary_mass(mass1, mass2) / (1 + redshift(distance)):srcmass1" \
 "secondary_mass(mass1, mass2) / (1 + redshift(distance)):srcmass2" \
 ra dec tc inclination coa_phase polarization distance \
--z-arg snr
