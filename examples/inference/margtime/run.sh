OMP_NUM_THREADS=1 pycbc_inference \
--config-file `dirname "$0"`/margtime.ini \
--nprocesses 2 \
--processing-scheme mkl \
--output-file marg_150914.hdf \
--seed 0 \
--force \
--verbose

# This reconstructs any marginalized parameters
OMP_NUM_THREADS=1 pycbc_inference_model_stats \
--input-file marg_150914.hdf \
--output-file demarg_150914.hdf \
--nprocesses 2 \
--reconstruct-parameters \
--force \
--verbose

pycbc_inference_plot_posterior \
--input-file demarg_150914.hdf \
--output-file demarg_150914.png \
--parameters \
 "primary_mass(mass1, mass2) / (1 + redshift(distance)):srcmass1" \
 "secondary_mass(mass1, mass2) / (1 + redshift(distance)):srcmass2" \
 ra dec tc inclination coa_phase polarization distance \
--z-arg snr
