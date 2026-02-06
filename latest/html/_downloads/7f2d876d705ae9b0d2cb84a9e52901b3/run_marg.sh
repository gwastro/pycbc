pycbc_inference \
--config-file `dirname "$0"`/single.ini \
--nprocesses=1 \
--output-file single_marg.hdf \
--seed 0 \
--force \
--verbose

pycbc_inference_plot_posterior \
--input-file single_marg.hdf \
--output-file single_marg.png \
--z-arg snr --vmin 31.85 --vmax 32.15

# This reconstructs any marginalized parameters
# and would be optional if you don't need them or
# have sampled over all parameters directly (see single.ini)
pycbc_inference_model_stats \
--input-file single_marg.hdf \
--output-file single_demarg.hdf \
--nprocesses 1 \
--reconstruct-parameters \
--force \
--verbose

pycbc_inference_plot_posterior \
--input-file single_demarg.hdf \
--output-file single_demarg.png \
--parameters distance inclination polarization coa_phase tc ra dec \
--z-arg snr --vmin 31.85 --vmax 32.15 \

