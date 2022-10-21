OMP_NUM_THREADS=1 python -m cProfile -o log `which pycbc_inference` \
--config-file `dirname "$0"`/reltime.ini \
--nprocesses=2 \
--output-file reltime.hdf \
--seed 0 \
--force \
--verbose
bash p.sh

# This reconstructs any marginalized parameters
# and would be optional if you don't need them or
# have sampled over all parameters directly (see reltime.ini)
OMP_NUM_THREADS=1 python -m cProfile -o log2 `which pycbc_inference_model_stats` \
--input-file reltime.hdf \
--output-file reltime2.hdf \
--nprocesses 2 \
--config-overrides model:marginalize_vector_samples:500 \
--reconstruct-parameters \
--force \
--verbose

pycbc_inference_plot_posterior \
--input-file reltime2.hdf \
--output-file reltime.png \
--parameters distance inclination polarization coa_phase tc ra dec mchirp eta spin \
--z-arg snr
