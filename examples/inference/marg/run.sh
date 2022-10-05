OMP_NUM_THREADS=1 python -m cProfile -o log `which pycbc_inference` \
--config-file `dirname "$0"`/single.ini \
--nprocesses=8 \
--output-file single.hdf \
--seed 0 \
--force \
--verbose
bash p.sh
# This reconstructs any marginalized parameters
# and would be optional if you don't need them or
# have sampled over all parameters directly (see single.ini)
OMP_NUM_THREADS=1 python -m cProfile -o log2 `which pycbc_inference_model_stats` \
--input-file single.hdf \
--output-file single2.hdf \
--nprocesses 8 \
--reconstruct-parameters \
--force \
--verbose

#pycbc_inference_plot_posterior \
#--input-file single2.hdf \
#--output-file single.png \
#--parameters distance inclination polarization coa_phase tc ra dec \
#--z-arg snr --vmin 31.85 --vmax 32.15 \

pycbc_inference_plot_posterior \
--input-file single2.hdf \
--output-file single.png \
--parameters polarization coa_phase distance inclination tc ra dec mchirp eta \
--z-arg snr
