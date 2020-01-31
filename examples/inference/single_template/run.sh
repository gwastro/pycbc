pycbc_inference \
--config-file single.ini \
--nprocesses=4 \
--output-file single.hdf \
--seed 0 \
--force

pycbc_inference_plot_posterior \
--input-file single.hdf \
--output-file single.png \
--z-arg snr \
--parameters distance inclination
