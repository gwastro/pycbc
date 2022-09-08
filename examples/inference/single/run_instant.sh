pycbc_inference \
--config-file `dirname "$0"`/single.ini \
--nprocesses=8 \
--output-file single.hdf \
--seed 0 \
--force \
--verbose

pycbc_inference_plot_posterior \
--input-file single.hdf \
--output-file single.png \
--parameters distance inclination polarization coa_phase tc ra dec \
--z-arg snr --vmin 31.85 --vmax 32.15 \

