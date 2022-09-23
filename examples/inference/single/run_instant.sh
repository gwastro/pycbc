pycbc_inference \
--config-file `dirname "$0"`/single_instant.ini \
--nprocesses=8 \
--output-file single_instant.hdf \
--seed 0 \
--force \
--verbose

pycbc_inference_plot_posterior \
--input-file single_instant.hdf \
--output-file single_instant.png \
--parameters distance inclination polarization coa_phase tc ra dec \
--z-arg snr --vmin 31.85 --vmax 32.15 \

