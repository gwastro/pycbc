#!/bin/sh
SAMPLERS="emcee_stub.ini emcee_pt_stub.ini dynesty_stub.ini ultranest_stub.ini epsie_stub.ini nessai_stub.ini snowline_stub.ini"
PLOT_INPUTS="emcee_stub.ini.hdf:emcee emcee_pt_stub.ini.hdf:emcee_pt dynesty_stub.ini.hdf:dynesty ultranest_stub.ini.hdf:ultranest epsie_stub.ini.hdf:epsie nessai_stub.ini.hdf:nessai snowline_stub.ini.hdf:snowline"

if [ "$(uname)" != "Darwin" ]; then
    SAMPLERS="cpnest_stub.ini $SAMPLERS"
    PLOT_INPUTS="$PLOT_INPUTS cpnest_stub.ini.hdf:cpnest"
fi

for f in $SAMPLERS; do
        echo $f
        pycbc_inference \
        --config-files `dirname $0`/simp.ini `dirname $0`/$f \
        --output-file $f.hdf \
        --nprocesses 1 \
        --seed 10 \
        --force
done

pycbc_inference_plot_posterior --input-file $PLOT_INPUTS \
--output-file sample.png \
--plot-contours \
--plot-marginal \
--no-contour-labels \
--no-marginal-titles
