#!/bin/sh
pycbc_inference_plot_posterior --verbose \
        --input-file eryn_normal2d.hdf \
        --output-file posterior-eryn_normal2d.png \
        --plot-scatter \
        --plot-contours \
        --plot-marginal \
        --z-arg 'loglikelihood:$\log p(h|\vartheta)$'
