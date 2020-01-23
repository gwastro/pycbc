#!/bin/sh
pycbc_inference_plot_posterior \
    --verbose \
    --input-file inference.hdf \
    --output-file posterior-scatter.png \
    --plot-scatter \
    --plot-marginal \
    --z-arg 'loglikelihood:$\log p(h|\vartheta)$' \
    --parameters "ra*12/pi:ra" \
                 "dec*180/pi:dec" \
                 "polarization*180/pi:polarization" \
                 mass1 mass2 spin1_a spin1_azimuthal spin1_polar \
                 spin2_a spin2_azimuthal spin2_polar \
                 "inclination*180/pi:$\iota$ (deg)" distance \
                 'delta_tc:$\Delta t_c$ (s)'
