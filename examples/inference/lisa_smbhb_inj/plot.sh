#!/bin/sh
pycbc_inference_plot_posterior \
    --input-file lisa_smbhb.hdf \
    --output-file lisa_smbhb_mass_tc.png \
    --z-arg snr --plot-scatter --plot-marginal \
    --plot-contours --contour-color black \
    --parameters \
        'mass1_from_mchirp_q(mchirp,q)':mass1 \
        'mass2_from_mchirp_q(mchirp,q)':mass2 \
        tc \
    --plot-injection-parameters injection_smbhb.hdf
