pycbc_inference_plot_posterior \
    --input-file lisa_smbhb.hdf \
    --output-file lisa_smbhb_mass_tc_0.png \
    --z-arg snr --plot-scatter --plot-marginal \
    --plot-contours --contour-color black \
    --parameters \
        'mass1_from_mchirp_q(mchirp,q)':mass1 \
        'mass2_from_mchirp_q(mchirp,q)':mass2 \
        tc \
    --expected-parameters \
        'mass1_from_mchirp_q(mchirp,q)':1015522.4376 \
        'mass2_from_mchirp_q(mchirp,q)':796849.1091 \
        tc:4799624.274911478 \
