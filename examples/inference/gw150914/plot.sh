ITER=-1
INPUT_FILE=gwin.hdf.checkpoint

OUTPUT_FILE=scatter_with_stats.png
gwin_plot_posterior \
    --verbose \
    --iteration ${ITER} \
    --input-file ${INPUT_FILE} \
    --output-file ${OUTPUT_FILE} \
    --plot-scatter \
    --plot-marginal \
    --z-arg 'mchirp_from_mass1_mass2(mass1, mass2):mchirp' \
    --parameters tc \
                 'snr_from_loglr(abs(H1_cplx_loglr)):$\rho_H$' \
                 'snr_from_loglr(abs(L1_cplx_loglr)):$\rho_L$' \
                 'snr_from_loglr(abs(loglr)):$\rho$' \

if false; then
OUTPUT_FILE=scatter.png
gwin_plot_posterior \
    --verbose \
    --iteration ${ITER} \
    --input-file ${INPUT_FILE} \
    --output-file ${OUTPUT_FILE} \
    --plot-scatter \
    --plot-marginal \
    --z-arg 'snr_from_loglr(loglr):$\rho$' \
    --parameters "ra*12/pi:$\alpha$ (h)" \
                 "dec*180/pi:$\delta$ (deg)" \
                 "polarization*180/pi:$\psi$ (deg)" \
                 mass1 mass2 spin1_a spin1_azimuthal spin1_polar \
                 spin2_a spin2_azimuthal spin2_polar \
                 "inclination*180/pi:$\iota$ (deg)" distance \
                 "coa_phase*180/pi:$\phi_0$ (deg)" tc
fi
