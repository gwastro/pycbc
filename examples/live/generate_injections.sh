#major work in progress

inj_snr=15
inj_time=1272790440
inj_mass1=7.2797217
inj_mass2=6.9102380
inj_spin1z=0.7189988
inj_spin2z=0.1991984
    
pycbc_generate_hwinj \
    --network-snr $inj_snr \
    --ra 45.0 \
    --dec 45.0 \
    --polarization 0.0 \
    --approximant SEOBNRv4 \
    --mass1 $inj_mass1 \
    --mass2 $inj_mass2 \
    --spin1z $inj_spin1z \
    --spin2z $inj_spin2z \
    --inclination 0.0 \
    --taper TAPER_START \
    --waveform-low-frequency-cutoff 10 \
    --geocentric-end-time 1272790440 \
    --instruments H1 L1 V1 \
    --low-frequency-cutoff 10 \
    --sample-rate H1:16384 L1:16384 V1:16384 \
    --gps-end-time 1272790500 \
    --gps-start-time 1272790000 \
    --psd-model H1:aLIGOMidLowSensitivityP1200087 L1:aLIGOMidLowSensitivityP1200087 V1:AdVEarlyLowSensitivityP1200087
