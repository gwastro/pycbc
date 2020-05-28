

pycbc_generate_hwinj \
  -- network_snr 20.0 \
  --geocentric_end_time 1272790440 \
  --ra 45.0 \
  --dec 45.0 \
  --polarization 0.0 \
  --approximant SEOBNRv4 \
  --mass1 25.0 \
  --mass2 25.0 \
  --inclination 0.0 \
  --taper TAPER_START \
  --instruments H1 L1 V1 \
  --low_frequency_cutoff 10 \
  --waveform_low_frequency_cutoff 10
  
  
