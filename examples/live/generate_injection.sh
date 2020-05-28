# set parameters for injection
NETWORK_SNR = 20.0
GEOCENTRIC_END_TIME = 1272790440
RA = 45.0
DEC = 45.0
POLARIZATION = 0.0
APPROXIMANT = 'SEOBNRv4'
MASS1 = 25.0
MASS2 = 25.0
INCLINATION = 0.0
TAPER = 'TAPER_START'
INSTRUMENTS = 'H1 L1 V1'
LOW_FREQUENCY_CUTOFF = 10.0


pycbc_generate_hwinj \
  -- network_snr NETWORK_SNR \
  --geocentric_end_time GEOCENTRIC_END_TIME \
  --ra RA \
  --dec DEC \
  --polarization POLARIZATION \
  --approximant APPROXIMANT \
  --mass1 MASS1 \
  --mass2 MASS2 \
  --inclination INCLINATION \
  --taper TAPER \
  --instruments INSTRUMENTS \
  --low_frequency_cutoff LOW_FREQUENCY_CUTOFF
  
  
