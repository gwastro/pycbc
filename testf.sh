python bin/pycbc_optimize_cbc \
\
--trigger-time 1187008882.43 \
--instruments H1 L1 \
--template-parameters mass1:1.37 mass2:1.37 f_lower:25 \
\
--autogating-pad 1.0 \
--autogating-width 0.25 \
--autogating-taper 0.25 \
--autogating-threshold 50 \
\
--hdf-store L1:L1.hdf H1:H1.hdf  \
--channel-name L1:GWOSC-16KHZ_R1_STRAIN H1:GWOSC-16KHZ_R1_STRAIN \
--low-frequency-cutoff 20.0 \
--strain-high-pass 15.0 \
--verbose
