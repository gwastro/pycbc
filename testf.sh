python bin/pycbc_optimize_cbc \
--trigger-time 1187008882.43 \
--instruments H1 L1 \
--frame-type LOSC  \
--verbose \
--pad-data 8 \
--autogating-pad 1.0 \
--autogating-width 0.25 \
--autogating-taper 0.25 \
--autogating-threshold 50 \
--template-parameters mass1:1.37 mass2:1.37 f_lower:25 \
--channel-name L1:GWOSC-16KHZ_R1_STRAIN H1:GWOSC-16KHZ_R1_STRAIN \
--strain-high-pass 15.0 \
--sample-rate 2048

#--trigger-time 1126259462.44 \
