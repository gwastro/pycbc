python bin/pycbc_optimize_cbc \
--trigger-time 1126259462 \
--instruments H1 L1 \
--frame-type LOSC  \
--verbose \
--pad-data 8 \
--template-parameters approximant:IMRPhenomD mass1:10 mass2:10 f_lower:25 \
--channel-name L1:LOSC-STRAIN H1:LOSC-STRAIN \
--strain-high-pass 15.0 \
--sample-rate 2048
