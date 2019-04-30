python bin/pycbc_optimize_cbc \
--trigger-time 1126259462.44 \
--instruments H1 L1 \
--frame-type LOSC  \
--verbose \
--pad-data 8 \
--template-parameters approximant:IMRPhenomD mass1:29 mass2:36 f_lower:25 \
--channel-name L1:LOSC-STRAIN H1:LOSC-STRAIN \
--strain-high-pass 15.0 \
--sample-rate 2048
