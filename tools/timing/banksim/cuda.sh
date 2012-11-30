pycbc_banksim --signal-file injection0.xml \
--template-file tmplt_bank.xml --asd-file ZERO_DET_high_P.txt \
--match-file cuda.dat --template-approximant="TaylorF2" \
--template-phase-order=7 --template-start-frequency=15 \
--signal-approximant="SpinTaylorT4" --signal-phase-order=7 \
--signal-start-frequency=14 --filter-low-frequency=15 \
--filter-sample-rate=4096 --filter-signal-length=1024 --use-cuda > cuda.out & 
pid=$!
sleep 3600
kill $pid

