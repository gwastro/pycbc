pycbc_banksim \
--signal-file ./injection0.hdf \
--template-file ./bank0.hdf \
--psd-model iLIGOModel \
--match-file cpu.dat \
--template-approximant="SpinTaylorT4" \
--template-phase-order=7 \
--template-amplitude-order=0 \
--template-start-frequency=15 \
--signal-approximant="SpinTaylorT4" \
--signal-phase-order=7 \
--signal-amplitude-order=0 \
--signal-start-frequency=15 \
--filter-low-frequency-cutoff=20 \
--filter-sample-rate=4096 \
--filter-signal-length=64

