pycbc_banksim \
--signal-file ./injection0.xml \
--template-file ./injection0.xml \
--psd iLIGOModel \
--match-file cpu.dat \
--template-approximant="SpinTaylorT4" \
--template-phase-order=7 \
--template-amplitude-order=7 \
--template-start-frequency=40 \
--signal-approximant="SpinTaylorT4" \
--signal-phase-order=7 \
--signal-amplitude-order=7 \
--signal-start-frequency=40 \
--filter-low-frequency=40 \
--filter-sample-rate=4096 \
--filter-signal-length=64

