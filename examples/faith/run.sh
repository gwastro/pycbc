#!/bin/bash
pycbc_faithsim \
--param-file injection0.xml \
--psd iLIGOModel \
--match-file match.dat \
--waveform1-approximant="SpinTaylorT4" \
--waveform1-start-frequency=40 \
--waveform2-approximant="SpinTaylorT4" \
--waveform2-start-frequency=40 \
--filter-low-frequency=40 \
--filter-sample-rate=4096 \
--filter-waveform-length=64

