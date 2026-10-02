#!/bin/sh
export OMP_NUM_THREADS=1
pycbc_inference \
--config-files `dirname "$0"`/lisa_smbhb_relbin.ini \
--output-file lisa_smbhb_inj_pe.hdf \
--force \
--nprocesses 1 \
--fft-backends fftw \
--verbose

# PLEASE NOTE: This example is currently forcing a FFTW backend because MKL
#              seems to fail for FFT lengths > 2^24. This is fine for most LIGO
#              applications, but an issue for most LISA applications.
