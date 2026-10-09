#!/bin/sh
export OMP_NUM_THREADS=1
pycbc_inference \
--config-files `dirname "$0"`/lisa_smbhb_relbin.ini \
--output-file lisa_smbhb_ldc_pe.hdf \
--nprocesses 1 \
--force \
--verbose
