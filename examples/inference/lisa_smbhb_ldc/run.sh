#!/bin/sh
pycbc_inference \
--config-files `dirname "$0"`/lisa_smbhb_relbin.ini \
--output-file lisa_smbhb.hdf \
--force \
--verbose
