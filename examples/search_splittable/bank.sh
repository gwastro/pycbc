#!/bin/bash
set -e

LDIR=`dirname -- "${BASH_SOURCE[0]}"`
echo $LDIR

export OMP_NUM_THREADS=1
pycbc_brute_bank \
--verbose \
--output-file bank.hdf \
--minimal-match 0.95 \
--tolerance .005 \
--buffer-length 2 \
--sample-rate 2048 \
--tau0-threshold 0.5 \
--tau0-crawl 5 \
--tau0-start 0 \
--tau0-end 15 \
--input-config $LDIR/bank.ini \
--psd-model aLIGOZeroDetLowPower \
--seed 1 \
--low-frequency-cutoff 20.0

pycbc_coinc_bank2hdf --bank-file bank.hdf --output-file H1L1V1-BANK2HDF-1186740100-3400.hdf

mkdir bank_full
pycbc_hdf5_splitbank \
--random-sort \
--bank-file H1L1V1-BANK2HDF-1186740100-3400.hdf \
--output-filenames bank_full/H1L1V1-BANK2HDF_SPLITBANK_BANK0_FULL_DATA-1186740100-3400.hdf bank_full/H1L1V1-BANK2HDF_SPLITBANK_BANK1_FULL_DATA-1186740100-3400.hdf \

mkdir bank_inj
pycbc_hdf5_splitbank \
--random-sort  \
--bank-file H1L1V1-BANK2HDF-1186740100-3400.hdf \
--output-filenames bank_inj/H1L1V1-BANK2HDF_SPLITBANK_BANK0_INJECTIONS-1186740100-3400.hdf bank_inj/H1L1V1-BANK2HDF_SPLITBANK_BANK1_INJECTIONS-1186740100-3400.hdf
