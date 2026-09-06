#!/bin/bash
set -e

# Runs the full example end to end: build the two banks, build the FIR
# ratio bank from them, run both searches over the same fake strain data,
# then compare their triggers.

LDIR=$(dirname -- "${BASH_SOURCE[0]}")
cd $LDIR

bash make_banks.sh
bash make_fir_bank.sh
bash run_reference_search.sh
bash run_fir_search.sh
python compare_triggers.py
