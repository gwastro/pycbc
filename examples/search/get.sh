#!/bin/bash
set -e

# Not downloading frames from dcc.ligo.org to avoid failures.
# DCC often is not responsive to queries from within the GitHub CI.
# The commented commands below are how to get the frame files from DCC if you
# wanted to verify they are the same.

#wget -nv https://dcc.ligo.org/public/0146/P1700341/001/H-H1_LOSC_CLN_4_V1-1186740069-3584.gwf
#wget -nv https://dcc.ligo.org/public/0146/P1700341/001/L-L1_LOSC_CLN_4_V1-1186740069-3584.gwf
#wget -nv https://dcc.ligo.org/public/0146/P1700341/001/V-V1_LOSC_CLN_4_V1-1186739813-4096.gwf

wget -nv https://media.githubusercontent.com/media/gwastro/pycbc_data/master/H-H1_LOSC_CLN_4_V1-1186740069-3584.gwf
wget -nv https://media.githubusercontent.com/media/gwastro/pycbc_data/master/L-L1_LOSC_CLN_4_V1-1186740069-3584.gwf
wget -nv https://media.githubusercontent.com/media/gwastro/pycbc_data/master/V-V1_LOSC_CLN_4_V1-1186739813-4096.gwf
