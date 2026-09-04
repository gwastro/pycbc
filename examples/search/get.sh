#!/bin/bash
set -e

# Not downloading frames from dcc.ligo.org to avoid failures.
# DCC often is not responsive to queries from within the GitHub CI.
# The commented commands below are how to get the frame files from DCC if you
# wanted to verify they are the same.

#wget -nv https://dcc.ligo.org/public/0146/P1700341/001/H-H1_LOSC_CLN_4_V1-1186740069-3584.gwf
#wget -nv https://dcc.ligo.org/public/0146/P1700341/001/L-L1_LOSC_CLN_4_V1-1186740069-3584.gwf
#wget -nv https://dcc.ligo.org/public/0146/P1700341/001/V-V1_LOSC_CLN_4_V1-1186739813-4096.gwf

# Files are hosted as GitHub release assets (not Git LFS) so downloads are not
# rate-limited. --retry* guards against transient GitHub/CDN errors.
BASE=https://github.com/gwastro/pycbc_data/releases/download/ci-data
WGET="wget -nv --tries=5 --retry-connrefused --waitretry=5 --timeout=30 --retry-on-http-error=429,500,502,503,504"
${WGET} ${BASE}/H-H1_LOSC_CLN_4_V1-1186740069-3584.gwf
${WGET} ${BASE}/L-L1_LOSC_CLN_4_V1-1186740069-3584.gwf
${WGET} ${BASE}/V-V1_LOSC_CLN_4_V1-1186739813-4096.gwf
