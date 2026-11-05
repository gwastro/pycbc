#!/bin/bash
set -e

# Not downloading frames from dcc.ligo.org to avoid failures.
# DCC often is not responsive to queries from within the GitHub CI.
# The URL below is the location of the frame files from DCC if you
# wanted to verify they are the same.
# https://dcc.ligo.org/public/0146/P1700341/001/


# List of frame files used by the example
FILES=(
	"H-H1_LOSC_CLN_4_V1-1186740069-3584.gwf"
	"L-L1_LOSC_CLN_4_V1-1186740069-3584.gwf"
	"V-V1_LOSC_CLN_4_V1-1186739813-4096.gwf"
)

# Keep frames alongside this script so CI can save them in its data cache.
cache_dir=$(dirname "$0")

BASE="https://github.com/gwastro/pycbc_data/releases/download/ci-data"

for f in "${FILES[@]}"; do
    # Prefer the working copy
    [ -f "$f" ] && continue
    # wget -nc leaves existing cache files alone
    wget -nv -nc -P "$cache_dir" \
        "${BASE}/$f"
    cp "$cache_dir/$f" .
done
