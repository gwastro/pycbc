#!/bin/bash
set -e

# List of frame files used by the example
FILES=(
	"H-H1_LOSC_CLN_4_V1-1186740069-3584.gwf"
	"L-L1_LOSC_CLN_4_V1-1186740069-3584.gwf"
	"V-V1_LOSC_CLN_4_V1-1186739813-4096.gwf"
)

for f in "${FILES[@]}"; do
	# If file is already in current working directory, skip
	if [ -f "./$f" ]; then
		echo "Found $f in working directory; skipping download."
		continue
	fi

	# If file is present in examples/search (i.e. the cache), then copy it
	if [ -f "examples/search/$f" ]; then
		echo "Found $f in examples/search (cache); copying."
		cp "examples/search/$f" ./
		continue
	fi

	# Otherwise, download from the DCC
	echo "Downloading $f from DCC..."
	wget -nv "https://dcc.ligo.org/public/0146/P1700341/001/$f"
done
