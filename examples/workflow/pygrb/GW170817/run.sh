#!/bin/sh

# Example of a PyGRB workflow analyzing data around the GW170817 event.

set -e

# Download strain data for the GW170817 event and prepare cache files for each
# detector.

for ifo in H-H1 L-L1
do
    ifo_short=${ifo#*-}
    cache_file=strain_${ifo_short}.lcf
    rm -f ${cache_file}
    for time in 1187000320 1187004416 1187008512 1187012608
    do
        file=${ifo}_GWOSC_O2_4KHZ_R1-${time}-4096.gwf
        echo ${ifo_short} GWOSC_O2_4KHZ_R1 ${time} 4096 file://localhost$(readlink -f ${file}) >> strain_${ifo_short}.lcf

        test -f ${file} && continue

        url=https://gwosc.org/archive/data/O2_4KHZ_R1/1186988032/${file}
        echo Downloading ${url}
        curl -O -L --show-error --silent ${url}
    done
done

# Download a basic template bank.

wget -nv -nc https://github.com/gwastro/pycbc-config/raw/master/test/multi_inspiral/gw170817_single_template.hdf

# Create the workflow but do not start it yet.

pycbc_make_offline_grb_workflow \
    --config-files \
        ext_trig.ini \
        analysis.ini \
        executables.ini \
        injections.ini \
        postprocessing.ini \
    --workflow-name pygrb-example \
    --output-dir output \
    --config-overrides \
        results_page:output-path:$(pwd)/html \
    --plan-now

# Start the workflow with the `start` script inside the `output/` directory.
