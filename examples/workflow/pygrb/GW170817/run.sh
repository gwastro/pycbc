#!/bin/sh

# Example of a PyGRB workflow analyzing data around the GW170817 event.

set -e

# Download strain data for the GW170817 event.

for ifo in H-H1 L-L1
do
    file=${ifo}_LOSC_CLN_4_V1-1187007040-2048.gwf
    test -f ${file} && continue

    # Not downloading frames from dcc.ligo.org to avoid failures.
    # DCC often is not responsive to queries from within the GitHub CI.
    # The commented command below is how to get the frame from DCC if you
    # wanted to verify they are the same.

    #curl -O -L --show-error --silent https://dcc.ligo.org/public/0146/P1700349/001/${file}

    curl -O -L --show-error --silent \
         https://github.com/gwastro/pycbc_data/releases/download/ci-data/${file}
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
        cats.ini \
        data.ini \
    --workflow-name pygrb-example \
    --output-dir output \
    --config-overrides \
        results_page:output-path:$(pwd)/html \
    --plan-now

# Start the workflow.

#cp *.gwf output
#cd output
#./start

# Wait for completion and check that everything completed well.

#python ../check_workflow.py