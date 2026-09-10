#!/bin/bash

set -e

CONFIG_URL=https://github.com/gwastro/pycbc-config/raw/master/test/multi_inspiral
BANK_FILE=gw170817_single_template.hdf
echo -e "\\n\\n>> [`date`] Getting template bank"
wget -nv -nc ${CONFIG_URL}/${BANK_FILE}

EVENT=1187008882
PAD=8
START_PAD=111
END_PAD=17
GPS_START=$((EVENT - 192 - PAD))
GPS_END=$((EVENT + 192 + PAD))
TRIG_START=$((GPS_START + START_PAD))
TRIG_END=$((GPS_END - END_PAD))
OUTPUT=GW170817_test_output.hdf

echo -e "\\n\\n>> [`date`] Running pycbc_multi_inspiral on GW170817 data"
pycbc_multi_inspiral \
    --verbose \
    --projection left+right \
    --instruments H1 \
    --channel-name H1:GWOSC-16KHZ_R1_STRAIN \
    --frame-type H1:GWOSC \
    --trigger-time ${EVENT} \
    --gps-start-time ${GPS_START} \
    --gps-end-time ${GPS_END} \
    --trig-start-time ${TRIG_START} \
    --trig-end-time ${TRIG_END} \
    --ra '3.44527994344 rad' \
    --dec '-0.408407044967 rad' \
    --bank-file ${BANK_FILE} \
    --approximant IMRPhenomD \
    --order -1 \
    --low-frequency-cutoff 30 \
    --sngl-snr-threshold 3.0 \
    --chisq-bins "0.9*get_freq('fSEOBNRv4Peak',params.mass1,params.mass2,params.spin1z,params.spin2z)**(2./3.)" \
    --pad-data 8 \
    --strain-high-pass 25 \
    --sample-rate 4096 \
    --autogating-threshold 100 \
    --autogating-cluster 0.5 \
    --autogating-width 0.25 \
    --autogating-taper 0.25 \
    --autogating-pad 0 \
    --cluster-method window \
    --cluster-window 0.1 \
    --segment-length 256 \
    --segment-start-pad ${START_PAD} \
    --segment-end-pad ${END_PAD} \
    --psd-estimation median \
    --psd-segment-length 32 \
    --psd-segment-stride 8 \
    --psd-num-segments 29 \
    --do-shortslides \
    --slide-shift 1 \
    --output ${OUTPUT}

echo -e "\\n\\n>> [`date`] Checking output files"
python check_gw170817_trigs.py
