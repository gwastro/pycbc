#!/bin/bash -e

CONFIG_URL=https://github.com/gwastro/pycbc-config/raw/master/test/multi_inspiral
BANK_FILE=gw170817_single_template.hdf
BANK_VETO_FILE=bank_veto_bank.xml
H1_FRAME=https://www.gw-openscience.org/eventapi/html/GWTC-1-confident/GW170817/v3/H-H1_GWOSC_4KHZ_R1-1187006835-4096.gwf
H1_CHANNEL=GWOSC-4KHZ_R1_STRAIN
L1_FRAME=https://dcc.ligo.org/public/0144/T1700406/003/L-L1_CLEANED_HOFT_C02_T1700406_v3-1187008667-4096.gwf
L1_CHANNEL=DCH-CLEAN_STRAIN_C02_T1700406_v3
V1_FRAME=https://www.gw-openscience.org/eventapi/html/GWTC-1-confident/GW170817/v3/V-V1_GWOSC_4KHZ_R1-1187006835-4096.gwf
V1_CHANNEL=GWOSC-4KHZ_R1_STRAIN

echo -e "\\n\\n>> [`date`] Getting template bank"
wget -nv -nc ${CONFIG_URL}/${BANK_FILE}
echo -e "\\n\\n>> [`date`] Bank veto bank"
wget -nv -nc ${CONFIG_URL}/${BANK_VETO_FILE}
for IFO in H1 L1 V1; do
    echo -e "\\n\\n>> [`date`] Getting ${IFO} frame"
    FRAME=${IFO}_FRAME
    wget -nv -nc ${!FRAME}
done

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
    --processing-scheme mkl \
    --instruments H1 L1 V1 \
    --trigger-time ${EVENT} \
    --gps-start-time ${GPS_START} \
    --gps-end-time ${GPS_END} \
    --trig-start-time ${TRIG_START} \
    --trig-end-time ${TRIG_END} \
    --ra 3.44527994344 \
    --dec -0.408407044967 \
    --bank-file ${BANK_FILE} \
    --approximant IMRPhenomD \
    --order -1 \
    --low-frequency-cutoff 30 \
    --snr-threshold 3.0 \
    --chisq-bins "0.9*get_freq('fSEOBNRv4Peak',params.mass1,params.mass2,params.spin1z,params.spin2z)**(2./3.)" \
    --pad-data 8 \
    --strain-high-pass 25 \
    --sample-rate 4096 \
    --channel-name H1:${H1_CHANNEL} L1:${L1_CHANNEL} V1:${V1_CHANNEL} \
    --frame-files \
        H1:`basename ${H1_FRAME}` \
        L1:`basename ${L1_FRAME}` \
        V1:`basename ${V1_FRAME}` \
    --cluster-method window \
    --cluster-window 0.1 \
    --segment-length 256 \
    --segment-start-pad ${START_PAD} \
    --segment-end-pad ${END_PAD} \
    --psd-estimation median \
    --psd-segment-length 32 \
    --psd-segment-stride 8 \
    --psd-num-segments 29 \
    --num-slides 1 \
    --slide-shift 1 \
    --output ${OUTPUT}

echo -e "\\n\\n>> [`date`] Checking output files"
python check_gw170817_trigs.py
