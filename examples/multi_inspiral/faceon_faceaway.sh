#!/bin/bash -e

EVENT=1187008882
STEP=300
START=$((EVENT - 2 * STEP))
END=$((EVENT + STEP + 64))
GPS_START=$((START - 32))
GPS_END=$((END + 32))
DUR=$((GPS_END - GPS_START))
CONFIG_URL=https://github.com/gwastro/pycbc-config/raw/master/test/multi_inspiral
INJ_FILE=faceon_faceaway_injs.hdf

echo -e "\\n\\n>> [`date`] Getting injection file"
wget -nv -nc ${CONFIG_URL}/${INJ_FILE}

# Generate mock data
for IFO in 'H1' 'L1' 'V1'; do
    echo -e "\\n\\n>> [`date`] Generating mock data for ${IFO}"
    pycbc_condition_strain \
        --fake-strain zeroNoise \
        --sample-rate 4096 \
        --gps-start-time ${GPS_START} \
        --gps-end-time ${GPS_END} \
        --channel-name ${IFO}:SIMULATED_GW170817 \
        --injection-file ${INJ_FILE} \
        --output-strain-file ${IFO}-SIMULATED_GW170817-${GPS_START}-${DUR}.gwf
done

TRIG_START=$((START + 256))
TRIG_END=$((END - 32))
PAD=8
BANK_FILE=gw170817_single_template.hdf
BANK_VETO_FILE=bank_veto_bank.xml
CHANNEL=SIMULATED_GW170817
RA=5.016076270234897
DEC=-0.408407044967

echo -e "\\n\\n>> [`date`] Getting template bank"
wget -nv -nc ${CONFIG_URL}/${BANK_FILE}
echo -e "\\n\\n>> [`date`] Bank veto bank"
wget -nv -nc ${CONFIG_URL}/${BANK_VETO_FILE}

for POL in 'standard' 'left' 'right' 'left+right'; do
    echo -e "\\n\\n>> [`date`] Running pycbc_multi_inspiral with ${POL} projection"
    pycbc_multi_inspiral \
        --projection ${POL} \
        --processing-scheme mkl \
        --instruments H1 L1 V1 \
        --trigger-time ${EVENT} \
        --gps-start-time $((GPS_START + PAD)) \
        --gps-end-time $((GPS_END - PAD)) \
        --trig-start-time ${TRIG_START} \
        --trig-end-time ${TRIG_END} \
        --ra ${RA} \
        --dec ${DEC} \
        --bank-file ${BANK_FILE} \
        --approximant IMRPhenomD \
        --order -1 \
        --low-frequency-cutoff 30 \
        --pad-data 8 \
        --strain-high-pass 25 \
        --sample-rate 4096 \
        --channel-name H1:${CHANNEL} L1:${CHANNEL} V1:${CHANNEL} \
        --frame-files \
            H1:H1-${CHANNEL}-${GPS_START}-${DUR}.gwf \
            L1:L1-${CHANNEL}-${GPS_START}-${DUR}.gwf \
            V1:V1-${CHANNEL}-${GPS_START}-${DUR}.gwf \
        --snr-threshold 4.0 \
        --chisq-bins "0.9*get_freq('fSEOBNRv4Peak',params.mass1,params.mass2,params.spin1z,params.spin2z)**(2./3.)" \
        --bank-veto-bank-file ${BANK_VETO_FILE} \
        --cluster-method window \
        --cluster-window 0.1 \
        --segment-length 256 \
        --segment-start-pad 111 \
        --segment-end-pad 17 \
        --psd-model aLIGOEarlyHighSensitivityP1200087 \
        --output ${POL}.hdf
done

echo -e "\\n\\n>> [`date`] Checking output files"
python ./check_faceon_faceaway_trigs.py
