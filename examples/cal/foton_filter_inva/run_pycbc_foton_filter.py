#! /bin/bash -e

# time options
GEOCENT_END_TIME=1126257416
GPS_START_TIME=1126256416
GPS_END_TIME=$((${GPS_START_TIME} + 2048))

# frame options note these options are for reading h(t) channels
IFO=H1
FRAME_TYPE=H1_HOFT_C00
CHANNEL_NAME=H1:GDS-CALIB_STRAIN

# injection options
SAMPLE_RATE=16384

# foton options
FILTER_FILE=/home/${USER}/src/calsvn/h1_archive/H1CALCS.txt
MODEL_NAME=CAL
FILTERBANK_NAME=INJ_BLIND

# output options
OUTPUT_FILE=output.txt

# where to run the executables
RUN_DIR=./test_case

# change into directory where executables will be run
mkdir -p ${RUN_DIR}
cd ${RUN_DIR}

# create an injection
pycbc_generate_hwinj --fftw-measure-level 0 --geocentric-end-time ${GEOCENT_END_TIME} --gps-start-time ${GPS_START_TIME} --gps-end-time ${GPS_END_TIME} --frame-type ${FRAME_TYPE} --channel-name ${CHANNEL_NAME} --approximant SpinTaylorT4 --order pseudoFourPN --mass1 1.4 --mass2 1.4 --inclination 1.04719755 --polarization 0.0 --ra 0.0 --dec 0.0 --taper TAPER_STARTEND --network-snr 28 --low-frequency-cutoff 40.0 --h1 --sample-rate 16384 --pad-data 8 --strain-high-pass 30.0 --psd-estimation median --psd-segment-length 16 --psd-segment-stride 8

# get injection single-column ASCII file
DATA_FILE=`ls ${IFO}-*.txt`

# frame options note these options are for reading SWSTAT channels
FRAME_TYPE=H1_R

# filter data and injection
python ../pycbc_foton_filter --model-name ${MODEL_NAME} --gps-start-time ${GPS_START_TIME} --gps-end-time ${GPS_END_TIME} --filterbank-name ${FILTERBANK_NAME} --data-file ${DATA_FILE} --filter-file ${FILTER_FILE} --output-file ${OUTPUT_FILE} --sample-rate ${SAMPLE_RATE} --frame-type ${FRAME_TYPE} --ifo ${IFO}

