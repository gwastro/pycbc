#! /bin/bash

# test analysis
GPS_START_TIME=1102660100
GPS_END_TIME=$((${GPS_START_TIME} + 1000))
IFO=L
FRAME_TYPE=L1_R
FRAME_CACHE=../calibrate_data/frame_caches/L1-DATAFIND-${GPS_START_TIME}-$((${GPS_END_TIME}-${GPS_START_TIME})).lcf
FILTER_FILE=../calibrate_data/filter_files/L1OAF_8379.txt

CHANNEL_NAME=L1:LSC-DARM_IN1_DQ
FILTER_NAME=CAL_SUM_DARM_ERR
SWSTAT_NAME=L1:OAF-CAL_SUM_DARM_ERR_SWSTAT

# run the datafind script to create a frame cache that contains
# the path of the frame files we requested
ligo_data_find --lal-cache --observatory ${IFO} --type ${FRAME_TYPE} --gps-start-time ${GPS_START_TIME} --gps-end-time ${GPS_END_TIME} --url-type file > ${FRAME_CACHE}

# run executable that will calibrate data in the time domain with a foton filter file
./pycbc_filter_data --gps-start-time ${GPS_START_TIME} --gps-end-time ${GPS_END_TIME} --frame-cache ${FRAME_CACHE} --filter-file ${FILTER_FILE} --channel-name ${CHANNEL_NAME} --swstat-name ${SWSTAT_NAME} --filter-name ${FILTER_NAME}
