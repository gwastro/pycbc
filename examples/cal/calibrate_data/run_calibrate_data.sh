#! /bin/bash

# test analysis
GPS_START_TIME=1102660100
GPS_END_TIME=$((${GPS_START_TIME} + 1000))
IFO=L
FRAME_TYPE=L1_R
FRAME_CACHE=frame_caches/L1-DATAFIND-${GPS_START_TIME}-$((${GPS_END_TIME}-${GPS_START_TIME})).lcf
FILTER_FILE=filter_files/L1OAF_8379.txt

# run the datafind script to create a frame cache that contains
# the path of the frame files we requested
ligo_data_find --lal-cache --observatory ${IFO} --type ${FRAME_TYPE} --gps-start-time ${GPS_START_TIME} --gps-end-time ${GPS_END_TIME} --url-type file > ${FRAME_CACHE}

# run executable that will calibrate data in the time domain with a foton filter file
./pycbc_calibrate_data --gps-start-time ${GPS_START_TIME} --gps-end-time ${GPS_END_TIME} --frame-cache ${FRAME_CACHE} --filter-file ${FILTER_FILE}
