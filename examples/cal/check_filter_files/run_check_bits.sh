#! /bin/bash

# test analysis
GPS_START_TIME=1102660100
GPS_END_TIME=$((${GPS_START_TIME} + 2048))
IFO=L
FRAME_TYPE=L1_R
FRAME_CACHE=../calibrate_data/frame_caches/L1-DATAFIND-${GPS_START_TIME}-$((${GPS_END_TIME}-${GPS_START_TIME})).lcf
FILTER_FILE=../calibrate_data/filter_files/L1OAF_8379.txt

# run executable that will calibrate data in the time domain with a foton filter file
./pycbc_check_bits --filter-file ${FILTER_FILE} --filterbank-name CAL_SUM_DARM_ERR
