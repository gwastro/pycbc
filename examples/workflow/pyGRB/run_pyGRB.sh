#!/bin/sh
export LOGPATH=/usr1/${USER}/log
export RUN_DIR=${PWD}
GPS_START_TIME=969645000
GPS_END_TIME=969685000
LOCAL_CONFIG_FILES="pyGRB.ini"
BANK_FILE=${PWD}/TMPLTBANKS/H1-TMPLTBANK_GRB100928A_DATAFIND-969673046-4992.xml
#GPS_START_TIME=969670108
#GPS_END_TIME=969681108
./pyGRB.py \
--local-config-files ${LOCAL_CONFIG_FILES} \
--config-overrides \
workflow:start-time:${GPS_START_TIME} \
workflow:end-time:${GPS_END_TIME} \
workflow-tmpltbank:tmpltbank-pregenerated-bank:${BANK_FILE}
