#!/bin/sh
export LOGPATH=/usr1/${USER}/log
mkdir -p ${LOGPATH}
export RUN_DIR=${PWD}
RA=0
DEC=0
GRB_NAME=
GRB_TIME=969650841
LOCAL_CONFIG_FILES="pygrb.ini"
BANK_FILE=${PWD}/TMPLTBANKS/H1-TMPLTBANK_GRB100928A_DATAFIND-969673046-4992.xml
./pygrb.py \
--local-config-files ${LOCAL_CONFIG_FILES} \
--config-overrides \
workflow:ra:${RA} \
workflow:dec:${DEC} \
workflow:trigger-name:${GRB_NAME} \
workflow:trigger-time:${GRB_TIME} \
workflow:start-time:$(( GRB_TIME - 4096 )) \
workflow:end-time:$(( GRB_TIME + 4096 )) \
workflow-tmpltbank:tmpltbank-pregenerated-bank:${BANK_FILE}
