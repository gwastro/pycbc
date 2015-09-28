#!/bin/sh
########
# THIS EXAMPLE SCRIPT GENERATES AN OFFLINE WORKFLOW FOR GRB 150906B
########
export LOGPATH=${LOCALDIR}/log
mkdir -p ${LOGPATH}
export OUT_DIR=${PWD}
export LAL_SRC=${LOCALDIR}/git/lalsuite
GRB_TIME=1125564162
RA=159.239
DEC=-25.603
SKY_ERROR=0
GRB_NAME=150906B
export HTML_DIR=${HOME}/public_html/LVC/ER8/pygrb/GRB${GRB_NAME}
CONFIG_FILES="er8_main_offline.ini er8_postprocessing.ini er8_injections_offline.ini"
pycbc_make_offline_grb_workflow \
--config-files ${CONFIG_FILES} \
--config-overrides \
workflow:output-directory:${OUT_DIR} \
workflow:ra:${RA} \
workflow:dec:${DEC} \
workflow:sky-error:${SKY_ERROR} \
workflow:trigger-name:${GRB_NAME} \
workflow:trigger-time:${GRB_TIME} \
workflow:start-time:$(( GRB_TIME - 5096 )) \
workflow:end-time:$(( GRB_TIME + 5096 )) \
workflow:html-dir:${HTML_DIR}
