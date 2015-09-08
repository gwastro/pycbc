#!/bin/sh
export LOGPATH=${LOCALDIR}/log
mkdir -p ${LOGPATH}
export OUT_DIR=${PWD}
export LAL_SRC=${LOCALDIR}/git/lalsuite
RA=72.440
DEC=-21.590
SKY_ERROR=4.450
GRB_NAME=150604
GRB_TIME=1117448699
export HTML_DIR=${HOME}/public_html/pygrb/GRB${GRB_NAME}
LOCAL_CONFIG_FILES="er7_offline_main.ini er7_postprocessing.ini er7_injections.ini"
pygrb_make_offline_workflow \
--local-config-files ${LOCAL_CONFIG_FILES} \
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
