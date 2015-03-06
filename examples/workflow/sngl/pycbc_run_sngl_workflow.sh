#! /bin/bash

#source ../../../etc/lscsoftrc
#source ../../../etc/pycbc-user-env.sh

export GPS_START_TIME=1102089616
export GPS_END_TIME=1102100616
export WORKFLOW_NAME=${GPS_START_TIME}-${GPS_END_TIME}

export HTMLDIR=/home/${USER}/public_html/daily_cbc_offline/${WORKFLOW_NAME}
export OUTPUT=/~${USER}/daily_cbc_offline/${WORKFLOW_NAME}

export INSTALLED_CONFIG_FILES="sngl/example_sngl.ini"

export LOGPATH=/usr1/${USER}/log
mkdir -p $LOGPATH

python /home/cbiwer/src/pycbc/bin/pycbc_make_sngl_workflow --name ${WORKFLOW_NAME} \
                          --installed-config-files ${INSTALLED_CONFIG_FILES} \
                          --config-overrides workflow:start-time:${GPS_START_TIME} \
                                             workflow:end-time:${GPS_END_TIME} \
                                             workflow:workflow-html-basedir:${HTMLDIR} \
                                             workflow:workflow-output-server:${OUTPUT} \
                                             html:analysis-subtitle:${GPS_START_TIME}-${GPS_END_TIME}

cd ${WORKFLOW_NAME}
pycbc_basic_pegasus_plan ${WORKFLOW_NAME}.dax ${LOGPATH}
