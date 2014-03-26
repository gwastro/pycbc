#!/bin/bash
GPS_START_TIME=961585543
GPS_END_TIME=961671943
export LOGPATH=/usr1/${USER}/log
mkdir -p $LOGPATH
python weekly_ahope.py --config-files weekly_ahope.ini pipedown.ini inj.ini \
                       --config-overrides ahope:start-time:${GPS_START_TIME} ahope:end-time:${GPS_END_TIME} \
                       --pipedown-log-dir ${LOGPATH}

# Move some files needed to do pegasus planning to the workspace folder
cp plan.sh ${GPS_START_TIME}-${GPS_END_TIME}/
cp pegasus.conf ${GPS_START_TIME}-${GPS_END_TIME}/

echo 'cat <<END_OF_TEXT' >  temp.sh
cat "site-local.xml"                 >> temp.sh
echo 'END_OF_TEXT'       >> temp.sh
bash temp.sh > "${GPS_START_TIME}-${GPS_END_TIME}/site-local.xml"

# Plan the workflow
cd ${GPS_START_TIME}-${GPS_END_TIME}/
sh plan.sh
