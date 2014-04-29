GPS_START_TIME=967593543
GPS_END_TIME=967679943
export LOGPATH=/usr1/${USER}/log
export PIPEDOWNLOG=/usr1/${USER}
export HTMLDIR=/home/${USER}/public_html/ahope/development/weekly_ahope/test

# Create the log directory
mkdir -p $LOGPATH

# Create the workflow
python weekly_ahope.py --config-files weekly_ahope.ini pipedown.ini inj_pycbc.ini \
--config-overrides ahope:start-time:${GPS_START_TIME} \
ahope:end-time:${GPS_END_TIME} \
ahope:ahope-html-basedir:${HTMLDIR} \
ahope:pipedown-tmp-space:${PIPEDOWNLOG} \
ahope:pipedown-log-path:${LOGPATH}

# Move some files needed to do pegasus planning to the workspace folder
#cp plan.sh ${GPS_START_TIME}-${GPS_END_TIME}/
#cp pegasus.conf ${GPS_START_TIME}-${GPS_END_TIME}/

#echo 'cat <<END_OF_TEXT' >  temp.sh
#cat "site-local.xml"                 >> temp.sh
#echo 'END_OF_TEXT'       >> temp.sh
#bash temp.sh > "${GPS_START_TIME}-${GPS_END_TIME}/site-local.xml"

# Plan the workflow
#cd ${GPS_START_TIME}-${GPS_END_TIME}/
#sh plan.sh weekly_ahope.dax
