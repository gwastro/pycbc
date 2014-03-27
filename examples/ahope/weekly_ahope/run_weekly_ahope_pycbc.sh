GPS_START_TIME=961585543
GPS_END_TIME=961671944
export LOGPATH=/usr1/${USER}/log
export PIPEDOWNLOG=/usr1/${USER}
export HTMLDIR=/home/${USER}/public_html/ahope/development/weekly_ahope/test
mkdir -p $LOGPATH
python weekly_ahope.py --config-files weekly_ahope_pycbc.ini pipedown.ini inj_pycbc.ini --config-overrides ahope:start-time:${GPS_START_TIME} ahope:end-time:${GPS_END_TIME} ahope:ahope-html-basedir:${HTMLDIR} ahope:pipedown-tmp-space:${PIPEDOWNLOG} ahope:pipedown-log-path:${LOGPATH}
