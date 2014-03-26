GPS_START_TIME=961585543
GPS_END_TIME=961671944
export LOGPATH=/usr1/${USER}/log
mkdir -p $LOGPATH
python weekly_ahope.py --config-files weekly_ahope_pycbc.ini pipedown.ini inj_pycbc.ini --config-overrides ahope:start-time:${GPS_START_TIME} ahope:end-time:${GPS_END_TIME} --pipedown-log-dir ${LOGPATH}
