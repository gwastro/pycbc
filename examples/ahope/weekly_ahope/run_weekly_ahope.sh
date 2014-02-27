GPS_START_TIME=961585543
GPS_END_TIME=961671943
export LOGPATH=/usr1/${USER}/log
mkdir -p $LOGPATH
python weekly_ahope.py --config-files weekly_ahope.ini pipedown.ini inj.ini --config-overrides ahope:start-time:${GPS_START_TIME} ahope:end-time:${GPS_END_TIME} --pipedown-log-dir ${LOGPATH}
