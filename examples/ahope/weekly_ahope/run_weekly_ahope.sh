GPS_START_TIME=961585543
GPS_END_TIME=961671943
export LOGPATH=/usr1/${USER}/log
mkdir -p $LOGPATH
python weekly_ahope.py --start-time ${GPS_START_TIME} --end-time ${GPS_END_TIME} --config-files weekly_ahope.ini,pipedown.ini,inj.ini --pipedown-log-dir ${LOGPATH}
