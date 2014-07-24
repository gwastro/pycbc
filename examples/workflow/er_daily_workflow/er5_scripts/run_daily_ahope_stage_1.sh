#! /bin/bash

# get yesterday at midnight in GPS time
YESTERDAY=`date -u --date="yesterday" +"%x"`
GPSSTART=`lalapps_tconvert $YESTERDAY`

MONTHDIR=`lalapps_tconvert -f "%Y%m" ${GPSSTART}`
DAYDIR=`lalapps_tconvert -f "%Y%m%d" ${GPSSTART}`

# source code
source /home/cbc/lalsuite/master_20140421_0bf46b5422f722c23a9b4b0b6e71dbae075d46a3/etc/lscsoftrc
source /home/cbc/pycbc/master_20140421_fd7360c846a3040a68b73e516eac98047ebc8c92/etc/pycbc-user-env.sh

#export X509_USER_PROXY='/home/cbc/daily_ihope_er5/grid-proxy-file_3'
export X509_USER_PROXY='/home/cbiwer/projects/daily_ihope_er5/grid-proxy-file_3'
export LIGO_DATAFIND_SERVER='10.13.20.73:80'
export S6_SEGMENT_SERVER='https://segdb-er.ligo.caltech.edu'

# run daily ahope
cd /home/cbc/daily_ahope_er5
python daily_ahope.py  --config-files daily_ahope.ini --start-time ${GPSSTART} -d ${PWD}
cd ${MONTHDIR}/${DAYDIR}
python ../../add_daily_ihope_classad.py
condor_submit_dag daily_ahope.dag
