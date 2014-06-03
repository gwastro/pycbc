#! /bin/bash

# get yesterday at midnight in GPS time
YESTERDAY=`date -u --date="yesterday" +"%x"`
GPSSTART=`lalapps_tconvert $YESTERDAY`
GPSEND=$(($GPSSTART + 60 * 60 * 24 + 72 + 8))

MONTHDIR=`lalapps_tconvert -f "%Y%m" ${GPSSTART}`
DAYDIR=`lalapps_tconvert -f "%Y%m%d" ${GPSSTART}`

IFO=L1

# make dir
TMPDIR=/home/cbc/public_html/daily_cbc_offline/${MONTHDIR}/${DAYDIR}
mkdir -p ${TMPDIR}

# create INSPIRAL cache
cd /home/cbc/daily_ahope_er5
ls ${MONTHDIR}/${DAYDIR}/${IFO}-INSPIRAL_30MILLISEC_CLUSTERED-*.xml.gz | lalapps_path2cache -o ${TMPDIR}/${IFO}-INSPIRAL_30MILLISEC_CLUSTERED.cache

# create TMPLTBANK cache
ls ${MONTHDIR}/${DAYDIR}/${IFO}-TMPLTBANK-*.xml.gz | lalapps_path2cache -o ${TMPDIR}/${IFO}-TMPLTBANK_30MILLISEC_CLUSTERED.cache

# set-up env for running gw_summary
export PYTHONPATH=$PYTHONPATH:/home/cbiwer/local/gwpy_dependancies/lib/python2.6/site-packages:/home/cbiwer/local/gwpy_dependancies/lib64/python2.6/site-packages
export PYTHONPATH=$PYTHONPATH:/home/cbiwer/local/gwpy_master/lib/python2.6/site-packages:/home/cbiwer/local/gwpy_master/lib64/python2.6/site-packages
export PYTHONPATH=$PYTHONPATH:/home/cbiwer/local/gwsumm_master/lib/python2.6/site-packages:/home/cbiwer/local/gwsumm_master/lib64/python2.6/site-packages
export PATH=$PATH:/home/cbiwer/local/gwpy_master/bin
export PATH=$PATH:/home/cbiwer/local/gwsumm_master/bin

export X509_USER_PROXY='/home/cbiwer/projects/daily_ihope_er5/grid-proxy-file_3'
export LIGO_DATAFIND_SERVER='10.13.20.73:80'
export S6_SEGMENT_SERVER='https://segdb-er.ligo.caltech.edu'

# run gw_summary
cd /home/cbc/daily_ahope_er5/
kinit christopher.biwer@LIGO.ORG -k -t cbc.keytab
gw_summary --verbose --ifo ${IFO} --output-dir /home/cbc/public_html/daily_cbc_offline/ahope_detchar --config-file defaults.ini --config-file ihope.ini --max-processes 12 --day ${DAYDIR}

