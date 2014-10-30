#! /bin/bash

# Copyright (C) 2014 Christopher M. Biwer
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

"""
A wrapper around pycbc_make_daily_workflow so that single-detector data
can be analyzed automatically. This is designed to be called within
a contab once a day, a couple hours after midnight UTC.
"""

# get yesterday at midnight in GPS time
YESTERDAY=`date -u --date="yesterday" +"%x"`
GPS_START_TIME=`lalapps_tconvert $YESTERDAY`
MONTHDIR=`lalapps_tconvert -f "%Y%m" ${GPS_START_TIME}`
DAYDIR=`lalapps_tconvert -f "%Y%m%d" ${GPS_START_TIME}`

# location of configuration file(s) and dirs
export INSTALLED_CONFIG_FILES='example_daily_pycbc_zpk.ini'
export HTMLDIR='/home/cbc/public_html/daily_cbc_offline'
export ASSETDIR='/home/cbc/lalsuite/lalapps/src/inspiral'

# use CBC account service certificate
export X509_USER_PROXY='/home/cbc/daily_ahope_er5/cert_daily_ahope/daily-ahope.pem'

# source code
source /home/cbc/lalsuite/master/etc/lscsoftrc
source /home/cbc/pycbc/master/etc/pycbc-user-env.sh

# source detchar
source /home/detchar/opt/gwpysoft/etc/gwpy-user-env.sh

# move to run dir
RUNDIR=/home/cbc/daily_ahope_er5
cd ${RUNDIR}

# run daily_ahope
pycbc_make_daily_workflow --installed-config-files ${INSTALLED_CONFIG_FILES} \
                          --start-time ${GPS_START_TIME} --output-dir ${RUNDIR} \
                          --config-overrides workflow:workflow-html-basedir:${HTMLDIR} \
                                             workflow:workflow-asset-dir:${ASSETDIR}

# plan the workflow
cd ${RUNDIR}/${MONTHDIR}/${DAYDIR}
pycbc_basic_pegasus_plan daily_ahope.dax ${RUNDIR}

