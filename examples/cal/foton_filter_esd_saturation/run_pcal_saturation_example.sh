#! /bin/bash

# Copyright (C) 2015  Christopher M. Biwer
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

set -e

# frame options note these options are for reading SWSTAT and GAIN channels
# first option is IFO, eg. H1 or L1
IFO=${1}
FRAME_TYPE=${IFO}_R

# get single-column ASCII file that contains h(t) time series
# second option on command line should be the *full* path to the data file
# eg. /home/user/path/to/file.txt
DATA_FILE=${2}

# output options
HTMLDIR=/home/${USER}/public_html/hwinj_check/doc/
mkdir -p ${HTMLDIR}

# checkout the calibration SVN to get foton filter files
svn co --username=christopher.biwer https://daqsvn.ligo-la.caltech.edu/svn/h1_filter_files/h1_archive
svn co --username=christopher.biwer https://daqsvn.ligo-la.caltech.edu/svn/l1_filter_files/l1_archive

# foton filter files
CALEX_FILTER_FILE=${IFO,,}_archive/${IFO}CALEX.txt

# look-back time options
# these are the time options to check the SWSTAT and GAIN channels
# note you must wait this much time before you can use any changes made to the
# digital controls systems (ie. filterbanks)
GPS_START_TIME=$((`lalapps_tconvert` - 900))
GPS_END_TIME=$((${GPS_START_TIME}+1))

# filter injection
sh pycbc_check_pcal_saturation.sh --data-file ${DATA_FILE} --gps-start-time ${GPS_START_TIME} --gps-end-time ${GPS_END_TIME} --ifo ${IFO} --frame-type ${FRAME_TYPE} --calex-filter-file ${PWD}/${CALEX_FILTER_FILE}

# plot strain
INPUT_FILE=${DATA_FILE}
TIMESERIES_FILE=${HTMLDIR}/${IFO}-TIMESERIES_HOFT.png
pycbc_plot_hwinj --input-file ${INPUT_FILE} --output-file ${TIMESERIES_FILE} \
    --title "HOFT" --y-label "Strain" --verbose

# plot filtered data
INPUT_FILE=`ls pcal_output/${IFO}-FILTER_PINJX_TRANSIENT-${GPS_START_TIME}.txt`
TIMESERIES_FILE=${HTMLDIR}/${IFO}-TIMESERIES_PINJX_TRANSIENT.png
pycbc_plot_hwinj --input-file ${INPUT_FILE} --output-file ${TIMESERIES_FILE} \
    --title "FILTERED HOFT" --y-label "Counts" --verbose
