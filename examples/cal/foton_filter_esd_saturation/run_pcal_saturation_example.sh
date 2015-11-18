#! /bin/bash -e

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

# frame options note these options are for reading SWSTAT and GAIN channels
IFO=H1
FRAME_TYPE=H1_R

# get single-column ASCII file that contains h(t) time series
DATA_FILE=`ls ${PWD}/H1-HWINJ_CBC-*-*.txt`

# injection options
SAMPLE_RATE=16384

# output options
HTMLDIR=/home/${USER}/public_html/hwinj_check/doc/
mkdir -p ${HTMLDIR}

# checkout the calibration SVN to get foton filter files
svn co --username=christopher.biwer https://daqsvn.ligo-la.caltech.edu/svn/h1_filter_files/h1_archive
svn co --username=christopher.biwer https://daqsvn.ligo-la.caltech.edu/svn/l1_filter_files/l1_archive

# foton filter files
CALEX_FILTER_FILE=h1_archive/H1CALEX.txt

# time options
GPS_START_TIME=$((`lalapps_tconvert` - 30000))
GPS_END_TIME=$((${GPS_START_TIME}+1))

# filter injection
sh pycbc_check_pcal_saturation.sh --data-file ${DATA_FILE} --gps-start-time ${GPS_START_TIME} --gps-end-time ${GPS_END_TIME} --ifo ${IFO} --frame-type ${FRAME_TYPE} --calex-filter-file ${PWD}/${CALEX_FILTER_FILE}

# source detchar env
# note that this path is specific to the CIT, LHO, and LLO clusters
source /home/detchar/opt/gwpysoft/etc/gwpy-user-env.sh

# output options
INPUT_FILE=${DATA_FILE}
TIMESERIES_FILE=${HTMLDIR}/${IFO}-TIMESERIES_HOFT.png
SPECTROGRAM_FILE=${HTMLDIR}/${IFO}-SPECTROGRAM_HOFT.png

# run plotting script
python gwpy_plot_hwinj ${INPUT_FILE} ${TIMESERIES_FILE} ${SPECTROGRAM_FILE} 0 1027 -1e-22 1e-22 1e-21 1e-27

# output options
INPUT_FILE=`ls pcal_output/${IFO}-FILTER_PINJX_HARDWARE-*.txt`
TIMESERIES_FILE=${HTMLDIR}/${IFO}-TIMESERIES_PINJX_HARDWARE.png
SPECTROGRAM_FILE=${HTMLDIR}/${IFO}-SPECTROGRAM_PINJX_HARDWARE.png

# run plotting script
python gwpy_plot_hwinj ${INPUT_FILE} ${TIMESERIES_FILE} ${SPECTROGRAM_FILE} 0 21 -30000 30000 1e1 1e-7
