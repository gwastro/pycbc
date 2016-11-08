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

# where to run the executables
RUN_DIR=./pcal_output

# default option values
DATA_FILE=""
GPS_START_TIME=""
GPS_END_TIME=""
IFO=""
FRAME_TYPE=""
SAMPLE_RATE=16384
CALEX_FILTER_FILE=""

# parse command line
GETOPT_CMD=`getopt -o d:c:a:F:h:l --long data-file:,gps-start-time:,gps-end-time:,ifo:,frame-type:,sample-rate:,calex-filter-file:,sus-filter-file:,help -n 'pycbc_check_esd_saturation.sh' -- "$@"`
eval set -- "$GETOPT_CMD"
while true ; do
  case "$1" in
    -d|--data-file)
      case "$2" in
        "") shift 2 ;;
        *) DATA_FILE=$2 ; shift 2 ;;
      esac ;;
    -s|--gps-start-time)
      case "$2" in
        "") shift 2 ;;
        *) GPS_START_TIME=$2 ; shift 2 ;;
      esac ;;
    -e|--gps-end-time)
      case "$2" in
        "") shift 2 ;;
        *) GPS_END_TIME=$2 ; shift 2 ;;
      esac ;;
    -i|--ifo)
      case "$2" in
        "") shift 2 ;;
        *) IFO=$2 ; shift 2 ;;
      esac ;;
    -f|--frame-type)
      case "$2" in
        "") shift 2 ;;
        *) FRAME_TYPE=$2 ; shift 2 ;;
      esac ;;
    -r|--sample-rate)
      case "$2" in
        "") shift 2 ;;
        *) SAMPLE_RATE=$2 ; shift 2 ;;
      esac ;;
    -c|--calex-filter-file)
      case "$2" in
        "") shift 2 ;;
        *) CALEX_FILTER_FILE=$2 ; shift 2 ;;
      esac ;;
    -h|--help)
      echo "usage: pycbc_check_esd_saturation.sh [-h] [--options]"
      echo
      echo "required arguments:"
      echo "  -d, --data-file DATA_FILE                       single-column ASCII file"
      echo "  -s, --gps-start-time GPS_START_TIME             time to start reading data for filterbank values"
      echo "  -e, --gps-end-time GPS_END_TIME                 time to stop reading data for filterbank values"
      echo "  -i, --ifo IFO                                   IFO, eg. H1 or L1"
      echo "  -f, --frame-type FRAME_TYPE                     frame type that has SWSTAT and GAIN channels"
      echo "  -r, --sample-rate SAMPLE_RATE                   sample rate of the input file"
      echo "  -c, --calex-filter-file CALEX_FILTER_FILE       path to CALEX.txt"
      echo
      echo "Filter a single-column ASCII file to see if it saturates the ETMY DAC."
      echo
      echo "You need to have ROOT and foton python packages in your environment."
      echo "To do this on the LLO or LHO clusters do:"
      echo 
      echo "NAME=/path/to/virtualenv"
      echo "cd ${NAME}/lib64/python2.6/site-packages"
      echo "ln -s /usr/lib64/python2.6/site-packages/libPyROOT.so"
      echo "ln -s /usr/lib64/python2.6/site-packages/ROOT.py"
      echo "ln -s /usr/lib64/python2.6/site-packages/ROOTwriter.py"
      echo "cd ${NAME}/lib/python2.6/site-packages"
      echo "ln -s /usr/lib/python2.6/site-packages/foton.py"
      echo
      exit 0 ;;
    --) shift ; break ;;
    *) echo "Internal error!" ; exit 1 ;;
  esac
done

# save current directory to get path to executables
EXE_DIR=${PWD}

# change into directory where executables will be run
mkdir -p ${RUN_DIR}
cd ${RUN_DIR}

# filter with CAL-INJ_HARDWARE
MODEL_NAME=CAL
FILTERBANK_NAME=PINJX_TRANSIENT
#DATA_FILE=`ls /home/cbiwer/src/pycbc_foton_dev/examples/cal/foton_filter_INJ_BLIND/hwinj/${IFO}-HWINJ_CBC-*-*.txt`
OUTPUT_FILE=${IFO}-FILTER_${FILTERBANK_NAME}-${GPS_START_TIME}.txt
python ${EXE_DIR}/pycbc_foton_filter \
   --model-name ${MODEL_NAME} \
   --gps-start-time ${GPS_START_TIME} \
   --gps-end-time ${GPS_END_TIME} \
   --filterbank-name ${FILTERBANK_NAME} \
   --filterbank-ignore-off \
   --filterbank-gain 1.0 \
   --data-file ${DATA_FILE} \
   --filter-file ${CALEX_FILTER_FILE} \
   --output-file ${OUTPUT_FILE} \
   --sample-rate ${SAMPLE_RATE} \
   --frame-type ${FRAME_TYPE} \
   --ifo ${IFO}
