#!/bin/bash
# Set up the site-local with the correct paths

echo 'cat <<END_OF_TEXT' >  temp.sh
cat "site-local.xml"                 >> temp.sh
echo 'END_OF_TEXT'       >> temp.sh
bash temp.sh > "${GPS_START_TIME}-${GPS_END_TIME}/site-local.xml"

# Plan the workflow
pegasus-plan --conf pegasus.conf -d $1 --sites local -o local --dir ${LOGPATH}
