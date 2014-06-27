#!/bin/bash
# Set up the site-local with the correct paths

echo 'cat <<END_OF_TEXT' >  temp.sh
cat "site-local.xml"                 >> temp.sh
echo 'END_OF_TEXT'       >> temp.sh
bash temp.sh > site-local-parsed.xml

# Plan the workflow
echo "Generating concrete workflow"
pegasus-plan --conf pegasus.conf -d $1 --sites local -o local --dir $2 --cleanup none

#echo "Generating the workflow shell script"
#pegasus-plan --conf pegasus.conf -d $1 --sites local -o local --dir $PWD --cleanup none
