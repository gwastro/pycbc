#!/bin/bash
# Set up the site-local with the correct paths

# Check I have two arguments supplied

if [ "x$1" == "x" ]
then
  echo "I take two arguments, the name of the dax file and the logpath. None supplied."
  exit 1
fi

if [ "x$2" == "x" ]
then
  echo "I take two arguments, the name of the dax file and the logpath. Only got one: $1"
  exit 1
fi


echo 'cat <<END_OF_TEXT' >  temp.sh
cat "site-local.xml"                 >> temp.sh
echo 'END_OF_TEXT'       >> temp.sh
bash temp.sh > site-local-parsed.xml

# Plan the workflow
echo "Generating concrete workflow"
pegasus-plan --conf pegasus.conf -d $1 --sites local -o local --dir $2 --nocleanup

#echo "Generating the workflow shell script"
#pegasus-plan --conf pegasus.conf -d $1 --sites local -o local --dir $PWD --nocleanup
