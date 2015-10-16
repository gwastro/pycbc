#!/bin/bash

if [ ! -e ../../bin/pycbc_inspiral ]
then
	echo "Please run this script from the tools/static directory of your source installation"
	exit 1
fi

cat > build_one.sub <<EOT
universe = vanilla
priority = 100
request_memory = 32000
executable = ${PWD}/build_one.sh
arguments = \$(prog)
output = pyinstaller_build.\$(cluster).\$(process).out
error = pyinstaller_build.\$(cluster).\$(process).err
log = /usr1/${USER}/log/pyinstaller_build_1.log
getenv = True
accounting_group = ligo.dev.o1.cbc.nsbh.pycbcoffline 
queue
EOT


for prog in `find ../../bin -type f `
do
	# don't try to pyinstall shell scripts
	if `head -1 ${prog} | grep -q python`
	then
		exename=`basename ${prog}`

		# Some programs can't be compiled.  At the moment these are
		# not needed in any workflow, but this list may need to be
		# revisited
		if ! grep -q ${exename} cant_be_built
		then
			NEW_UUID=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 32 | head -n 1)

			echo "JOB ${NEW_UUID} ${PWD}/build_one.sub DIR ${PWD}"
			echo "VARS ${NEW_UUID} prog=\"${prog}\""
			echo
		fi
	fi
done > build_static.dag


