#!/bin/bash

if [ ! -e ../../bin/pycbc_inspiral ]
then
	echo "Please run this script from the tools/static directory of your source installation"
	exit 1
fi

for prog in `find ../../bin -type f`
do
	# don't try to pyinstall shell scripts
	if `head -1 ${prog} | grep -q python`
	then 
		exename=`basename ${prog}`
		target=${exename}_static

		# don't rebuild if the executable is already
		# present
		if [ ! -e dist/${target} ]
		then
			# remove any cached information
			rm -rf ${target}.spec build/${target}

			# This will be used by hooks/hook-pycbc.py to
			# determine what needs to be included
			export NOW_BUILDING=${exename}

			pyinstaller ${prog} \
			  --additional-hooks-dir ./hooks/ \
			  --runtime-hook runtime-scipy.py \
			  --name ${target} \
			  --strip \
			  --onefile
		fi
	fi
done

