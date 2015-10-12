#!/bin/bash

if [ ! -e ../../bin/pycbc_inspiral ]
then
	echo "Please run this script from the tools/static directory of your source installation"
	exit 1
fi

#for prog in `find ../../bin -type f | egrep -v '(gstlal|mvsc|hwinj)'`
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
			# don't rebuild if the executable is already
			# present
			if [ ! -e dist/${exename} ]
			then
				# remove any cached information
				rm -rf ${exename}.spec build/${exename}

				if grep -q ^${exename}$ needs_full_build 
				then
					# This will be used by hooks/hook-pycbc.py to
					# determine what needs to be included
					export NOW_BUILDING=${exename}

					pyinstaller ${prog} \
					  --additional-hooks-dir ./hooks/ \
					  --runtime-hook runtime-scipy.py \
					  --name ${exename} \
					  --strip \
					  --onefile
				else
					pyinstaller ${prog} --name ${exename} --strip --onefile
				fi

				if ! dist/${exename} --help
				then
					echo "Build of ${exename} failed"
					exit 1
				fi
			fi
		fi
	fi
done

