#!/bin/bash

source ../../../../bin/activate
export PYINSTALLER_CONFIG_DIR=/usr1/${USER}/`hostname`/${BASHPID}

prog=$1

exename=`basename ${prog}`

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
