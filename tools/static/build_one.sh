#!/bin/bash

source ../../../../bin/activate
export PYINSTALLER_CONFIG_DIR=`mktemp --tmpdir -d -t pyinstaller-XXXXXXXXXX`

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

dist/${exename} --help >&  ${exename}.out

if [ $? != 0 ]
then
	echo "Build of ${exename} failed"
        rm -rf ${PYINSTALLER_CONFIG_DIR}
	exit 1
fi

if grep -q 'No module named cython_blas' ${exename}.out
then
	echo "Build of ${exename} is missing hidden dependancies"
        rm -rf ${PYINSTALLER_CONFIG_DIR}
	exit 1
fi

rm -rf ${PYINSTALLER_CONFIG_DIR}
