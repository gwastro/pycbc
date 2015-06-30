#!/bin/bash

if [ ! -e ../../bin/pycbc_inspiral ]
then
	echo "Please run this script from the tools/static directory of your source installation"
	exit 1
fi

pyinstaller ../../bin/pycbc_inspiral \
--additional-hooks-dir ./hooks/ \
--hidden-import scipy.integrate \
--runtime-hook runtime-scipy.py \
--name pycbc_inspiral_static \
--strip \
--onefile


if [ ! -e dist/pycbc_inspiral_static ]
then
	echo "===== pyinstaller failed ====="
	echo "If you installed pyinstaller from pip, you may need a more recent version"
	echo "which can be obtained via "
	echo "git clone https://github.com/pyinstaller/pyinstaller.git"
fi

