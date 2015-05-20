#!/bin/bash
pyinstaller `which pycbc_inspiral` \
--additional-hooks-dir ./hooks/ \
--runtime-hook runtime-scipy.py \
--name pycbc_inspiral_static \
--strip \
--onefile
