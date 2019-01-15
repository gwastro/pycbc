#!/bin/bash
set -e

# FIXME this is a workaround for a bug in psycopg2 2.6 (required by pegasus)
# see e.g. https://stackoverflow.com/questions/47044854/error-installing-psycopg2-2-6-2
# We should ask pegasus to update their requirement to psycopg2 2.7 which fixes
# this bug
echo -e "Trying to get rid of pg_config"
sudo apt-get -y purge libpq-dev
echo -e "Making sure it is really gone..."
if [ -n "`which pg_config`" ]
then
    echo -e "...still here:"
    which pg_config
    sudo rm -f `which pg_config`
else
    echo -e "...seems gone"
fi

echo -e ">> [`date`] upgrading setuptools and pip"
pip install --upgrade setuptools pip

echo -e ">> [`date`] installing requirements"
pip install -r requirements.txt

echo -e ">> [`date`] installing pycbc"
pip install .

# LAL extra data files
# FIXME, should be a way to make reduced package (with subset of data files)
git clone -n git clone -n https://git.ligo.org/lscsoft/lalsuite-extra.git
git lfs fetch -I data/lalsimulation/SEOBNRv2ROM_*.dat
git lfs fetch -I data/lalsimulation/*ChirpTime*.dat
git lfs fetch -I data/lalsimulation/SEOBNRv4ROM_v2.0.hdf5
mv data/lalsimulation/* ./
