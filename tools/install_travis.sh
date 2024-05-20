#!/bin/bash
set -e

########################Workarounds and Hacks #################################
# FIXME this is a workaround for a bug in psycopg2 2.6 (required by pegasus)
# see e.g. https://stackoverflow.com/questions/47044854/error-installing-psycopg2-2-6-2
# We should ask pegasus to update their requirement to psycopg2 2.7 which fixes
# this bug
echo -e "Trying to get rid of pg_config"
sudo apt-get -o Acquire::Retries=3 -y purge libpq-dev
echo -e "Making sure it is really gone..."
if [ -n "`which pg_config`" ]
then
    echo -e "...still here:"
    which pg_config
    sudo rm -f `which pg_config`
else
    echo -e "...seems gone"
fi

# Fixme replace with apt-get from intel repo
# get library needed to build documentation
wget_opts="-c --passive-ftp --no-check-certificate --tries=5 --timeout=30"
primary_url="https://www.atlas.aei.uni-hannover.de/~dbrown/cea5bd67440f6c3195c555a388def3cc6d695a5c/x86_64/composer_xe_2015.0.090"
p="libmkl_rt.so"
test -r $p || wget $wget_opts ${primary_url}/${p}
chmod +x $p

# LAL extra data files
# FIXME, should be a way to make reduced package (with subset of data files)
GIT_CLONE_PROTECTION_ACTIVE=false GIT_LFS_SKIP_SMUDGE=1 git clone https://git.ligo.org/lscsoft/lalsuite-extra
cd lalsuite-extra
git lfs pull -I "data/lalsimulation/SEOBNRv2ROM_*.dat"
git lfs pull -I "data/lalsimulation/*ChirpTime*.dat"
git lfs pull -I "data/lalsimulation/SEOBNRv4ROM_v2.0.hdf5"
mv data/lalsimulation/* ../
cd ../

###############################################################################

echo -e ">> [`date`] upgrading setuptools and pip"
pip install --upgrade setuptools pip

echo -e ">> [`date`] installing requirements"
pip install -r requirements.txt

echo -e ">> [`date`] installing companion components"
pip install -r companion.txt

echo -e ">> [`date`] installing mpi4py"
pip install mpi4py

echo -e ">> [`date`] installing pycbc"
pip install .
