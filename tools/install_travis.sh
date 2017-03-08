#!/bin/bash

set -ev

# determine the git branch and origin
git branch -vvv
PYCBC_MERGE_COMMIT=`git rev-parse HEAD`

# store the travis test directory
LOCAL=${PWD}

# create working dir for build script
BUILD=${HOME}/build
mkdir -p ${BUILD}

# run the einstein at home build and test script
pushd ${BUILD}
${LOCAL}/tools/einsteinathome/pycbc_build_eah.sh --lalsuite-commit=a2a5a476d33f169b8749e2840c306a48df63c936 --clean-pycbc --pycbc-commit=${PYCBC_MERGE_COMMIT}
popd

# setup the pycbc environment to run the additional travis tests
source ${BUILD}/pycbc-build/environment/etc/lalsuite-user-env.sh
source ${BUILD}/pycbc-build/environment/bin/activate
export PYTHONUSERBASE=${BUILD}/.local
export XDG_CACHE_HOME=${BUILD}/.cache

# update setuptools
pip install --upgrade pip setuptools

# needed by mock 
pip install 'setuptools==18.2' --upgrade

# install pegasus
pip install http://download.pegasus.isi.edu/pegasus/4.7.2/pegasus-python-source-4.7.4.tar.gz

# re-install pycbc
python setup.py install
