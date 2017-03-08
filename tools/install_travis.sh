#!/bin/bash

set -ev

# determine the git branch and origin
GIT_BRANCH=`git name-rev --name-only HEAD`
echo "Testing branch ${GIT_BRANCH}"
GIT_ORIGIN=`git config --get remote.origin.url`
echo "Testing from ${GIT_ORIGIN}"

git remote -vvv
git branch -vvv

# store the travis test directory
LOCAL=${PWD}

# create working dir for build script
BUILD=${HOME}/build
mkdir -p ${BUILD}

# run the einstein at home build and test script
pushd ${BUILD}
${LOCAL}/tools/einsteinathome/pycbc_build_eah.sh --lalsuite-commit=a2a5a476d33f169b8749e2840c306a48df63c936 --clean-pycbc
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
