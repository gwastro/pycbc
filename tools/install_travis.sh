#!/bin/bash

set -ev

# determine the pycbc git branch and origin
git branch -vvv
if test x$TRAVIS_PULL_REQUEST = "xfalse" ; then
    PYCBC_CODE="--pycbc-commit=${TRAVIS_COMMIT}"
else
    PYCBC_CODE="--pycbc-fetch-ref=refs/pull/${TRAVIS_PULL_REQUEST}/merge"
fi

# set the lalsuite checkout to use
LALSUITE_CODE="--lalsuite-commit=a2a5a476d33f169b8749e2840c306a48df63c936"
# LALSUITE_CODE="--lalsuite-commit=master" --clean-lalsuite

# store the travis test directory
LOCAL=${PWD}

# create working dir for build script
BUILD=${HOME}/build
mkdir -p ${BUILD}
export PYTHONUSERBASE=${BUILD}/.local
export XDG_CACHE_HOME=${BUILD}/.cache

# run the einstein at home build and test script
pushd ${BUILD}
${LOCAL}/tools/einsteinathome/pycbc_build_eah.sh ${LALSUITE_CODE} ${PYCBC_CODE} --clean-pycbc --silent-build
popd

# setup the pycbc environment to install the additional software needed
. tools/travis_env.sh

# update setuptools
pip install --upgrade pip setuptools

# needed by mock 
pip install 'setuptools==18.2' --upgrade

# install pegasus
pip install http://download.pegasus.isi.edu/pegasus/4.7.4/pegasus-python-source-4.7.4.tar.gz

# install M2Crypto
SWIG_FEATURES="-cpperraswarn -includeall -I/usr/include/openssl" pip install M2Crypto

# install the segment database tools
pip install git+https://github.com/ligovirgo/dqsegdb@clean_pip_install_1_4_1#egg=dqsegdb

# re-install pycbc
python setup.py install
