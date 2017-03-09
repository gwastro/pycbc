#!/bin/bash

set -ev

# determine the git branch and origin
git branch -vvv
PYCBC_MERGE_REF=`cut -f3 .git/FETCH_HEAD | cut -d " " -f1 | tr -d "'"`

env | grep TRAVIS
exit 1

# store the travis test directory
LOCAL=${PWD}

# create working dir for build script
BUILD=${HOME}/build
mkdir -p ${BUILD}
export PYTHONUSERBASE=${BUILD}/.local
export XDG_CACHE_HOME=${BUILD}/.cache

# run the einstein at home build and test script
pushd ${BUILD}
${LOCAL}/tools/einsteinathome/pycbc_build_eah.sh --lalsuite-commit=a2a5a476d33f169b8749e2840c306a48df63c936 --pycbc-fetch-ref=${PYCBC_MERGE_REF} --clean-pycbc --silent-build
popd

# setup the pycbc environment to run the additional travis tests
BUILDDIRNAME="pycbc-build"
PYCBC="$BUILD/$BUILDDIRNAME"
PYTHON_PREFIX="$PYCBC"
ENVIRONMENT="$PYCBC/environment"
PREFIX="$ENVIRONMENT"
PATH="$PREFIX/bin:$PYTHON_PREFIX/bin:$PATH"
export LD_LIBRARY_PATH="$PREFIX/lib:$PREFIX/bin:$PYTHON_PREFIX/lib:/usr/local/lib:$LD_LIBRARY_PATH"
export PKG_CONFIG_PATH="$PREFIX/lib/pkgconfig:$PYTHON_PREFIX/lib/pkgconfig:/usr/local/lib/pkgconfig:$PKG_CONFIG_PATH"
source ${BUILD}/pycbc-build/environment/etc/lalsuite-user-env.sh
source ${BUILD}/pycbc-build/environment/bin/activate

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
