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
# FIXME: For now dqsegdb needs a glue, which it does not declare. Please remove
#        when no longer needed.
pip install pycbc-glue

# install the packges needed to build the documentation
pip install "Sphinx>=1.4.2"
pip install numpydoc
pip install sphinx-rtd-theme
pip install git+https://github.com/ligo-cbc/sphinxcontrib-programoutput.git#egg=sphinxcontrib-programoutput

# get library needed to build documentation
wget_opts="-c --passive-ftp --no-check-certificate --tries=5 --timeout=30"
url="https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-software/download/b3680bfb627a7350f29d31c7d91c4e09ff8a9fdc/x86_64/composer_xe_2015.0.090"
p="libmkl_rt.so"
pushd ${BUILD}/pycbc-sources
test -r $p || wget $wget_opts ${url}/${p}
cp -v $p $PREFIX/lib/$p
chmod +x $PREFIX/lib/$p
popd

# re-install pycbc
python setup.py install
