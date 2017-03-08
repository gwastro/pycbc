#!/bin/bash

set -ev

LOCAL=${PWD}

# working dir for build script
BUILD=${HOME}/build
mkdir -p ${BUILD}

pushd ${BUILD}
${LOCAL}/tools/einsteinathome/pycbc_build_eah.sh --lalsuite-commit=a2a5a476d33f169b8749e2840c306a48df63c936 --clean-pycbc
popd

# Update setuptools
#pip install --upgrade pip setuptools

# Needed by mock 
#pip install 'setuptools==18.2' --upgrade

#python setup.py install
