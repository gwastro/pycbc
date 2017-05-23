#!/bin/bash

set -ev

LOCAL=${PWD}

# working dir for downloaded dependencies
SRC=${HOME}/src
mkdir -p ${SRC}

# install dir for dependencies
INST=${HOME}/inst

export LD_LIBRARY_PATH=${INST}/lib:${INST}/lib64
export PKG_CONFIG_PATH=${INST}/lib/pkgconfig
export PATH=/usr/lib/ccache:${PATH}:${INST}/bin

# Update setuptools
pip install --upgrade pip setuptools

# Install needed version of numpy
pip install 'numpy==1.9.3' --upgrade 

if [ -f ${INST}/dep_install.done ]
then
    echo "Found cache of installed dependencies, using it"
else

    # install the version of swig that for some reason we have to use

    cd ${SRC}
    wget -q http://download.sourceforge.net/project/swig/swig/swig-2.0.11/swig-2.0.11.tar.gz
    tar -xzf swig-2.0.11.tar.gz
    cd swig-2.0.11
    ./configure -q --prefix=${INST}
    make -j
    make install

    # Install metaio

    cd ${SRC}
    wget -q https://www.lsc-group.phys.uwm.edu/daswg/download/software/source/metaio-8.2.tar.gz
    tar -xzf metaio-8.2.tar.gz
    cd metaio-8.2
    CPPFLAGS=-std=gnu99 ./configure -q --prefix=${INST}
    make -j
    make install

    # install framel

    cd ${SRC}
    wget -q http://lappweb.in2p3.fr/virgo/FrameL/v8r26.tar.gz
    tar -xzf v8r26.tar.gz
    cd v8r26
    autoreconf
    ./configure -q --prefix=${INST}
    make -j
    make install

    # Install lalsuite

    cd ${SRC}
    git clone -q https://github.com/lscsoft/lalsuite.git
    cd lalsuite
    # This sets the test release to https://versions.ligo.org/cgit/lalsuite/commit/?id=a2a5a476d33f169b8749e2840c306a48df63c936
    git checkout a2a5a476d33f169b8749e2840c306a48df63c936

    ./00boot
    ./configure -q --prefix=${INST} --enable-swig-python \
        --disable-lalstochastic --disable-lalinference --disable-laldetchar \
        --disable-lalxml --disable-lalburst --disable-lalapps
    make -j
    make install

    # run lalsimulation tests
    cd lalsimulation
    make check

    touch ${INST}/dep_install.done

    cd ${LOCAL}
fi

source ${INST}/etc/lal-user-env.sh

# Scipy would be required to build scikit-learn but that does not work with travis currently
#pip install 'scipy==0.16.0' --upgrade

# Install Pegasus
pip install http://download.pegasus.isi.edu/pegasus/4.7.2/pegasus-python-source-4.7.2.tar.gz

# Needed by mock 
pip install 'setuptools==18.2' --upgrade

python setup.py install
