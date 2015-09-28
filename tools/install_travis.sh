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
    git checkout lalsuite_o1_branch

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

# Install needed version of numpy
pip install 'numpy==1.9.3' --upgrade 

# Install Pegasus
pip install http://download.pegasus.isi.edu/pegasus/4.5.2/pegasus-python-source-4.5.2.tar.gz


python setup.py install
