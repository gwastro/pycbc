#!/bin/bash

set -ev

LOCAL=${PWD}

# working dir for downloaded dependencies
SRC=${HOME}/src
mkdir -p ${SRC}

# install dir for dependencies
# FIXME try caching this!
INST=${HOME}/inst

export LD_LIBRARY_PATH=${INST}/lib:${INST}/lib64
export PKG_CONFIG_PATH=${INST}/lib/pkgconfig
export PATH=/usr/lib/ccache:${PATH}:${INST}/bin

# install the version of swig that for some reason we have to use

cd ${SRC}
# FIXME SF mirror hardcoded to heanet to work around occasional failures
wget -q http://heanet.dl.sourceforge.net/project/swig/swig/swig-2.0.11/swig-2.0.11.tar.gz
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
git clone -q https://github.com/ahnitz/lalsuite.git
cd lalsuite
./00boot
./configure -q --prefix=${INST} --enable-swig-python \
--disable-lalstochastic --disable-lalinference --disable-laldetchar \
--disable-lalxml --disable-lalburst --disable-lalapps
make -j
make install

cd ${LOCAL}

source ${INST}/etc/lal-user-env.sh

python setup.py install
