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

touch ${INST}/dep_install.done

if [ -f ${INST}/dep_install.done ]
then
    echo "Found cache of installed dependencies, using it"
else
    # Install needed version of numpy
    pip install 'numpy==1.9.3' --upgrade 

    pip install http://download.pegasus.isi.edu/pegasus/4.5.2/pegasus-python-source-4.5.2.tar.gz
fi

/bin/true

#source ${INST}/etc/lal-user-env.sh

#python setup.py install
