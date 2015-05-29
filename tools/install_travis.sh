#!/bin/bash
set -e

LOCAL=$PWD

mkdir -p $LOCAL/src

pip install distribute --upgrade

pipi() {
    pip install --ignore-installed $@
}

pipn() {
    pip install --upgrade --no-deps $@
    pip install $@
}

# Standard python dependencies
pipi 'Mako>=1.0.1' 
pipi pillow 
pipi 'decorator>=3.4.2' 
pipi 'numpy>=1.6.4'
pipi jinja2
pipi 'scipy>=0.13.0'
pipi 'matplotlib>=1.3.0'
pipi -e 'git+http://github.com/jakevdp/mpld3.git#egg=mpld3-0.3'
pipi 'argparse>=1.3.0'
SWIG_FEATURES="-cpperraswarn -includeall -I/usr/include/openssl" pip install M2Crypto

# Install lalsuite itself
cd $LOCAL/src/
wget https://github.com/ahnitz/lalsuite/tarball/master
tar -xvf *lalsuite*.tar.gz
cd lalsuite
./00boot
./configure --prefix=$LOCAL --enable-swig-python
make -j install
source $LOCAL/etc/lal-user-env.sh

# LAL python dependencies
pipn git+https://github.com/ahnitz/glue.git#egg=glue
pipn git+https://github.com/ahnitz/pylal.git#egg=pylal
cd $LOCAL

echo source $LOCAL/etc/glue-user-env.sh >> source
echo source $LOCAL/etc/lal-user-env.sh >> source
echo source $LOCAL/etc/pycbc-user-env.sh >> source
chmod 755 source
