#!/bin/bash
set -e

LOCAL=$PWD
mkdir -p $LOCAL/src

sudo apt-get install \
libfftw3-dev \
python-decorator \
python-jinja2 \
liblapack-dev \
python-scipy \
python-numpy \
gfortran \
libgsl0-dev

pip install decorator --upgrade
pip install matplotlib --upgrade

# install the version of swig that for some reason we have to use
wget http://downloads.sourceforge.net/project/swig/swig/swig-2.0.11/swig-2.0.11.tar.gz
tar -xvf swig-2.0.11.tar.gz
cd swig-2.0.11; ./configure --prefix=$LOCAL; make -j; make install; cd ..

export LD_LIBRARY_PATH=$LOCAL/lib:$LOCAL/lib64
export PKG_CONFIG_PATH=$LOCAL/lib/pkgconfig
export PATH=$PATH:$LOCAL/bin


# Standard python dependencies
pip install pillow 
pip install -e 'git+http://github.com/jakevdp/mpld3.git#egg=mpld3-0.3'
SWIG_FEATURES="-cpperraswarn -includeall -I/usr/include/openssl" pip install M2Crypto

# Install metaio
wget https://www.lsc-group.phys.uwm.edu/daswg/download/software/source/metaio-8.2.tar.gz
tar -zxvf metaio-8.2.tar.gz
cd metaio-8.2; CPPFLAGS=-std=gnu99 ./configure --prefix=$LOCAL; make -j; make install; cd ..

# install framel
wget http://lappweb.in2p3.fr/virgo/FrameL/v8r26.tar.gz
tar -xvf v8r26.tar.gz 
cd v8r26; autoreconf; ./configure --prefix=$LOCAL;make -j; make install; cd ..


# Install lalsuite itself
cd $LOCAL/src/
git clone https://github.com/ahnitz/lalsuite.git
cd lalsuite
./00boot; ./configure --prefix=$LOCAL --enable-swig-python; make -j; make install
source $LOCAL/etc/lal-user-env.sh
cd pylal; python setup.py install --prefix=$LOCAL; cd ..
cd glue; python setup.py install --prefix=$LOCAL; cd ..
cd ..

echo source $LOCAL/etc/glue-user-env.sh >> $TRAVIS_BUILD_DIR/source
echo source $LOCAL/etc/lal-user-env.sh >> $TRAVIS_BUILD_DIR/source
echo $PWD
chmod 755 $TRAVIS_BUILD_DIR/source
