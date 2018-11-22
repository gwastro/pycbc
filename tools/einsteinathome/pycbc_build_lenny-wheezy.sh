#!/bin/bash

# script to build PyCBC suite on Debian 7.0 Wheezy
# and probably any Linux system where Python 2.7 can be installed from packages
# With additions to work on 5.0 (Lenny) as well, where Python 2.7 needs to be compiled

# packages to install on the system:
# screen emacs git-core ant gcc g++ gfortran automake autoconf make libtool pkg-config
# libfftw3-dev libgsl0-dev libpcre3-dev libfreetype6-dev libjpeg-dev libpng-dev libhdf5-dev libmysqlclient-dev libpq-dev liblapack-dev
# python python-pip python-virtualenv python-dev
# libhdf5-dev is not available on Debian Lenny

set -e

if test "v`cat /etc/debian_version 2>/dev/null`" = "v4.0"; then
  echo -e "\\n\\n>> [`date`] Using Debian 5.0 (Lenny) settings"
  lenny=true
else
  lenny=false
fi

echo "[`date`] Start $0"

PYCBC="$PWD/pycbc"
SOURCE="$PYCBC/source"
PYTHON_PREFIX="$PYCBC"
ENVIRONMENT="$PYCBC/environment"
PREFIX="$ENVIRONMENT"
PATH="$PYTHON_PREFIX/bin:$PREFIX/bin:$PATH"
export FC=gfortran
libgfortran="`$FC -print-file-name=libgfortran.so|sed 's%/[^/]*$%%'`"
export LD_LIBRARY_PATH="$PYTHON_PREFIX/lib:$PREFIX/lib:$libgfortran:$LD_LIBRARY_PATH"
export PKG_CONFIG_PATH="$PYTHON_PREFIX/lib/pkgconfig:$PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH"
export LIBS="$LIBS -lgfortran"

#export CPPFLAGS="-I$PREFIX/include"
#export LDFLAGS="-L$PREFIX/lib -L$libgfortran"
#export LDFLAGS="-L$libgfortran $LDFLAGS"
# -static-libgfortran

# don't check certificates, these are usually outdated on older systems
export GIT_SSL_NO_VERIFY=true
wget_opts="-c --passive-ftp --no-check-certificate"
pip_install="install --trusted-host pypi.python.org"

# try to find lalsuite root
LALSUITE="`echo '$PWD/$0'|sed 's%.*//%/%;s%/[^/]*$%%'`/../../../.."
if ls "$LALSUITE/lal/00boot" >/dev/null 2>&1 ; then
  echo found lalsuite at $LALSUITE
else
  LALSUITE=""
fi

mkdir -p "$SOURCE"
cd "$SOURCE"

# git://versions.ligo.org/lalsuite.git
# test ".$LALSUITE" = "." || 
if test -d lalsuite/.git; then
  :
else
  git clone git://gitmaster.atlas.aei.uni-hannover.de/einsteinathome/lalsuite.git
  cd lalsuite
  git checkout -b eah_cbc origin/eah_cbc
  cd ..
fi

if $lenny; then

  # build a local python-2.7 and set up a virtal environment
  # build hdfs (not available), gsl (too old) and zlib (missing zlib.pc) while at it

  # PYTHON
  v=2.7.10
  p=Python-$v
  echo "[`date`] building $p"
  test -r $p.tgz || wget $wget_opts http://www.python.org/ftp/python/$v/$p.tgz
  rm -rf $p
  tar -xzf $p.tgz
  cd $p
  ./configure --enable-shared --prefix="$PYTHON_PREFIX"
  make
  make install
  cd ..
  python -m ensurepip
  echo "[`date`] pip install --upgrade pip"
  pip install --upgrade pip
  echo "[`date`] pip install virtualenv"
  pip $pip_install virtualenv
  unset PYTHONPATH
  rm -rf "$ENVIRONMENT"
  virtualenv "$ENVIRONMENT"
  source "$ENVIRONMENT/bin/activate"
  # workaround to make the virtualenv accept .pth files
  export PYTHONPATH="$PREFIX/lib/python2.7/site-packages:$PYTHONPATH"

  # GSL
  p=gsl-1.16
  echo "[`date`] building $p"
  test -r $p.tar.gz || wget $wget_opts ftp://ftp.fu-berlin.de/unix/gnu/gsl/$p.tar.gz
  rm -rf $p
  tar -xzf $p.tar.gz
  cd $p
  ./configure --enable-shared --enable-static --prefix="$PREFIX"
  make
  make install
  cd ..

  # FFTW
  p=fftw-3.3.3
  echo "[`date`] building $p"
  test -r $p.tar.gz || wget $wget_opts http://www.aei.mpg.de/~bema/$p.tar.gz
  rm -rf $p
  tar -xzf $p.tar.gz
  cd $p
  ./configure --enable-shared --enable-static --prefix="$PREFIX"
  make
  make install
  ./configure --enable-shared --enable-static --prefix="$PREFIX" --enable-float
  make
  make install
  cd ..

  # ZLIB
  p=zlib-1.2.8
  echo "[`date`] building $p"
  test -r $p.tar.gz || wget $wget_opts http://www.aei.mpg.de/~bema/$p.tar.gz
  rm -rf $p
  tar -xzf $p.tar.gz
  cd $p
  ./configure --prefix=$PREFIX
  make
  make install
  echo 'prefix=/Users/bema/EinsteinAtHome/EinsteinAtHome/install
exec_prefix=${prefix}
libdir=${exec_prefix}/lib
includedir=${prefix}/include
Name: ZLIB
Description: zlib Compression Library
Version: 1.2.3
Libs: -L${libdir} -lz
Cflags: -I${includedir}' |
  sed "s%^prefix=.*%prefix=$PREFIX%;s/^Version: .*/Version: $p/;s/^Version: zlib-/Version: /" > "$PREFIX/lib/pkgconfig/zlib.pc"
  cd ..

  # HDF5
  p=hdf5-1.8.12
  echo "[`date`] building $p"
  test -r $p.tar.gz || wget $wget_opts https://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.12/src/$p.tar.gz
  rm -rf $p
  tar -xzf $p.tar.gz
  cd $p
  ./configure --enable-shared --enable-static --prefix="$PREFIX"
  make
  make install
  echo 'prefix=
exec_prefix=${prefix}
includedir=${prefix}/include
libdir=${exec_prefix}/lib
Name: hdf5
Description: HDF5
Version: 1.8.12
Requires.private: zlib
Cflags: -I${includedir}
Libs: -L${libdir} -lhdf5' |
  sed "s%^prefix=.*%prefix=$PREFIX%" > "$PREFIX/lib/pkgconfig/hdf5.pc"
  cd ..

else # if $lenny

  # set up virtual environment
  unset PYTHONPATH
  rm -rf "$ENVIRONMENT"
  virtualenv "$ENVIRONMENT"
  source "$ENVIRONMENT/bin/activate"
  # workaround to make the virtualenv accept .pth files
  export PYTHONPATH="$PREFIX/lib/python2.7/site-packages:$PYTHONPATH"

  # upgrade pip in virtual environment, don't use $pip_install options here
  echo "[`date`] pip install --upgrade pip"
  pip install --upgrade pip

fi # if $lenny

echo "[`date`] pip install --upgrade distribute"
pip $pip_install --upgrade distribute

# SCIPY
echo "[`date`] pip install numpy==1.9.3"
pip $pip_install numpy==1.9.3
echo "[`date`] pip install scipy==0.16.0"
pip $pip_install scipy==0.16.0

# LIBFRAME
p=libframe-8.21
echo "[`date`] building $p"
test -r $p.tar.gz || wget $wget_opts http://lappweb.in2p3.fr/virgo/FrameL/$p.tar.gz
rm -rf $p
tar -xzf $p.tar.gz
cd $p
./configure --enable-shared --enable-static --prefix="$PREFIX"
make
make install
sed "s%^prefix=.*%prefix=$PREFIX%" src/libframe.pc > $PREFIX/lib/pkgconfig/libframe.pc
cd ..

# METAIO
p=metaio-8.3.0
echo "[`date`] building $p"
test -r $p.tar.gz || wget $wget_opts https://www.lsc-group.phys.uwm.edu/daswg/download/software/source/$p.tar.gz
rm -rf $p
tar -xzf $p.tar.gz
cd $p
./configure --enable-shared --enable-static --prefix="$PREFIX"
make
make install
cd ..

# SWIG
p=swig-3.0.7
echo "[`date`] building $p"
test -r $p.tar.gz || wget $wget_opts http://atlas3.atlas.aei.uni-hannover.de/~bema/tarballs/$p.tar.gz
rm -rf $p
tar -xzf $p.tar.gz
cd $p
./configure --prefix=$PREFIX --without-tcl --with-python --without-python3 --without-perl5 --without-octave --without-scilab --without-java --without-javascript --without-gcj --without-android --without-guile --without-mzscheme --without-ruby --without-php --without-ocaml --without-pike --without-chicken --without-csharp --without-lua --without-allegrocl --without-clisp --without-r --without-go --without-d 
make
make install
cd ..

# LALSUITE
echo "[`date`] building lalsuite"
cd lalsuite
echo "[`date`] git HEAD: `git log -1 --pretty=oneline --abbrev-commit`"
./00boot
cd ..
rm -rf lalsuite-build
mkdir lalsuite-build
cd lalsuite-build
../lalsuite/configure --disable-gcc-flags --enable-shared --enable-static --enable-swig-python --prefix="$PREFIX" --disable-lalxml --disable-lalpulsar --disable-laldetchar --disable-lalstochastic --disable-lalinference
make
make install
for i in $PREFIX/etc/*-user-env.sh; do
  source "$i"
done
echo "[`date`] building PyLAL"
cd ../lalsuite/pylal
python setup.py install --prefix="$PREFIX"
echo "[`date`] building GLUE"
cd ../glue
python setup.py install --prefix="$PREFIX"
cd ../..
test -r "$PREFIX/etc/pylal-user-env.sh" && source "$PREFIX/etc/pylal-user-env.sh"
test -r "$PREFIX/etc/glue-user-env.sh" && source "$PREFIX/etc/glue-user-env.sh"

# MPLD
p=mpld3-0.3git
echo "[`date`] installing $p"
pip $pip_install "https://github.com/ligo-cbc/mpld3/tarball/master#egg=$p"

# PEGASUS
v=4.5.2
p=pegasus-source-$v
echo "[`date`] installing $p"
test -r $p.tar.gz || wget $wget_opts http://download.pegasus.isi.edu/pegasus/$v/$p.tar.gz
tar xzf $p.tar.gz
cd $p
ant dist-python-source 
pushd lib/pegasus/python/
echo "[`date`] installing dependencies for $p"
pip $pip_install -r pegasus_wms.egg-info/requires.txt
echo "[`date`] building $p"
python setup.py install --prefix="$PREFIX"
popd
cd ..

# PYCBC
# pip install git+https://github.com/ligo-cbc/pycbc@v1.1.0#egg=pycbc --process-dependency-links
echo "[`date`] pip install pycbc==1.2.0"
pip $pip_install pycbc==1.2.0

test -r "$PREFIX/etc/pycbc-user-env.sh" && source "$PREFIX/etc/pycbc-user-env.sh"

# PYINSTALLER
echo "[`date`] pip install pyinstaller"
pip $pip_install pyinstaller

# env|tee pycbc.env

# TEST
cd ..
echo "[`date`] testing"
cd $PREFIX
./bin/pycbc_inspiral --help

# BUNDLE DIR
echo "[`date`] running pyinstaller"
rm -rf dist/pycbc_inspiral_d dist/pycbc_inspiral
pyinstaller --hidden-import=scipy.linalg.cython_blas --hidden-import=scipy.linalg.cython_lapack --hidden-import=pkg_resources --onedir ./bin/pycbc_inspiral
cd dist/pycbc_inspiral

# TEST BUNDLE
echo "[`date`] testing"
./pycbc_inspiral --help
cd ..
mv pycbc_inspiral pycbc_inspiral_d
cd ..

# BUNDLE FILE
echo "[`date`] running pyinstaller"
pyinstaller --hidden-import=scipy.linalg.cython_blas --hidden-import=scipy.linalg.cython_lapack --hidden-import=pkg_resources --onefile ./bin/pycbc_inspiral

# TEST BUNDLE
echo "[`date`] testing"

# leave virtual environment
deactivate

./dist/pycbc_inspiral --help

echo "[`date`] Success $0"
