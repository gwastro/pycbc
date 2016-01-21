#!/bin/bash

# script to build PyCBC suite on Debian Etch and Windows Cygwin

# FIXME/todo:
# - 3 [main] python 2664 child_info_fork::abort: couldn't allocate memory 0x1F3D90000(0) for '_rank.dll' alignment, Win32 error 87
#   - _rank.dll is a scipy library
#   - Cygwin's "rebaseall" doesn't help as suggested, because _rank.dll is built afterwards
#   - Cygwin "rebase" may help http://pipeline.lbl.gov/code/3rd_party/licenses.win/Cygwin/rebase-3.0.README
# - pycbc-glue build: ".tp_base = NULL"

set -e

echo -e "\\n\\n>> [`date`] Start $0"

if [ ".$DEBUG" = "." ]; then
  cleanup=true
else
  cleanup=false
fi

if test "v`cat /etc/debian_version 2>/dev/null`" = "v4.0"; then
  echo -e "\\n\\n>> [`date`] Using Debian 4.0 (etch) settings"
  shared="--enable-shared"
  build_dlls=false
  build_ssl=true
  build_pyssl=true
  build_python=false
  build_lapack=true
  compile_numpy=false
  compile_scipy=false
  build_hdf5=true
  build_freetype=true
  build_libpq=false
  build_gsl=true
  build_swig=true
  fake_psycopg26=true
  compile_pycbc_glue=false
  build_pegasus_source=true
  pip_install_pycbc=true
  build_pyinstaller=false
  build_pyinstaller_git=false
  build_preinst_before_lalsuite=false
  add_python_runtime_options=false
else
  echo -e "\\n\\n>> [`date`] Using Cygwin settings"
  shared="--enable-shared"
  build_dlls=true
  build_ssl=false
  build_pyssl=true
  build_python=false
  build_lapack=true
  compile_numpy=true
  compile_scipy=true
  build_hdf5=false
  build_freetype=false
  build_libpq=false
  build_gsl=false
  build_swig=false
  compile_pycbc_glue=true
  fake_psycopg26=true
  build_pegasus_source=false
  pip_install_pycbc=true
  build_pyinstaller=true
  build_pyinstaller_git=false
  build_preinst_before_lalsuite=true
  add_python_runtime_options=false
fi

PYCBC="$PWD/pycbc"
SOURCE="$PYCBC/source"
PYTHON_PREFIX="$PYCBC"
ENVIRONMENT="$PYCBC/environment"
PREFIX="$ENVIRONMENT"
PATH="$PREFIX/bin:$PYTHON_PREFIX/bin:$PATH"
export FC=gfortran
libgfortran="`$FC -print-file-name=libgfortran.so|sed 's%/[^/]*$%%'`"
export LD_LIBRARY_PATH="$PREFIX/lib:$PREFIX/bin:$PYTHON_PREFIX/lib:$libgfortran:/usr/local/lib:$LD_LIBRARY_PATH"
export PKG_CONFIG_PATH="$PREFIX/lib/pkgconfig:$PYTHON_PREFIX/lib/pkgconfig:/usr/local/lib/pkgconfig:$PKG_CONFIG_PATH"
export LIBS="$LIBS -lgfortran"

#export CPPFLAGS="-I$PREFIX/include"
#export LDFLAGS="-L$PREFIX/lib -L$libgfortran"
#export LDFLAGS="-L$libgfortran $LDFLAGS"
# -static-libgfortran

export GIT_SSL_NO_VERIFY=true
wget_opts="-c --passive-ftp --no-check-certificate"
pip_install="install --trusted-host pypi.python.org --trusted-host github.com"

# try to find lalsuite root
LALSUITE="`echo '$PWD/$0'|sed 's%.*//%/%;s%/[^/]*$%%'`/../../../.."
if ls "$LALSUITE/lal/00boot" >/dev/null 2>&1 ; then
  echo found lalsuite at $LALSUITE
else
  LALSUITE=""
fi

echo "export PATH='$PATH'"
echo "export LD_LIBRARY_PATH='$LD_LIBRARY_PATH'"
echo "export PKG_CONFIG_PATH='$PKG_CONFIG_PATH'"

if test -r pycbc-preinst.tgz -o -r pycbc-preinst-lalsuite.tgz; then

  rm -rf pycbc
  if test -r pycbc-preinst-lalsuite.tgz; then
    echo -e "\\n\\n>> [`date`] using pycbc-preinst-lalsuite.tgz"
    tar -xzf pycbc-preinst-lalsuite.tgz
  else
    echo -e "\\n\\n>> [`date`] using pycbc-preinst.tgz"
    tar -xzf pycbc-preinst.tgz
  fi
  # set up virtual environment
  unset PYTHONPATH
  source "$ENVIRONMENT/bin/activate"
  # workaround to make the virtualenv accept .pth files
  export PYTHONPATH="$PREFIX/lib/python2.7/site-packages:$PYTHONPATH"
  cd "$SOURCE"
  if ls $PREFIX/etc/*-user-env.sh >/dev/null 2>&1; then
    for i in $PREFIX/etc/*-user-env.sh; do
      source "$i"
    done
  fi

else # if pycbc-preinst.tgz

  mkdir -p "$SOURCE"
  cd "$SOURCE"

  if test -d lalsuite/.git; then
    :
  else
    git clone git://gitmaster.atlas.aei.uni-hannover.de/einsteinathome/lalsuite.git
    cd lalsuite
    git checkout -b eah_cbc origin/eah_cbc
    cd ..
  fi

  # OpenSSL
  if $build_ssl; then
    # p=openssl-1.0.2e # compile error on pyOpenSSL 0.13:
    # pycbc/include/openssl/x509.h:751: note: previous declaration X509_REVOKED_ was here
    p=openssl-1.0.1p
    echo -e "\\n\\n>> [`date`] building $p"
    test -r $p.tar.gz || wget $wget_opts http://www.openssl.org/source/$p.tar.gz
    rm -rf $p
    tar -xzvf $p.tar.gz  &&
    cd $p &&
    ./config shared -fPIC "--prefix=$PYTHON_PREFIX" #  no-shared no-sse2 CFLAGS=-fPIC
    make
    make install
    cd ..
    $cleanup && rm -rf $p
  fi

  # PYTHON
  if $build_python; then
    v=2.7.10
    p=Python-$v
    echo -e "\\n\\n>> [`date`] building $p"
    test -r $p.tgz || wget $wget_opts http://www.python.org/ftp/python/$v/$p.tgz
    rm -rf $p
    tar -xzf $p.tgz
    cd $p
    ./configure $shared --prefix="$PYTHON_PREFIX"
    make
    make install
    cd ..
    $cleanup && rm -rf $p
    python -m ensurepip
    echo -e "\\n\\n>> [`date`] pip install --upgrade pip"
    pip install --upgrade pip
    echo -e "\\n\\n>> [`date`] pip install virtualenv"
    pip $pip_install virtualenv
  fi

  # set up virtual environment
  unset PYTHONPATH
  rm -rf "$ENVIRONMENT"
  virtualenv "$ENVIRONMENT"
  source "$ENVIRONMENT/bin/activate"
  # workaround to make the virtualenv accept .pth files
  export PYTHONPATH="$PREFIX/lib/python2.7/site-packages:$PYTHONPATH"

  # pyOpenSSL-0.13
  if $build_pyssl; then
    p=pyOpenSSL-0.13
    echo -e "\\n\\n>> [`date`] building $p"
    test -r $p.tar.gz || wget $wget_opts "http://pypi.python.org/packages/source/p/pyOpenSSL/$p.tar.gz"
    rm -rf $p
    tar -xzf $p.tar.gz
    cd $p
    sed -i~ 's/X509_REVOKED_dup/X509_REVOKED_dup_static/' OpenSSL/crypto/crl.c
    python setup.py build_ext "-I$PYTHON_PREFIX/include" "-L$PYTHON_PREFIX/lib"
    python setup.py build
    python setup.py install --prefix="$PREFIX"
    cd ..
    $cleanup && rm -rf $p
  else
    echo -e "\\n\\n>> [`date`] pip install pyOpenSSL==0.13"
    pip $pip_install pyOpenSSL==0.13
  fi

  # http://sourceforge.net/projects/math-atlas/files/Stable/3.10.2/atlas3.10.2.tar.bz2

  # LAPACK & BLAS
  if $build_lapack; then
    p=lapack-3.6.0
    echo -e "\\n\\n>> [`date`] building $p"
    test -r $p.tgz || wget $wget_opts http://www.netlib.org/lapack/$p.tgz
    rm -rf $p
    tar -xzf $p.tgz
    cd $p
    # configure: compile with -fPIC, remove -frecoursive, build deprecated functions
    sed "s/ *-frecursive//;s/gfortran/$FC -fPIC/;s/^#MAKEDEPRECATED.*/BUILD_DEPRECATED = Yes/" make.inc.example > make.inc
    make lapack_install lib blaslib
    mkdir -p "$PREFIX/lib"
    cp lib*.a "$PREFIX/lib"
    cp librefblas.a "$PREFIX/lib/libblas.a"
    cd ..
    $cleanup && rm -rf $p
  fi
  
  # NUMPY
  if $compile_numpy; then
    p=numpy-1.9.3
    echo -e "\\n\\n>> [`date`] building $p"
    test -r $p.tar.gz || wget $wget_opts https://pypi.python.org/packages/source/n/numpy/$p.tar.gz
    rm -rf $p
    tar -xzf $p.tar.gz 
    cd $p
    python setup.py build --fcompiler=$FC
    python setup.py install --prefix=$PREFIX
    cd ..
    $cleanup && rm -rf $p
  else
    echo -e "\\n\\n>> [`date`] pip install numpy==1.9.3"
    pip $pip_install numpy==1.9.3
  fi

  echo -e "\\n\\n>> [`date`] pip install nose"
  pip $pip_install nose
  echo -e "\\n\\n>> [`date`] pip install Cython==0.23.2"
  pip $pip_install Cython==0.23.2

  # SCIPY
  if $compile_scipy; then
    p=scipy-0.16.0
    echo -e "\\n\\n>> [`date`] building $p"
    if test -d scipy/.git; then
      cd scipy
    else
      git clone https://github.com/scipy/scipy.git
      cd scipy
      git checkout v0.16.0
      git cherry-pick 832baa20f0b5d521bcdf4784dda13695b44bb89f
    fi
    python setup.py build --fcompiler=$FC
    python setup.py install --prefix=$PREFIX
    cd ..
  else
    echo -e "\\n\\n>> [`date`] pip install scipy==0.16.0"
    pip $pip_install scipy==0.16.0
  fi

  # this test will catch scipy build errors that would not emerge before running pycbc_inspiral
  echo -e "\\n\\n>> [`date`] Testing: python -c 'from scipy.io.wavfile import write as write_wav'"
  python -c 'from scipy.io.wavfile import write as write_wav'
# python -c 'import scipy; scipy.test(verbose=2);'

  # LIBPQ
  if $build_libpq; then
    v=9.2.0 # 8.4.22 # 9.3.10
    p=postgresql-$v
    echo -e "\\n\\n>> [`date`] building $p"
    test -r $p.tar.gz || wget $wget_opts "https://ftp.postgresql.org/pub/source/v$v/$p.tar.gz"
    rm -rf $p
    tar -xzf $p.tar.gz
    cd $p
    ./configure --without-readline --prefix="$PREFIX"
    cd src/interfaces/libpq
    make
    make install
    mkdir -p "$PREFIX/lib/pkgconfig"
    echo 'prefix=
exec_prefix=${prefix}
includedir=${prefix}/include
libdir=${exec_prefix}/lib
Name: libpq
Description: Postgres client lib
Version: 8.4.22
Cflags: -I${includedir}
Libs: -L${libdir} -lpq' |
  sed "s%^prefix=.*%prefix=$PREFIX%;s/^Version: .*/Version: $v/" > "$PREFIX/lib/pkgconfig/libpq.pc"
    cd ../../../..
    $cleanup && rm -rf $p
  fi

  # GSL
  if $build_gsl; then
    p=gsl-1.16
    echo -e "\\n\\n>> [`date`] building $p"
    test -r $p.tar.gz || wget $wget_opts ftp://ftp.fu-berlin.de/unix/gnu/gsl/$p.tar.gz
    rm -rf $p
    tar -xzf $p.tar.gz
    cd $p
    ./configure $shared --enable-static --prefix="$PREFIX"
    make
    make install
    cd ..
    $cleanup && rm -rf $p
  fi

  # FFTW
  p=fftw-3.3.3
  echo -e "\\n\\n>> [`date`] building $p"
  test -r $p.tar.gz || wget $wget_opts http://www.aei.mpg.de/~bema/$p.tar.gz
  rm -rf $p
  tar -xzf $p.tar.gz
  cd $p
  ./configure $shared --enable-static --prefix="$PREFIX"
  make
  make install
  ./configure $shared --enable-static --prefix="$PREFIX" --enable-float
  make
  make install
  cd ..
  $cleanup && rm -rf $p

  # ZLIB
  p=zlib-1.2.8
  echo -e "\\n\\n>> [`date`] building $p"
  test -r $p.tar.gz || wget $wget_opts http://www.aei.mpg.de/~bema/$p.tar.gz
  rm -rf $p
  tar -xzf $p.tar.gz
  cd $p
  ./configure --prefix=$PREFIX
  make
  make install
  mkdir -p "$PREFIX/lib/pkgconfig"
  echo 'prefix=
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
  $cleanup && rm -rf $p

  # HDF5
  if $build_hdf5; then
    p=hdf5-1.8.12
    echo -e "\\n\\n>> [`date`] building $p"
    test -r $p.tar.gz || wget $wget_opts https://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.12/src/$p.tar.gz
    rm -rf $p
    tar -xzf $p.tar.gz
    cd $p
    ./configure $shared --enable-static --prefix="$PREFIX"
    make
    make install
    mkdir -p "$PREFIX/lib/pkgconfig"
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
    $cleanup && rm -rf $p
  fi

  # FREETYPE
  if $build_freetype; then
    p=freetype-2.3.0
    echo -e "\\n\\n>> [`date`] building $p"
    test -r $p.tar.gz || wget $wget_opts http://download.savannah.gnu.org/releases/freetype/freetype-old/$p.tar.gz
    rm -rf $p
    tar -xzf $p.tar.gz
    cd $p
    ./configure $shared --enable-static --prefix="$PREFIX"
    make
    make install
    cd ..
    $cleanup && rm -rf $p
  fi

  echo -e "\\n\\n>> [`date`] pip install --upgrade distribute"
  pip $pip_install --upgrade distribute

  # LIBFRAME
  p=libframe-8.21
  echo -e "\\n\\n>> [`date`] building $p"
  test -r $p.tar.gz || wget $wget_opts http://lappweb.in2p3.fr/virgo/FrameL/$p.tar.gz
  rm -rf $p
  tar -xzf $p.tar.gz
  cd $p
  if $build_dlls; then
    for i in src/Makefile*; do
      echo 'libFrame_la_LDFLAGS += -no-undefined' >> $i
    done
  fi
  ./configure $shared --enable-static --prefix="$PREFIX"
  make
  make install
  mkdir -p "$PREFIX/lib/pkgconfig"
  sed "s%^prefix=.*%prefix=$PREFIX%" src/libframe.pc > $PREFIX/lib/pkgconfig/libframe.pc
  cd ..
  $cleanup && rm -rf $p

  # METAIO
  p=metaio-8.3.0
  echo -e "\\n\\n>> [`date`] building $p"
  test -r $p.tar.gz || wget $wget_opts https://www.lsc-group.phys.uwm.edu/daswg/download/software/source/$p.tar.gz
  rm -rf $p
  tar -xzf $p.tar.gz
  cd $p
  if $build_dlls; then
    for i in src/Makefile*; do
      echo 'libmetaio_la_LDFLAGS += -no-undefined' >> $i
    done
  fi
  ./configure $shared --enable-static --prefix="$PREFIX"
  make
  make install
  cd ..
  $cleanup && rm -rf $p

  # SWIG
  if $build_swig; then
    p=swig-3.0.7
    echo -e "\\n\\n>> [`date`] building $p"
    test -r $p.tar.gz || wget $wget_opts http://atlas3.atlas.aei.uni-hannover.de/~bema/tarballs/$p.tar.gz
    rm -rf $p
    tar -xzf $p.tar.gz
    cd $p
    ./configure --prefix=$PREFIX --without-tcl --with-python --without-python3 --without-perl5 --without-octave --without-scilab --without-java --without-javascript --without-gcj --without-android --without-guile --without-mzscheme --without-ruby --without-php --without-ocaml --without-pike --without-chicken --without-csharp --without-lua --without-allegrocl --without-clisp --without-r --without-go --without-d 
    make
    make install
    cd ..
    $cleanup && rm -rf $p
  fi

if $build_preinst_before_lalsuite; then
  pushd $PYCBC/..
  tar -czf pycbc-preinst.tgz pycbc
  popd
fi

fi # if pycbc-preinst.tgz

if test -r $PYCBC/../pycbc-preinst-lalsuite.tgz; then
  :
else

  # LALSUITE
  echo -e "\\n\\n>> [`date`] building lalsuite"
  cd lalsuite
  git reset --hard HEAD
  echo -e "\\n\\n>> [`date`] git HEAD: `git log -1 --pretty=oneline --abbrev-commit`"
  sed -i~ s/func__fatal_error/func_fatal_error/ */gnuscripts/ltmain.sh
  if $build_dlls; then
    fgrep -l lib_LTLIBRARIES `find . -name Makefile.am` | while read i; do
      sed -n 's/.*lib_LTLIBRARIES *= *\(.*\).la/\1_la_LDFLAGS += -no-undefined/p' $i >> $i
    done
    sed -i~ 's/\(swiglal_python_la_LDFLAGS = .*\)$/\1 -no-undefined/;
             s/\(swiglal_python_la_LIBADD = .*\)$/\1 -lpython2.7/;
             s/swiglal_python\.la/libswiglal_python.la/g;
             s/swiglal_python_la/libswiglal_python_la/g;
             s/mv -f swiglal_python/mv -f cygswiglal_python/;' gnuscripts/lalsuite_swig.am
    shared="$shared --enable-win32-dll"
  fi
  ./00boot
  cd ..
  rm -rf lalsuite-build
  mkdir lalsuite-build
  cd lalsuite-build
  ../lalsuite/configure --disable-gcc-flags $shared --enable-static --enable-swig-python --prefix="$PREFIX" --disable-silent-rules \
                        --disable-lalxml --disable-lalpulsar --disable-laldetchar --disable-lalstochastic --disable-lalinference
  make
  make install
  for i in $PREFIX/etc/*-user-env.sh; do
    source "$i"
  done
  cd ..
  $cleanup && rm -rf lalsuite-build
  echo -e "\\n\\n>> [`date`] building PyLAL"
  cd lalsuite/pylal
  python setup.py install --prefix="$PREFIX"
#  echo -e "\\n\\n>> [`date`] building GLUE"
#  cd ../glue
#  python setup.py install --prefix="$PREFIX"
  cd ../..
  test -r "$PREFIX/etc/pylal-user-env.sh" && source "$PREFIX/etc/pylal-user-env.sh"
  test -r "$PREFIX/etc/glue-user-env.sh" && source "$PREFIX/etc/glue-user-env.sh"

  pushd $PYCBC/..
  tar -czf pycbc-preinst-lalsuite.tgz pycbc
  popd

fi # if pycbc-preinst.tgz

echo 'Flask==0.10
Flask-Cache==0.13.1
Flask-SQLAlchemy==0.16
Jinja2==2.7
MarkupSafe==0.18
MySQL-python==1.2.5
SQLAlchemy==0.8.0
WTForms==1.0.3
Werkzeug==0.9.3
boto==2.5.2
itsdangerous==0.21
pam==0.1.4
requests==1.2.3
decorator==4.0.4
pycbc-pylal==0.9.5
pyRXP==2.1.0
nose==1.3.7
pkgconfig==1.1.0
six==1.9.0
linecache2==1.0.0
traceback2==1.4.0
unittest2==1.1.0
cycler==0.9.0
python-dateutil==2.4.2
Babel==2.1.1
M2Crypto==0.22.3
Mako==1.0.2
Pillow==2.9.0
Pygments==2.0.2
Sphinx==1.3.1
alabaster==0.7.6
argparse==1.3.0
docutils==0.12
funcsigs==0.4
matplotlib==1.4.3
mock==1.3.0
numpydoc==0.5
pbr==1.8.0
pyparsing==2.0.3
python-cjson==1.1.0
pytz==2015.6
snowballstemmer==1.2.0
sphinx-rtd-theme==0.1.9
sphinxcontrib-programoutput==0.8' > requirements.txt
pip $pip_install -r requirements.txt
# don't downgrade to setuptools==18.2 here yet

if $compile_pycbc_glue; then
  p=pycbc-glue-0.9.6
  echo -e "\\n\\n>> [`date`] building $p"
  test -r $p.tar.gz || wget $wget_opts "http://pypi.python.org/packages/source/p/pycbc-glue/$p.tar.gz"
  rm -rf $p
  tar -xzf $p.tar.gz
  cd $p
  # This is pretty ugly just to get things compiled. It won't work when used.
  # This avoids error: initializer element is not constant
  sed -i~ 's/\.tp_base = &PyTuple_Type/.tp_base = NULL/' src/segments/segment.c
  sed -i~ 's/\.tp_base = &PyList_Type/.tp_base = NULL/' src/segments/segmentlist.c
  python setup.py install --prefix="$PREFIX"
  cd ..
  $cleanup && rm -rf $p
else
  echo -e "\\n\\n>> [`date`] pip install pycbc-glue==0.9.6"
  pip $pip_install pycbc-glue==0.9.6
fi

echo -e "\\n\\n>> [`date`] pip install h5py==2.5.0"
pip $pip_install h5py==2.5.0

# This is a pretty dirty hack faking psycopg2-2.5.5 to be v2.6
# pegasus-wms is pinned to it but I couldn't get 2.6 to compile
# see https://github.com/psycopg/psycopg2/issues/305
# psycopgmodule.c _IS_ compiled with -fPIC, but still doesn't link
# tried with newer gcc-4.4.4 and newer binutils
if $fake_psycopg26; then
  p=psycopg2-2.5.5
  echo -e "\\n\\n>> [`date`] building $p"
  test -r $p.tar.gz || wget $wget_opts "http://pypi.python.org/packages/source/p/psycopg2/$p.tar.gz"
  rm -rf $p
  tar -xzf $p.tar.gz
  cd $p
  sed -i~ 's/2\.5\.5/2.6/' setup.py PKG-INFO
  cd ..
  rm -rf psycopg2-2.6
  mv psycopg2-2.5.5 psycopg2-2.6
  tar -czf psycopg2-2.6f.tar.gz psycopg2-2.6
  cd psycopg2-2.6
# LDFLAGS="-L$PREFIX/lib -fPIC" CPPFLAGS="-I$PREFIX/include" python setup.py build
  python setup.py install --prefix="$PREFIX"
  cd ..
  $cleanup && rm -rf psycopg2-2.6
else
  echo -e "\\n\\n>> [`date`] pip install psycopg2==2.6"
  pip $pip_install psycopg2==2.6
fi

# PEGASUS
v=4.5.2
p=pegasus-source-$v
if test -r $p-lib-pegasus-python.tgz; then
  :
elif $build_pegasus_source; then
  test ".$LC_CTYPE" = ".UTF-8" && export LC_ALL=en_US.UTF-8
  echo -e "\\n\\n>> [`date`] building $p-lib-pegasus-python.tgz"
  test -r $p.tar.gz || wget $wget_opts "http://download.pegasus.isi.edu/pegasus/$v/$p.tar.gz"
  rm -rf $p
  tar -xzf $p.tar.gz
  if $fake_psycopg26; then
    cp psycopg2-2.6f.tar.gz $p/src/externals/psycopg2-2.6.tar.gz
  fi
  cd $p
  ant dist-python-source
  cd ..
  tar -czf $p-lib-pegasus-python.tgz $p/lib/pegasus/python $p/release-tools $p/build.properties
  $cleanup && rm -rf $p
else
  wget $wget_opts http://atlas3.atlas.aei.uni-hannover.de/~bema/tarballs/$p-lib-pegasus-python.tgz
fi
echo -e "\\n\\n>> [`date`] building $p"
tar -xzf $p-lib-pegasus-python.tgz
pushd $p/lib/pegasus/python/
# echo -e "\\n\\n>> [`date`] installing dependencies for $p"
# pip $pip_install -r pegasus_wms.egg-info/requires.txt
python setup.py install --prefix="$PREFIX"
popd
$cleanup && rm -rf $p

# MPLD
p=mpld3-0.3git
# pip $pip_install "https://github.com/ligo-cbc/mpld3/tarball/master#egg=$p"
echo -e "\\n\\n>> [`date`] building $p"
test -r $p.tar.gz || wget $wget_opts -O $p.tar.gz "https://github.com/ligo-cbc/mpld3/tarball/master#egg=$p"
tar -xzf $p.tar.gz
cd ligo-cbc-mpld3-25aee65/
python setup.py install --prefix="$PREFIX"
cd ..
$cleanup && rm -rf ligo-cbc-mpld3-25aee65

# PYINSTALLER
if $build_pyinstaller_git; then
  p=pyinstaller
  echo -e "\\n\\n>> [`date`] building pyinstaller"
  if test -d pyinstaller/.git; then
    cd pyinstaller
  else
    git clone https://github.com/pyinstaller/pyinstaller.git
    cd pyinstaller
    git checkout 3.0
    sed -i~ 's%# Python library NOT found. Return just None.%return "/usr/bin/libpython2.7.dll"%' ./PyInstaller/depend/bindepend.py
  fi
  cd bootloader
  python ./waf distclean all
  cd ..
  python setup.py install --prefix="$PREFIX"
  cd ..
elif $build_pyinstaller; then
  p=PyInstaller-3.0
  echo -e "\\n\\n>> [`date`] building $p"
  test -r $p.tar.gz || wget $wget_opts "https://pypi.python.org/packages/source/P/PyInstaller/$p.tar.gz"
  rm -rf $p
  tar -xzf $p.tar.gz
  cd $p
  # patch PyInstaller to return /usr/bin/libpython2.7.dll when no other Python lib was found
  sed -i~ 's%# Python library NOT found. Return just None.%return "/usr/bin/libpython2.7.dll"%' ./PyInstaller/depend/bindepend.py
  cd bootloader
  python ./waf distclean all
  cd ..
  python setup.py install --prefix="$PREFIX"
  cd ..
  $cleanup && rm -rf $p
else
  echo -e "\\n\\n>> [`date`] pip install pyinstaller==3.0"
  pip $pip_install pyinstaller==3.0
fi

# PYCBC
pip $pip_install --upgrade setuptools==18.2
# NOTE: None of the installation methods listed below works for a local git reposity such that a bundle created with PyInstaller runs successfully,
#       it will always fail with "pycbc distribution not found" (the "pip install" apparently doesn't even install the binaries!?):
#  python setup.py install --prefix="$PREFIX"
#  easy_install .
#  pip install $PWD
#       So currently we're patching libutils after install
if $pip_install_pycbc; then
  echo -e "\\n\\n>> [`date`] pip install pycbc==1.2.0"
  pip $pip_install pycbc==1.2.0
  if $build_dlls; then
    sed -i~ "s/fnmatch.fnmatch(libfile,'lib'+libname+'.dylib\\*')/& or fnmatch.fnmatch(libfile,libname+'.dll') or fnmatch.fnmatch(libfile,'cyg'+libname+'-*.dll')/" \
      $PREFIX/lib/python2.7/site-packages/pycbc/libutils.py
    rm -f $PREFIX/lib/python2.7/site-packages/pycbc/libutils.pyc
  fi
else
  echo -e "\\n\\n>> [`date`] building pycbc"
  if test -d pycbc/.git; then
    cd pycbc
  else
    rm -rf pycbc
    git clone git://github.com/ligo-cbc/pycbc
    cd pycbc
    git checkout v1.2.0
    # make get_libpath_from_dirlist() aware of dlls
    sed -i~ "s/fnmatch.fnmatch(libfile,'lib'+libname+'.dylib\\*')/& or fnmatch.fnmatch(libfile,libname+'.dll') or fnmatch.fnmatch(libfile,'cyg'+libname+'-*.dll')/" pycbc/libutils.py
  fi
  python setup.py install --prefix="$PREFIX"
  cd ..
fi
test -r "$PREFIX/etc/pycbc-user-env.sh" && source "$PREFIX/etc/pycbc-user-env.sh"

# log environment
echo -e "\\n\\n>> [`date`] ENVIRONMENT ..."
env
echo -e "... ENVIRONMENT"

# prepare for Cygwin rebaseall in case the following test fails with a DLL relocation/alignment error
# from https://cygwin.com/ml/cygwin/2009-12/msg00168.html:
# /bin/rebase -d -b 0x61000000 -o 0x20000 -v -T <file with list of dll and so files> > rebase.out
if $build_dlls; then
  find /home/bema/pycbc/environment -name \*.dll > /dlls.txt
  rebase -d -b 0x61000000 -o 0x20000 -v -T /dlls.txt
fi

# TEST
echo -e "\\n\\n>> [`date`] testing local executable"
cd $PREFIX
./bin/pycbc_inspiral --help

# BUNDLE DIR
echo -e "\\n\\n>> [`date`] running pyinstaller"
rm -rf dist/pycbc_inspiral_d dist/pycbc_inspiral
# patch spec file to add "-v" to python interpreter options
pyi-makespec --hidden-import=scipy.linalg.cython_blas --hidden-import=scipy.linalg.cython_lapack --hidden-import=pkg_resources --onedir ./bin/pycbc_inspiral
if $add_python_runtime_options; then
  sed -i~ 's%exe = EXE(pyz,%options = [ ("v", None, "OPTION"), ("W error", None, "OPTION") ]\
exe = EXE(pyz, options,%' pycbc_inspiral.spec
fi
pyinstaller pycbc_inspiral.spec

# TEST BUNDLE
echo -e "\\n\\n>> [`date`] testing"
cd dist/pycbc_inspiral
./pycbc_inspiral --help

# end here for Windows
if $build_dlls; then
  exit
fi

# BUNDLE FILE
cd ..
mv pycbc_inspiral pycbc_inspiral_d
cd ..
echo -e "\\n\\n>> [`date`] running pyinstaller"
pyinstaller --hidden-import=scipy.linalg.cython_blas --hidden-import=scipy.linalg.cython_lapack --hidden-import=pkg_resources --onefile ./bin/pycbc_inspiral

# TEST BUNDLE
echo -e "\\n\\n>> [`date`] testing"

# leave virtual environment
deactivate

./dist/pycbc_inspiral --help

echo -e "\\n\\n>> [`date`] Success $0"
