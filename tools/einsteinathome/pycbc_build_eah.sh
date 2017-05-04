#!/bin/bash

# script to build PyCBC suite on Debian Etch and Windows Cygwin

# FIXME/todo:
#

# make a pipe fail if any process fails
set -o pipefail

# check an md5 sum of a file, independent of platform
# $1 file to check, $2 md5 checksum
# returns 0 for fail, 1 for success(!)
check_md5() {
    md5s=`( md5sum $1 2>/dev/null || md5 $1 2>/dev/null) | sed 's/.* = //;s/  .*//;s/^\(................................\).*/\1/'`
    test ".$2" != ".$md5s"
}

function exit_on_error {
  rm -f "$PYCBC/lock"
  if $silent_build ; then
      if [ -f $LOG_FILE ] ; then
          echo "--- Error or interrupt: dumping last 100 lines of log ---------" >&4
          tail -n 100 $LOG_FILE >&4
          echo "---------------------------------------------------------------" >&4
          rm -f $LOG_FILE
      fi
  fi
  exit 1
}
trap exit_on_error ERR INT

echo -e ">> [`date`] Start $0 $*"

test ".$LANG" = "." && export LANG="en_US.UTF-8"
test ".$LC_ALL" = "." && export LC_ALL="$LANG"
export CC=gcc
export CXX=g++
export FC=gfortran

# compilation environment
BUILDDIRNAME="pycbc-build"

# defaults, possibly overwritten by command-line arguments
cleanup=true # usually, build directories are removed after a successful build
verbose_pyinstalled_python=false
pycbc_branch=master
pycbc_remote=ligo-cbc
pycbc_fetch_ref=""
scratch_pycbc=false
libgfortran=libgfortran.so
extra_libs=""
extra_bank=""
extra_approx=""
processing_scheme=""
lal_data_path="."

# defaults, possibly overwritten by OS-specific settings
build_ssl=false
build_python=false
fftw_flags=--enable-avx
shared="--enable-shared"
static="--disable-static"
build_dlls=false
rebase_dlls_before_pycbc=false
build_lapack=true
pyssl_from="tarball" # "pip-install"
numpy_from="pip-install" # "tarball"
scipy_from="pip-install" # "git"
scipy_version=0.19.0
build_gsl=true
build_swig=true
build_pcre=false
build_fftw=true
build_framecpp=false
build_preinst_before_lalsuite=true
build_subprocess32=false
build_hdf5=true
build_freetype=true
build_zlib=true
build_pegasus=true
build_wrapper=false
build_fstab=true
pyinstaller_version=v3.2.1 # 9d0e0ad4, v2.1 .. v3.2.1 -> git, 2.1 .. 3.2.1 -> pypi
patch_pyinstaller_bootloader=true
pyinstaller21_hacks=false # use hacks & workarounds necessary for PyInstaller <3.0
use_pycbc_pyinstaller_hooks=true
build_onefile_bundles=false
run_analysis=true
silent_build=false

if echo ".$WORKSPACE" | grep CYGWIN64_FRONTEND >/dev/null; then
    # hack to use the script as a frontend for a Cygwin build slave for a Jenkins job
    # example: WORKSPACE='/Users/jenkins/workspace/workspace/EAH_PyCBC_Master/label/OSX107'
    test ".$CYGWIN_HOST" = "."      && CYGWIN_HOST=moss-cygwin64
    test ".$CYGWIN_HOST_USER" = "." && CYGWIN_HOST_USER=jenkins
    test ".$CYGWIN_HOST_PORT" = "." && CYGWIN_HOST_PORT=2222
    echo -e "\\n\\n>> [`date`] running remotely at $CYGWIN_HOST_USER@$CYGWIN_HOST:$CYGWIN_HOST_PORT"
    echo "WORKSPACE='$WORKSPACE'" # for Jenkins jobs
    unset WORKSPACE # avoid endless recoursion
    # copy the script
    scp "-P$CYGWIN_HOST_PORT" "$0" "$CYGWIN_HOST_USER@$CYGWIN_HOST:."
    # run it remotely
    ssh "-p$CYGWIN_HOST_PORT" "$CYGWIN_HOST_USER@$CYGWIN_HOST" bash -l `basename $0` "$@"
    # fetch the artifacts to local workspace
    dist="$BUILDDIRNAME/environment/dist"
    rm -rf "$dist"
    mkdir -p "$dist"
    scp "-P$CYGWIN_HOST_PORT" "$CYGWIN_HOST_USER@$CYGWIN_HOST:$dist/*.exe" "$dist"
    scp "-P$CYGWIN_HOST_PORT" "$CYGWIN_HOST_USER@$CYGWIN_HOST:$dist/*.zip" "$dist"
    exit 0
elif test ".$1" = ".--force-debian4" ||
  test "v`cat /etc/debian_version 2>/dev/null`" = "v4.0" ||
  lsb_release -a 2>/dev/null | grep 'Ubuntu 6.06' >/dev/null; then
    # detected a Debian 4.0 (Etch) or Ubuntu 6.06 installation, expect a prepared build system
    echo -e "\\n\\n>> [`date`] Using Debian 4.0 (etch) settings"
    test ".$LC_ALL" = "." && export LC_ALL="$LANG"
    link_gcc_version=4.8.5
    gcc_path="/usr/local/bin"
    build_ssl=true
    pyinstaller_lsb="--no-lsb"
    $pyinstaller21_hacks || build_subprocess32=true
    build_onefile_bundles=true
    appendix="_Linux64"
elif grep -q "Scientific Linux release 6" /etc/redhat-release 2>/dev/null; then # SL6
    echo -e "\\n\\n>> [`date`] Using Scientific Linux release 6 (Carbon) settings"
    test ".$LC_ALL" = "." && export LC_ALL="$LANG"
    link_gcc_version=4.4.7
    gcc_path="/usr/bin"
    build_ssl=true
    build_python=true
    build_pegasus=false
    pyinstaller_lsb="--no-lsb"
    build_onefile_bundles=true
    appendix="_Linux64"
elif grep -q "Scientific Linux CERN SLC release 6" /etc/redhat-release 2>/dev/null; then # SL6
    echo -e "\\n\\n>> [`date`] Using Scientific Linux CERN SLC release 6 (Carbon) settings"
    test ".$LC_ALL" = "." && export LC_ALL="$LANG"
    link_gcc_version=4.4.7
    gcc_path="/usr/bin"
    build_python=true
    build_hdf5=true
    build_pegasus=false
    build_fftw=false
    build_gsl=false
    build_ssl=false
    build_lapack=false
    build_freetype=false
    build_zlib=false
    build_wrapper=false
    build_fstab=false
    pyinstaller_lsb="--no-lsb"
    build_onefile_bundles=true
    appendix="_Linux64"
elif grep -q "Ubuntu 12" /etc/issue 2>/dev/null; then
    link_gcc_version=4.6
    gcc_path="/usr/bin"
    build_python=true
    build_pcre=true
    pyinstaller_lsb="--no-lsb"
    build_onefile_bundles=false
    appendix="_Linux64"
    if test x$TRAVIS_OS_NAME = xlinux ; then
        build_fftw=false
        build_hdf5=false
        build_ssl=false
        build_lapack=false
        build_gsl=false
        build_freetype=false
        build_zlib=false
        build_wrapper=false
        build_fstab=false
    fi
elif test "`uname -s`" = "Darwin" ; then # OSX
    echo -e "\\n\\n>> [`date`] Using OSX 10.7 settings"
    export FC=gfortran-mp-4.8
    export CC=gcc-mp-4.8
    export CXX=g++-mp-4.8
    export FFLAGS="$FFLAGS -m64"
    export CFLAGS="$CFLAGS -m64"
    export CXXFLAGS="$CXXFLAGS -m64"
    export LDFLAGS="$LDFLAGS -m64"
#    libframe_debug_level=3
    gcc_path="/opt/local/bin"
    libgfortran=libgfortran.dylib
    fftw_cflags="-Wa,-q"
    build_framecpp=true
    appendix="_OSX64"
elif uname -s | grep ^CYGWIN >/dev/null; then # Cygwin (Windows)
    echo -e "\\n\\n>> [`date`] Using Cygwin settings"
    lal_cppflags="-D_WIN32"
    build_dlls=true
    rebase_dlls_before_pycbc=false
    numpy_from="tarball" # "pip-install"
    build_hdf5=false
    build_freetype=false
    build_gsl=false
    build_swig=false
#    pyinstaller_version=develop # HEAD # 9d0e0ad4
#    pyinstaller21_hacks=true
    patch_pyinstaller_bootloader=false
    appendix="_Windows64"
else
    echo ERROR: unknown OS, exiting
    exit 1
fi

# compilation environment
PYCBC="$PWD/$BUILDDIRNAME"
SOURCE="$PWD/pycbc-sources"
PYTHON_PREFIX="$PYCBC"
ENVIRONMENT="$PYCBC/environment"
PREFIX="$ENVIRONMENT"
PATH="$PREFIX/bin:$PYTHON_PREFIX/bin:$PATH:/opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin"

# log environment
if [ ".$1" == ".--print-env" ]; then
    echo -e "\\n\\n>> [`date`] ENVIRONMENT ..."
    env
    echo -e "... ENVIRONMENT"
fi

# locking
if [ -r "$PYCBC/lock" ]; then
    for pid in `cat "$PYCBC/lock"`; do
        while ps -p "$pid" >/dev/null; do
            echo -e ">> [`date`] waiting for build with PID $pid to finish"
            sleep 30
        done
    done
else
    mkdir -p "$PYCBC"
fi
echo "$$" > "$PYCBC/lock"

# make sure the sources directory exixts
mkdir -p "$SOURCE"

# usage message
usage="
    --force-debian4                 force Debian 4.0 build, must be first option in command-line

    --print-env                     dump environment at beginning, must be first option in command-line

    --help                          print this messge and exit

    --clean                         perform a clean build (takes quite some time); delete ~/.cache and
                                    tarballs containing precompiled libraries (lalsuite, scipy etc.)

    --clean-lalsuite                checkout and build lalsuite from scratch

    --clean-sundays                 perform a clean-lalsuite build on sundays

    --clean-pycbc                   check out pycbc git repo from scratch

    --clean-weave-cache             clean weave code cache before running analysis

    --lalsuite-commit=<commit>      specify a commit (or tag or branch) of lalsuite to build from

    --blessed-lalsuite              get the lalsuite commit to build from
                                    https://github.com/ligo-cbc/pycbc/releases/latest

    --pycbc-commit=<commit>         specify a commit or tag of pycbc to build from (specifying a
                                    branch will only work reliably in conjunction with --clean-pycbc)

    --pycbc-remote=<username>       add pycbc repository github.com/username as remote

    --pycbc-branch=<branch>         checkout branch before building

    --no-pycbc-update               don't update local pycbc repo

    --no-lalsuite-update            don't update local lalsuite repo

    --bema-testing                  use einsteinathome_testing branch of bema-ligo/pycbc repo

    --no-cleanup                    keep build directories after successful build for later inspection

    --with-extra-libs=<url>         add extra files from a tar file at <url> to the bundles

    --with-extra-bank=<file>        run pycbc_inspiral again with an extra template bank

    --with-extra-approximant=<file> run pycbc_inspiral again with an extra approximant

    --with-lal-data-path=<path>     run test job using ROM data from <path>

    --processing-scheme=<scheme>    run test job using processing scheme <scheme>

    --verbose-python                run PyInstalled Python in verbose mode, showing imports

    --no-analysis                   for testing, don't run analysis, assume weave cache is already there

    --silent-build                  do not brint build messages unless there is an error

    --pycbc-fetch-ref               fetch and use a specific reference for pycbc
"

# circumvent old certificate chains on old systems
wget_opts="-c --passive-ftp --no-check-certificate --tries=5 --timeout=30"
export GIT_SSL_NO_VERIFY=true
export PIP_TRUSTED_HOST="pypi.python.org github.com"

# handle command-line arguments, possibly overriding above settings
for i in $*; do
    case $i in
        --force-debian4) ;;
        --print-env) ;;
        --no-pycbc-update) pycbc_branch="HEAD";;
        --no-lalsuite-update) no_lalsuite_update=true;;
        --bema-testing)
            pycbc_branch=einsteinathome_testing
            pycbc_remote=bema-ligo;;
        --pycbc-remote=*) pycbc_remote="`echo $i|sed 's/^--pycbc-remote=//'`";;
        --pycbc-branch=*) pycbc_branch="`echo $i|sed 's/^--pycbc-branch=//'`";;
        --no-cleanup) cleanup=false;;
        --no-analysis) run_analysis=false;;
        --verbose-python) verbose_pyinstalled_python=true;;
        --clean) rm -rf "$HOME/.cache" "$HOME/Library/Caches/pip" "$SOURCE/$BUILDDIRNAME-preinst.tgz" "$SOURCE/$BUILDDIRNAME-preinst-lalsuite.tgz" "$PYCBC";;
        --clean-lalsuite) rm -rf "$SOURCE/lalsuite" "$SOURCE/$BUILDDIRNAME-preinst-lalsuite.tgz";;
        --lalsuite-commit=*) lalsuite_branch="`echo $i|sed 's/^--lalsuite-commit=//'`";;
        --blessed-lalsuite) lalsuite_branch=`wget $wget_opts -O- https://github.com/ligo-cbc/pycbc/releases/latest |
            awk -F'>' '(p==0) && /This release has been tested against LALSuite with the hash:/ {p=1}; (p==1) && /<code>[0-9a-f]*$/{print $NF; p=2}'`;
            if ! echo "0$lalsuite_branch"|egrep '^[0-9a-f]*$' >/dev/null; then
                echo "ERROR: blessed LALSuite version couldn't be determined" >&2
                exit 1;
            fi;;
        --pycbc-commit=*) pycbc_commit="`echo $i|sed 's/^--pycbc-commit=//'`";;
        --pycbc-fetch-ref=*) pycbc_fetch_ref="`echo $i|sed 's/^--pycbc-fetch-ref=//'`";;
        --clean-pycbc) scratch_pycbc=true;;
        --clean-weave-cache) rm -rf "$SOURCE/test/pycbc_inspiral";;
        --clean-sundays)
            if [ `date +%u` -eq 7 ]; then
                if [ -r "$SOURCE/last_sunday_build" ]; then
                    last_build=`cat "$SOURCE/last_sunday_build"`
                else
                    last_build=0
                fi
                now=`date +%s`
                d2_ago=`expr $now - 172800` # two days ago
                if [  $last_build -le $d2_ago ]; then # last 'clean-sundays' build was two days ago or older
                    rm -rf "$HOME/.cache" "$SOURCE/$BUILDDIRNAME-preinst-lalsuite.tgz" "$SOURCE/lalsuite/configure"
                    echo $now > "$SOURCE/last_sunday_build"
                fi
            fi ;;
        --with-extra-libs=*) extra_libs="`echo $i|sed 's/^--with-extra-libs=//'`";;
        --with-extra-bank=*) extra_bank="$extra_bank `echo $i|sed 's/^--with-extra-bank=//'`";;
        --with-extra-approximant=*) extra_approx="${extra_approx}`echo $i|sed 's/^--with-extra-approximant=//'` ";;
        --with-lal-data-path=*) lal_data_path="`echo $i|sed 's/^--with-lal-data-path=//'`";;
        --processing-scheme=*) processing_scheme="--processing-scheme `echo $i|sed 's/^--processing-scheme=//'`";;
        --silent-build) silent_build=true;;
        --help) echo -e "Options:\n$usage">&2; exit 0;;
        *) echo -e "unknown option '$i', valid are:\n$usage">&2; exit 1;;
    esac
done


# compilation environment
if [ ".$link_gcc_version" != "." ]; then
    mkdir -p $PYTHON_PREFIX/bin
    ( cd $PYTHON_PREFIX/bin &&
        for i in gcc g++ gfortran; do
            rm -f $i &&
                if test -x "${gcc_path}/$i-$link_gcc_version"; then
                    ln -s "${gcc_path}/$i-$link_gcc_version" $i
                else
                    echo ERROR; "${gcc_path}/$i-$link_gcc_version" not found
                fi
        done
    )
fi
libgfortran_dir="`$FC -print-file-name=$libgfortran|sed 's%/[^/]*$%%'`"
export LD_LIBRARY_PATH="$PREFIX/lib:$PREFIX/bin:$PYTHON_PREFIX/lib:$libgfortran_dir:/usr/local/lib:$LD_LIBRARY_PATH"
export CPPFLAGS="-I$PREFIX/include -I$PYTHON_PREFIX/include $CPPFLAGS"
export PKG_CONFIG_PATH="$PREFIX/lib/pkgconfig:$PYTHON_PREFIX/lib/pkgconfig:/usr/local/lib/pkgconfig:$PKG_CONFIG_PATH"
export LIBS="$LIBS -lgfortran"

# log compilation environment
echo "export PATH='$PATH'"
echo "export LD_LIBRARY_PATH='$LD_LIBRARY_PATH'"
echo "export PKG_CONFIG_PATH='$PKG_CONFIG_PATH'"
echo "WORKSPACE='$WORKSPACE'" # for Jenkins jobs

# obsolete, kept for reference
#export CPPFLAGS="-I$PREFIX/include"
#export LDFLAGS="-L$PREFIX/lib -L$libgfortran"
#export LDFLAGS="-L$libgfortran $LDFLAGS"
# -static-libgfortran

# URL abbrevations
pypi="https://pypi.python.org/packages"
gitlab="https://gitlab.aei.uni-hannover.de/einsteinathome"
atlas="https://www.atlas.aei.uni-hannover.de/~bema"
albert="http://albert.phys.uwm.edu/download"
duncan="https://www.atlas.aei.uni-hannover.de/~dbrown"
aei="http://www.aei.mpg.de/~bema"

if $silent_build ; then
    LOG_FILE=$(mktemp -t pycbc-build-log.XXXXXXXXXX)
    echo -e "\\n\\n>> [`date`] writing build logs to $LOG_FILE"

    # make a copy of stdin and stdout and close them
    exec 3>&1-
    exec 4>&2-

    # open stdout as $LOG_FILE file for read and write.
    exec 1<>$LOG_FILE

    # redirect stderr to stdout
    exec 2>&1
else
    # make a copy of stdin and stdout
    exec 3>&1
    exec 4>&2
fi

# use previously compiled scipy, lalsuite etc. if available
if test -r "$SOURCE/$BUILDDIRNAME-preinst.tgz" -o -r "$SOURCE/$BUILDDIRNAME-preinst-lalsuite.tgz"; then

    rm -rf "$PYCBC"
    if test -r "$SOURCE/$BUILDDIRNAME-preinst-lalsuite.tgz"; then
        echo -e "\\n\\n>> [`date`] using $BUILDDIRNAME-preinst-lalsuite.tgz" >&3
        tar -xzf "$SOURCE/$BUILDDIRNAME-preinst-lalsuite.tgz"
    else
        echo -e "\\n\\n>> [`date`] using $BUILDDIRNAME-preinst.tgz" >&3
        tar -xzf "$SOURCE/$BUILDDIRNAME-preinst.tgz"
    fi
    # set up virtual environment
    unset PYTHONPATH
    source "$ENVIRONMENT/bin/activate"
    echo -e "\\n\\n>> [`date`] pip install --upgrade pip" >&3
    pip install --upgrade pip
    # workaround to make the virtualenv accept .pth files
    export PYTHONPATH="$PREFIX/lib/python2.7/site-packages:$PYTHONPATH"
    cd "$SOURCE"
    if ls $PREFIX/etc/*-user-env.sh >/dev/null 2>&1; then
        for i in $PREFIX/etc/*-user-env.sh; do
            source "$i"
        done
    fi

else # if $BUILDDIRNAME-preinst.tgz

    cd "$SOURCE"

    # OpenSSL
    if $build_ssl; then
    # p=openssl-1.0.2e # compile error on pyOpenSSL 0.13:
    # pycbc/include/openssl/x509.h:751: note: previous declaration X509_REVOKED_ was here
	p=openssl-1.0.1p
	echo -e "\\n\\n>> [`date`] building $p" >&3
	test -r $p.tar.gz || wget $wget_opts $aei/$p.tar.gz
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
	echo -e "\\n\\n>> [`date`] building $p" >&3
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
	echo -e "\\n\\n>> [`date`] pip install --upgrade pip" >&3
	pip install --upgrade pip
	echo -e "\\n\\n>> [`date`] pip install virtualenv" >&3
	pip install virtualenv
    fi

    # set up virtual environment
    unset PYTHONPATH
    rm -rf "$ENVIRONMENT"
    virtualenv "$ENVIRONMENT"
    source "$ENVIRONMENT/bin/activate"
    # workaround to make the virtualenv accept .pth files
    export PYTHONPATH="$PREFIX/lib/python2.7/site-packages:$PYTHONPATH"

    # pyOpenSSL-0.13
    if [ "$pyssl_from" = "pip-install" ] ; then
	echo -e "\\n\\n>> [`date`] pip install pyOpenSSL==0.13" >&3
	pip install pyOpenSSL==0.13
    else
	p=pyOpenSSL-0.13
	echo -e "\\n\\n>> [`date`] building $p" >&3
	test -r $p.tar.gz || wget $wget_opts "$pypi/source/p/pyOpenSSL/$p.tar.gz"
	rm -rf $p
	tar -xzf $p.tar.gz
	cd $p
	sed -i~ 's/X509_REVOKED_dup/X509_REVOKED_dup_static/' OpenSSL/crypto/crl.c
	python setup.py build_ext "-I$PYTHON_PREFIX/include" "-L$PYTHON_PREFIX/lib"
	python setup.py build
	python setup.py install --prefix="$PREFIX"
	cd ..
	$cleanup && rm -rf $p
    fi

    # LAPACK & BLAS
    if $build_lapack; then
	p=lapack-3.6.0
	echo -e "\\n\\n>> [`date`] building $p" >&3
	test -r $p.tgz || wget $wget_opts http://www.netlib.org/lapack/$p.tgz
	rm -rf $p
	tar -xzf $p.tgz
	cd $p
        # configure: compile with -fPIC, remove -frecoursive, build deprecated functions
	sed "s/ *-frecursive//;s/gfortran/$FC -fPIC -m64/;s/^#MAKEDEPRECATED.*/BUILD_DEPRECATED = Yes/" make.inc.example > make.inc
	make lapack_install lib blaslib
	mkdir -p "$PREFIX/lib"
	cp lib*.a "$PREFIX/lib"
	cp librefblas.a "$PREFIX/lib/libblas.a"
	cd ..
	$cleanup && rm -rf $p
    fi
    
    # NUMPY
    if [ "$numpy_from" = "tarball" ]; then
	p=numpy-1.9.3
	echo -e "\\n\\n>> [`date`] building $p" >&3
	test -r $p.tar.gz || wget $wget_opts $pypi/source/n/numpy/$p.tar.gz
	rm -rf $p
	tar -xzf $p.tar.gz 
	cd $p
	python setup.py build --fcompiler=$FC
	python setup.py install --prefix=$PREFIX
	cd ..
	$cleanup && rm -rf $p
    else
	echo -e "\\n\\n>> [`date`] pip install numpy==1.9.3" >&3
	pip install numpy==1.9.3
    fi

    echo -e "\\n\\n>> [`date`] pip install nose, Cython-0.23.2" >&3
    pip install nose Cython==0.23.2

    echo -e "\\n\\n>> [`date`] make sure dbhash and shelve exist" >&3
    python -c "import dbhash, shelve"

    # SCIPY
    if [ "$scipy_from" = "pip-install" ] ; then
	echo -e "\\n\\n>> [`date`] pip install scipy==$scipy_version" >&3
	pip install scipy==$scipy_version
    else
        echo -e "\\n\\n>> [`date`] building scipy-$scipy_version from git" >&3
        if test -d scipy/.git; then
            rm -rf scipy
        fi
        git clone https://github.com/scipy/scipy.git
        cd scipy
        git checkout v$scipy_version
        if test "v$scipy_version" = "v0.16.0"; then
            git config user.name "Dummy"
            git config user.email "dummy@dummy.net"
            git cherry-pick 832baa20f0b5d521bcdf4784dda13695b44bb89f
        fi
        python setup.py build --fcompiler=$FC
        python setup.py install --prefix=$PREFIX
        cd ..
    fi

    # this test will catch scipy build errors that would not emerge before running pycbc_inspiral
    echo -e "\\n\\n>> [`date`] Testing: python -c 'from scipy.io.wavfile import write as write_wav'" >&3
    python -c 'from scipy.io.wavfile import write as write_wav'
    # python -c 'import scipy; scipy.test(verbose=2);'

    # GSL
    if $build_gsl; then
	p=gsl-1.16
	echo -e "\\n\\n>> [`date`] building $p" >&3
	test -r $p.tar.gz || wget $wget_opts ftp://ftp.fu-berlin.de/unix/gnu/gsl/$p.tar.gz
	rm -rf $p
	tar -xzf $p.tar.gz
	cd $p
	./configure $shared $static --prefix="$PREFIX"
	make
	make install
	cd ..
	$cleanup && rm -rf $p
    fi

    # FFTW
    if $build_fftw ; then
        p=fftw-3.3.5
        echo -e "\\n\\n>> [`date`] building $p" >&3
        test -r $p.tar.gz ||
            wget $wget_opts $aei/$p.tar.gz ||
            wget $wget_opts ftp://ftp.fftw.org/pub/fftw/$p.tar.gz
        rm -rf $p
        tar -xzf $p.tar.gz
        cd $p
        if test ".$fftw_cflags" = "."; then
            ./configure $shared $static --prefix="$PREFIX" --enable-sse2 $fftw_flags
        else
            ./configure CFLAGS="$CFLAGS $fftw_cflags" $shared $static --prefix="$PREFIX" --enable-sse2 $fftw_flags
        fi
        make
        make install
        if test ".$fftw_cflags" = "."; then
            ./configure $shared $static --prefix="$PREFIX" --enable-float --enable-sse $fftw_flags
        else
            ./configure CFLAGS="$CFLAGS $fftw_cflags" $shared $static --prefix="$PREFIX" --enable-float --enable-sse $fftw_flags
        fi
        make
        make install
        cd ..
        $cleanup && rm -rf $p
    fi

    # ZLIB
    if $build_zlib ; then
        p=zlib-1.2.8
        echo -e "\\n\\n>> [`date`] building $p" >&3
        test -r $p.tar.gz || wget $wget_opts $aei/$p.tar.gz
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
    fi

    # HDF5
    if $build_hdf5; then
	p=hdf5-1.8.13
	echo -e "\\n\\n>> [`date`] building $p" >&3
	test -r $p.tar.gz ||
            wget $wget_opts $aei/$p.tar.gz ||
            wget $wget_opts https://support.hdfgroup.org/ftp/HDF5/releases/$p/src/$p.tar.gz
	rm -rf $p
	tar -xzf $p.tar.gz
	cd $p
	./configure $shared $static --prefix="$PREFIX"
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
	echo -e "\\n\\n>> [`date`] building $p" >&3
	test -r $p.tar.gz || wget $wget_opts http://download.savannah.gnu.org/releases/freetype/freetype-old/$p.tar.gz
	rm -rf $p
	tar -xzf $p.tar.gz
	cd $p
	./configure $shared $static --prefix="$PREFIX"
	make
	make install
	cd ..
	$cleanup && rm -rf $p
    fi

    echo -e "\\n\\n>> [`date`] pip install --upgrade distribute" >&3
    pip install --upgrade distribute

    if $build_framecpp; then

        # FrameCPP
        p=ldas-tools-al-2.5.6
        echo -e "\\n\\n>> [`date`] building $p" >&3
        test -r $p.tar.gz || wget $wget_opts http://software.ligo.org/lscsoft/source/$p.tar.gz
        rm -rf $p
        tar -xzf $p.tar.gz
        cd $p
        ./configure --disable-latex --disable-swig --disable-python --disable-tcl --enable-64bit --disable-warnings-as-errors $shared $static --prefix="$PREFIX" # --disable-cxx11
        make
        make install
        cd ..
        $cleanup && rm -rf $p
        p=ldas-tools-framecpp-2.5.5
        test -r $p.tar.gz || wget $wget_opts http://software.ligo.org/lscsoft/source/$p.tar.gz
        rm -rf $p
        tar -xzf $p.tar.gz
        cd $p
        ./configure --disable-latex --disable-swig --disable-python --disable-tcl --enable-64bit --disable-warnings-as-errors $shared $static --prefix="$PREFIX" # --disable-cxx11
        make
        make install
        cd ..
        $cleanup && rm -rf $p

    fi # else # build_framecpp

        # LIBFRAME / FrameL
        p=libframe-8.30
        echo -e "\\n\\n>> [`date`] building $p" >&3
        test -r $p.tar.gz || wget $wget_opts http://lappweb.in2p3.fr/virgo/FrameL/$p.tar.gz
        rm -rf $p
        tar -xzf $p.tar.gz
        cd $p
        if test -n "$libframe_debug_level"; then
            sed -i~ "s/^FILE *\\*FrFOut *= *NULL;/#define FrFOut stderr/;
                s/FrDebugLvl *= *[^;]*;/FrDebugLvl = $libframe_debug_level;/" src/FrameL.c
        fi
        # make sure frame files are opened in binary mode
        sed -i~ 's/\([Oo]pen.*"r\)"/\1b"/;' src/FrameL.c # `egrep -r '[Oo]pen.*"r"' .`
        if $build_dlls; then
            for i in src/Makefile*; do
                echo 'libFrame_la_LDFLAGS += -no-undefined' >> $i
            done
        fi
        ./configure $shared $static --prefix="$PREFIX"
        make
        make install
        mkdir -p "$PREFIX/lib/pkgconfig"
        sed "s%^prefix=.*%prefix=$PREFIX%" src/libframe.pc > $PREFIX/lib/pkgconfig/libframe.pc
        cd ..
        $cleanup && rm -rf $p

    # fi # build_framecpp

    # METAIO
    p=metaio-8.3.0
    echo -e "\\n\\n>> [`date`] building $p" >&3
    test -r $p.tar.gz || wget $wget_opts http://software.ligo.org/lscsoft/source/$p.tar.gz
    rm -rf $p
    tar -xzf $p.tar.gz
    cd $p
    if $build_dlls; then
	for i in src/Makefile*; do
	    echo 'libmetaio_la_LDFLAGS += -no-undefined' >> $i
	done
    fi
    ./configure $shared $static --prefix="$PREFIX"
    make
    make install
    cd ..
    $cleanup && rm -rf $p

    # SWIG
    if $build_swig; then
	p=swig-3.0.7
	echo -e "\\n\\n>> [`date`] building $p" >&3
	test -r $p.tar.gz ||
            wget $wget_opts "$aei/$p.tar.gz" ||
            wget $wget_opts "$atlas/tarballs/$p.tar.gz"
	rm -rf $p
	tar -xzf $p.tar.gz
	cd $p
        if $build_pcre; then
            test -r pcre-8.40.tar.gz || wget $wget_opts "https://ftp.pcre.org/pub/pcre/pcre-8.40.tar.gz"
            ./Tools/pcre-build.sh
        fi
	./configure --prefix=$PREFIX --without-tcl --with-python --without-python3 --without-perl5 --without-octave \
            --without-scilab --without-java --without-javascript --without-gcj --without-android --without-guile \
            --without-mzscheme --without-ruby --without-php --without-ocaml --without-pike --without-chicken \
            --without-csharp --without-lua --without-allegrocl --without-clisp --without-r --without-go --without-d 
	make
	make install
	cd ..
	$cleanup && rm -rf $p
    fi

    if $build_preinst_before_lalsuite; then
        echo -e "\\n\\n>> [`date`] building $BUILDDIRNAME-preinst.tgz"
        pushd $PYCBC/..
        tar -czf "$SOURCE/$BUILDDIRNAME-preinst.tgz" "$BUILDDIRNAME"
        popd
    fi

fi # if $BUILDDIRNAME-preinst.tgz

if test -r "$SOURCE/$BUILDDIRNAME-preinst-lalsuite.tgz"; then
    :
else

    # LALSUITE
    if [ ".$no_lalsuite_update" != "." ]; then
        echo -e "\\n\\n>> [`date`] Not updating lalsuite" >&3
	cd lalsuite
    elif test -d lalsuite/.git; then
        echo -e "\\n\\n>> [`date`] Updating lalsuite" >&3
        cd lalsuite
        git reset --hard HEAD
        if [ ".$lalsuite_branch" != ".HEAD" ]; then
            git checkout master
            git pull
        fi
        if [ ".$lalsuite_branch" = ".eah_cbc" ]; then
            git remote update gitlab
            git branch -D eah_cbc
            git checkout -b eah_cbc gitlab/eah_cbc
        elif [ ".$lalsuite_branch" != "." ]; then
            git checkout "$lalsuite_branch"
        fi
    else
        echo -e "\\n\\n>> [`date`] Cloning lalsuite" >&3
        git clone git://versions.ligo.org/lalsuite.git
        cd lalsuite
        git remote add gitlab $gitlab/lalsuite.git
        if [ ".$lalsuite_branch" != "." ]; then
            git checkout "$lalsuite_branch"
        fi
    fi
    echo -e ">> [`date`] git HEAD: `git log -1 --pretty=oneline --abbrev-commit`" >&3
    # fix typo in ltmain
    # only helpful if something goes wrong, normally not needed
    # sed -i~ s/func__fatal_error/func_fatal_error/ */gnuscripts/ltmain.sh
    if $build_dlls; then
	git apply <<'EOF' || true
From accb37091abbc8d8776edfb3484259f6059c4e25 Mon Sep 17 00:00:00 2001
From: Karl Wette <karl.wette@ligo.org>
Date: Mon, 6 Mar 2017 20:18:07 +0100
Subject: [PATCH] SWIG: add Python libraries to linker flags

- Some platforms require them, e.g. cygwin
---
 gnuscripts/lalsuite_swig.m4 | 4 +++-
 1 file changed, 3 insertions(+), 1 deletion(-)

diff --git a/gnuscripts/lalsuite_swig.m4 b/gnuscripts/lalsuite_swig.m4
index 831ff4a..2a2695e 100644
--- a/gnuscripts/lalsuite_swig.m4
+++ b/gnuscripts/lalsuite_swig.m4
@@ -493,6 +493,8 @@ sys.stdout.write(' -L' + cfg.get_python_lib())
 sys.stdout.write(' -L' + cfg.get_python_lib(plat_specific=1))
 sys.stdout.write(' -L' + cfg.get_python_lib(plat_specific=1,standard_lib=1))
 sys.stdout.write(' -L' + cfg.get_config_var('LIBDIR'))
+sys.stdout.write(' -lpython%i.%i' % (sys.version_info.major, sys.version_info.minor))
+sys.stdout.write(' ' + cfg.get_config_var('LIBS'))
 EOD`]
     AS_IF([test $? -ne 0],[
       AC_MSG_ERROR([could not determine Python linker flags])
-- 
2.7.4

EOF
	fgrep -l lib_LTLIBRARIES `find . -name Makefile.am` | while read i; do
	    sed -n 's/.*lib_LTLIBRARIES *= *\(.*\).la/\1_la_LDFLAGS += -no-undefined/p' $i >> $i
	done
	sed -i~ 's/\(swiglal_python_la_LDFLAGS = .*\)$/\1 -no-undefined/' gnuscripts/lalsuite_swig.am
        rm -f configure # force running ./00boot if lalsuite was patched before 'make' triggers it
        shared="$shared --enable-win32-dll"
    fi
    # if we are using FrameCPP, enable it (and disable FrameLib) for LALFrame
    if $build_framecpp; then
        shared="$shared --enable-framec --disable-framel"
    fi
    if [ ! -x configure ]; then
        git clean -xdff
        echo -e "\\n\\n>> [`date`] Creating lalsuite configure scripts" >&3
        ./00boot
    fi
    cd ..
    rm -rf lalsuite-build
    mkdir lalsuite-build
    cd lalsuite-build
    echo -e "\\n\\n>> [`date`] Configuring lalsuite" >&3
    ../lalsuite/configure CPPFLAGS="$lal_cppflags $CPPFLAGS" --disable-gcc-flags $shared $static --prefix="$PREFIX" --disable-silent-rules \
        --disable-all-lal --enable-lalframe --enable-lalmetaio --enable-lalsimulation --enable-lalinspiral --enable-lalpulsar --enable-swig-python
    if $build_dlls; then
	echo '#include "/usr/include/stdlib.h"
extern int setenv(const char *name, const char *value, int overwrite);
extern int unsetenv(const char *name);' > lalsimulation/src/stdlib.h
    fi
    echo -e "\\n\\n>> [`date`] Building lalsuite" >&3
    if $silent_build ; then
        make -j 2 2>&1 | tee >(grep Entering >&3 1>&3 2>&3) 
    else
        make
    fi
    echo -e "\\n\\n>> [`date`] Installing lalsuite" >&3
    make install
    for i in $PREFIX/etc/*-user-env.sh; do
        source "$i"
    done
    cd ..
    $cleanup && rm -rf lalsuite-build

    echo -e "\\n\\n>> [`date`] building $BUILDDIRNAME-preinst-lalsuite.tgz" >&3
    pushd "$PYCBC/.."
    tar -czf "$SOURCE/$BUILDDIRNAME-preinst-lalsuite.tgz" "$BUILDDIRNAME"
    popd

fi # if $BUILDDIRNAME-preinst.tgz

# Pegasus
if $build_pegasus ; then
    v=4.7.4
    p=pegasus-python-source-$v
    echo -e "\\n\\n>> [`date`] building $p" >&3
    test -r $p.tar.gz ||
      wget $wget_opts "$aei/$p.tar.gz" ||
      wget $wget_opts http://download.pegasus.isi.edu/pegasus/$v/$p.tar.gz
    pip install --no-deps $p.tar.gz
fi

# PyInstaller 9d0e0ad4 crashes with newer Jinja2,
# so install this old version before it gets pulled in from PyCBC
if $pyinstaller21_hacks; then
    echo -e "\\n\\n>> [`date`] pip install Jinja2==2.8.1" >&3
    pip install Jinja2==2.8.1
fi

# manually build subprocess32 with -D_GNU_SOURCE
# would be needed on old Linux with matplotlib>=2.0.0
if $build_subprocess32; then
    p=subprocess32-3.2.7
    echo -e "\\n\\n>> [`date`] building $p" >&3
    test -r $p.tar.gz ||
        wget $wget_opts "$aei/$p.tar.gz" ||
        wget $wget_opts $pypi/b8/2f/49e53b0d0e94611a2dc624a1ad24d41b6d94d0f1b0a078443407ea2214c2/$p.tar.gz
    rm -rf $p
    tar -xzf $p.tar.gz
    cd $p
    sed -i~ '/^# *define *HAVE_PIPE2 *1/d' _posixsubprocess.c
    python setup.py install --prefix="$PREFIX"
    cd ..
    $cleanup && rm -rf $p
fi

# on Windows, rebase DLLs
# doing this here _might_ fix a recurring problem with fork+git+PyCBC
# will be done again after building PyCBC
if $rebase_dlls_before_pycbc; then
    echo -e "\\n\\n>> [`date`] Rebasing DLLs" >&3
    find "$ENVIRONMENT" -name \*.dll > "$PREFIX/dlls.txt"
    rebase -d -b 0x61000000 -o 0x20000 -v -T "$PREFIX/dlls.txt"
fi

if $silent_build ; then
    # redirect stdout and stderr back to the screen
    exec 1>&- 
    exec 2>&- 
    exec 1>&3
    exec 2>&4
fi

# PyCBC
echo -e "\\n\\n>> [`date`] building pycbc"
if $scratch_pycbc || ! test -d pycbc/.git ; then
    # clone
    rm -rf pycbc
    git clone -n -o ligo-cbc git://github.com/ligo-cbc/pycbc
    cd pycbc
    git branch|grep ' master$'>/dev/null ||
        git checkout -b master ligo-cbc/master
    cd ..
fi
cd pycbc
if test ".$pycbc_remote" = ".ligo-cbc" ; then
    :
else
    git remote rm $pycbc_remote || true
    git remote add $pycbc_remote git://github.com/${pycbc_remote}/pycbc.git
    git remote update
fi
if test ".$pycbc_branch" = ".HEAD" ; then
    :
elif test ".$pycbc_branch" = ".master" ; then
    git checkout master
    git pull
    test ".$pycbc_commit" != "." &&
        git checkout $pycbc_commit
else
    # checkout branch from scratch, master must and should exist
    git reset --hard HEAD
    git checkout master
    git branch -D $pycbc_branch || true
    git remote update
    git checkout -b $pycbc_branch $pycbc_remote/$pycbc_branch
fi
if test "x$pycbc_fetch_ref" != "x" ; then
    git fetch $pycbc_remote +${pycbc_fetch_ref}
    git checkout FETCH_HEAD
fi

echo -e "[`date`] install pkgconfig beforehand"
pip install `grep -w ^pkgconfig requirements.txt||echo pkgconfig==1.1.0`
if $pyinstaller21_hacks; then
    echo -e "[`date`] install six and matplotlib beforehand"
    pip install `grep -w ^six requirements.txt||echo 'six>=1.9.0'`
    pip install `grep ^matplotlib== requirements.txt||echo matplotlib==1.4.3`
    p="`grep -w ^setuptools requirements.txt||true`"
    if test ".$p" != "."; then
        echo -e "[`date`] downgrade setuptools"
        pip install --upgrade $p
    fi
fi
echo -e "[`date`] git HEAD: `git log -1 --pretty=oneline --abbrev-commit`"
pycbc_tag="`git describe --tags --exact-match HEAD 2>/dev/null||true`"
pip install .
hooks="$PWD/tools/static"
cd ..
test -r "$PREFIX/etc/pycbc-user-env.sh" && source "$PREFIX/etc/pycbc-user-env.sh"

# on Windows, rebase DLLs
# from https://cygwin.com/ml/cygwin/2009-12/msg00168.html:
# /bin/rebase -d -b 0x61000000 -o 0x20000 -v -T <file with list of dll and so files> > rebase.out
if $build_dlls; then
    echo -e "\\n\\n>> [`date`] Rebasing DLLs" >&3
    find "$ENVIRONMENT" -name \*.dll > "$PREFIX/dlls.txt"
    rebase -d -b 0x61000000 -o 0x20000 -v -T "$PREFIX/dlls.txt"
fi

# TEST
echo -e "\\n\\n>> [`date`] testing local executable" >&3
$PREFIX/bin/pycbc_inspiral --help

if $silent_build ; then
    # close stdin and stdout
    exec 1>&-
    exec 2>&-

    # open stdout as $LOG_FILE file for read and write.
    exec 1<>$LOG_FILE

    # redirect stderr to stdout
    exec 2>&1
fi

# clean dist directory
rm -rf "$ENVIRONMENT/dist"
mkdir -p "$ENVIRONMENT/dist"

# if the build machine has dbhash & shelve, scipy weave will use bsddb
# make sure these exist and get added to the bundle(s)
python -c "import dbhash, shelve"
hidden_imports="--hidden-import=dbhash --hidden-import=shelve"

# PyInstaller
if echo "$pyinstaller_version" | egrep '^[0-9]\.[0-9][0-9]*$' > /dev/null; then
    # regular release version, get source tarball from pypi
    p=PyInstaller-$pyinstaller_version
    echo -e "\\n\\n>> [`date`] building $p" >&3
    test -r $p.tar.gz || wget $wget_opts "$pypi/source/P/PyInstaller/$p.tar.gz"
    rm -rf $p
    tar -xzf $p.tar.gz
    cd $p
else
    p=pyinstaller
    echo -e "\\n\\n>> [`date`] building $p" >&3
    if test -d pyinstaller/.git; then
        cd $p
        git remote update
        git reset --hard HEAD
        git clean -xdf
    else
        git clone git://github.com/$p/$p.git
        cd $p
    fi
    if test "$pyinstaller_version" = "v3.0"; then
        git checkout 3.0 # workaround for misnamed tag
    else
        git checkout $pyinstaller_version
        if test "$pyinstaller_version" = "develop"; then
            git pull
        fi
    fi
fi

# if we are to patch pyinstaller, save an unpatched version for later use
if $patch_pyinstaller_bootloader && $build_onefile_bundles; then
    python setup.py install --prefix="$PREFIX" --record installed-files.txt
    rm -f ../pyinstaller-clean-installed.tar ../pyinstaller-clean-installed.tar.gz
    xargs tar -rPf ../pyinstaller-clean-installed.tar < installed-files.txt
    gzip  ../pyinstaller-clean-installed.tar
    xargs rm -f < installed-files.txt
    rm installed-files.txt
fi

# patch PyInstaller to find the Python library on Cygwin
# This needs to be applied to compat.py for all PyI versions < 3.1.2
if $build_dlls; then
    if ! git log --pretty=oneline --abbrev-commit|grep "^e8c08c" >/dev/null; then
        sed -i~ "s|'libpython%d%d.dll'|'libpython%d.%d.dll'|" `find PyInstaller -name bindepend.py -or -name compat.py`
    fi
fi

# build bootloader (in any case: for Cygwin it wasn't precompiled, for Linux it was patched)
cd bootloader
# patch PyInstaller bootloader to not fork a second process
if $patch_pyinstaller_bootloader; then
    sed -i~ 's/ pid *= *fork()/ pid = 0/' */pyi_utils.c
fi
if $pyinstaller21_hacks; then
    python waf configure $pyinstaller_lsb build install
else
    test "$appendix" = "_OSX64" &&
        sed -i~ /-Wdeclaration-after-statement/d wscript
    pip install future
    python waf distclean configure $pyinstaller_lsb all
fi
cd ..
python setup.py install --prefix="$PREFIX" --record installed-files.txt
cd ..
test "$p" = "PyInstaller-$pyinstaller_version" && cleanup && rm -rf "$p"

if $build_wrapper; then
# on Linux, build "wrapper"
    echo -e "\\n\\n>> [`date`] Building 'BOINC wrapper', 'fstab' and 'fstab_test'" >&3
    if test -d boinc/.git ; then
	cd boinc
	git pull
    else
	# clone
	rm -rf boinc
	git clone $gitlab/boinc.git
	cd boinc
	git checkout -b eah_wrapper_improvements origin/eah_wrapper_improvements
	./_autosetup
	./configure LDFLAGS=-static-libgcc --disable-client --disable-manager --disable-server --enable-apps --disable-shared
    fi
    echo -e ">> [`date`] git HEAD: `git log -1 --pretty=oneline --abbrev-commit`" >&3
    make
    cp samples/wrapper/wrapper "$ENVIRONMENT/dist/wrapper$appendix"
    cd ..
    gcc -o "$ENVIRONMENT/dist/fstab" $SOURCE/pycbc/tools/einsteinathome/fstab.c
    gcc -DTEST_WIN32 -o "$ENVIRONMENT/dist/fstab_test" $SOURCE/pycbc/tools/einsteinathome/fstab.c
elif $build_fstab ; then
    echo -e "\\n\\n>> [`date`] Building 'fstab.exe'" >&3
    if $build_dlls; then
        x86_64-w64-mingw32-gcc -o "$ENVIRONMENT/dist/fstab$appendix.exe" $SOURCE/pycbc/tools/einsteinathome/fstab.c
    else
        gcc -o "$ENVIRONMENT/dist/fstab$appendix" $SOURCE/pycbc/tools/einsteinathome/fstab.c
    fi
fi

# log environment
echo -e "\\n\\n>> [`date`] ENVIRONMENT ..." >&3
env
echo -e "... ENVIRONMENT"

cd $PREFIX

# BUNDLE DIR
echo -e "\\n\\n>> [`date`] building pyinstaller spec" >&3
# don't use UPX with PyInstaller
upx="--noupx"
if $build_dlls && ! $pyinstaller21_hacks; then
    ext=".exe"
fi
# create spec file
if $use_pycbc_pyinstaller_hooks; then
    export NOW_BUILDING=NULL
    export PYCBC_HOOKS_DIR="$hooks"
    pyi-makespec $upx --additional-hooks-dir $hooks/hooks --runtime-hook $hooks/runtime-tkinter.py $hidden_imports --hidden-import=pkg_resources --onedir --name pycbc_inspiral$ext ./bin/pycbc_inspiral
else
    # find hidden imports (pycbc CPU modules)
    hidden_imports=`find $PREFIX/lib/python2.7/site-packages/pycbc/ -name '*_cpu.py' | sed 's%.*/site-packages/%%;s%\.py$%%;s%/%.%g;s%^% --hidden-import=%' | tr -d '\012'`
    pyi-makespec $upx $hidden_imports --hidden-import=scipy.linalg.cython_blas --hidden-import=scipy.linalg.cython_lapack --hidden-import=pkg_resources --onedir --name pycbc_inspiral$ext ./bin/pycbc_inspiral
fi
# patch spec file to add "-v" to python interpreter options
if $verbose_pyinstalled_python; then
    sed -i~ 's%exe = EXE(pyz,%options = [ ("v", None, "OPTION"), ("W error", None, "OPTION") ]\
exe = EXE(pyz, options,%' pycbc_inspiral$ext.spec
fi
echo -e "\\n\\n>> [`date`] running pyinstaller" >&3
pyinstaller $upx pycbc_inspiral$ext.spec

test ".$ext" != "." &&
    mv dist/pycbc_inspiral$ext dist/pycbc_inspiral
cd dist/pycbc_inspiral

# fix some libraries manually
if test -r "libgcc_s.1.dylib"; then
    cp "`$CC -print-file-name=libgcc_s.1.dylib`" .
fi
if test -r /usr/bin/cyggomp-1.dll; then
    cp /usr/bin/cyggomp-1.dll .
else
    libgomp=`gcc -print-file-name=libgomp.so.1`
    if test "$libgomp" != "libgomp.so.1"; then
        cp "$libgomp" .
    fi
fi

# OSX doesn't have a GNU error C extension, so drop an "error.h" header
# with a fake 'error' function somewhere for scipy wave to pick it up
if test ".$appendix" = "._OSX64"; then
    if test -d scipy/weave; then
        f=scipy/weave/error.h
    else
        f=weave/error.h
    fi
    echo '#define error(status, errnum, errstr, ...) fprintf(stderr,"pycbc_inspiral: %d:%d:" errstr, status, errnum, ##__VA_ARGS__)' > $f
fi

# TEST BUNDLE
echo -e "\\n\\n>> [`date`] Testing pycbc_inspiral bundle" >&3
./pycbc_inspiral --help
cd ..

# build zip file from dir
echo -e "\\n\\n>> [`date`] Creating zip archive from directory" >&3
zip -r pycbc_inspiral$appendix.zip pycbc_inspiral

# if the executable is "pycbc_inspiral.exe", add a "XML soft link" "pycbc_inspiral" to the bundle for the wrapper
if $build_dlls; then
    mkdir -p tmp/pycbc_inspiral
    echo '<soft_link>pycbc_inspiral/pycbc_inspiral.exe<soft_link/>' > tmp/pycbc_inspiral/pycbc_inspiral
    cd tmp
    zip ../pycbc_inspiral$appendix.zip pycbc_inspiral/pycbc_inspiral
fi

# run 10min self-test, build wave cache
cd "$SOURCE"
mkdir -p test
cd test

if $run_analysis; then
echo -e "\\n\\n>> [`date`] running analysis" >&3
p="H-H1_LOSC_4_V1-1126257414-4096.gwf"
md5="a7d5cbd6ef395e8a79ef29228076d38d"
if check_md5 "$p" "$md5"; then
    rm -f "$p"
    wget $wget_opts "$albert/$p"
    if check_md5 "$p" "$md5"; then
        echo "can't download $p - md5 mismatch"
        exit 1
    fi
fi
frames="$PWD/$p"

# free some space by removing a possible old ROM tarball
rm -f lal-data-r7.tar.gz

failed=false
while read f md5; do
    if check_md5 "$f" "$md5"; then
	failed=true
	break
    else
	roms="$f $roms"
    fi
done <<EOF
SEOBNRv2ChirpTimeSS.dat                  7b7dbadacc3f565fb2c8e6971df2ab74
SEOBNRv2ROM_DS_sub1_AmpPrefac_ci.dat     62afa5351d6b775ac33cb4d898f0016b
SEOBNRv2ROM_DS_sub1_Amp_ciall.dat        f82ddc5dc0b6fdc75122e767bd5e78c8
SEOBNRv2ROM_DS_sub1_Bamp_bin.dat         a6829fa05437cc0aad81e3f8dae839cc
SEOBNRv2ROM_DS_sub1_Bphase_bin.dat       98ea14b01e729d15ff666caa25afaed6
SEOBNRv2ROM_DS_sub1_Phase_ciall.dat      b41f0f7fbaf8be1d1848de7ee702bc67
SEOBNRv2ROM_DS_sub2_AmpPrefac_ci.dat     96c384617edd8375ceaa03f9b7456467
SEOBNRv2ROM_DS_sub2_Amp_ciall.dat        20ee260c870109766a6f048e20c7e10f
SEOBNRv2ROM_DS_sub2_Bamp_bin.dat         67d4f206fe19104fbc98b923b37318bb
SEOBNRv2ROM_DS_sub2_Bphase_bin.dat       d0bf97b4e17b5c9a7cfd222aaaafd742
SEOBNRv2ROM_DS_sub2_Phase_ciall.dat      c2ea5d296fee01abe16c0dd9e5f71f04
SEOBNRv2ROM_DS_sub3_AmpPrefac_ci.dat     4d5378935a7fba5e96f671581bce99fb
SEOBNRv2ROM_DS_sub3_Amp_ciall.dat        412953726ca4bc72a810b27b810831c7
SEOBNRv2ROM_DS_sub3_Bamp_bin.dat         31f48cb651a60837a3e99ee050aa9bc2
SEOBNRv2ROM_DS_sub3_Bphase_bin.dat       727d31f6dc678aba8539817c8d0ae930
SEOBNRv2ROM_DS_sub3_Phase_ciall.dat      d0e1601c7cf4bd727d03e6cf7d2f722b
SEOBNRv2ROM_SS_AmpPrefac_ci.dat          08186a21682d2e73cb00a3ef35aa5c9c
SEOBNRv2ROM_SS_Amp_ciall.dat             e6c243f76609cada55612cfe53f82e41
SEOBNRv2ROM_SS_Bamp_bin.dat              1ef7953a977a1fb551f59585c5d63d7a
SEOBNRv2ROM_SS_Bphase_bin.dat            b5923860bf021e6a2a23d743e5724bee
SEOBNRv2ROM_SS_Phase_ciall.dat           2947032d0ad7ffde9704e24bf9e676f5
SEOBNRv4ROM_v2.0.hdf5                    256f33aecada770b871750227fa87749
EOF
if $failed; then
    p="lal-data-r11.tar.gz"
    md5="d4009ac3328fa55b654c9fafc1006a5c"
    rm -f "$p"
    wget $wget_opts "$albert/$p"
    if check_md5 "$p" "$md5"; then
        echo "can't download $p - md5 mismatch"
        exit 1
    fi
    roms=`tar -zxvf $p`
fi

#fb5ec108c69f9e424813de104731370c  H1L1-PREGEN_TMPLTBANK_SPLITBANK_BANK16-1126051217-3331800-short2k.xml.gz
p="H1L1-SBANK_FOR_GW150914ER10.xml.gz"
md5="c24f5513d3066b4f637daffb6aa20fec"
if check_md5 "$p" "$md5"; then
    rm -f "$p"
    wget $wget_opts "$albert/$p"
    if check_md5 "$p" "$md5"; then
        echo "can't download $p - md5 mismatch"
        exit 1
    fi
fi

# Fetch any extra libraries specified on the command line
if [ ! -z ${extra_libs} ] ; then
  echo -e "\\n\\n>> [`date`] installing extra libraries from ${extra_libs} to $PREFIX/lib" >&3
  curl --ftp-pasv --insecure -o extra_libs.tar.gz ${extra_libs}
  tar -C $PREFIX/lib -zxvf extra_libs.tar.gz
fi

if $build_framecpp; then
    export LAL_FRAME_LIBRARY=FrameC
else
    export LAL_FRAME_LIBRARY=FrameL
fi

gw150914_bank=$p
gw150914_approx='SPAtmplt:mtotal<4 SEOBNRv4_ROM:mtotal<20 SEOBNRv2_ROM_DoubleSpin:else'
bank_array=( "$gw150914_bank" )
approx_array=( "$gw150914_approx" )

if ! test -z "$extra_approx" || ! test -z "$extra_bank" ; then
    if test -z "$extra_bank" ; then
        bank_array=( "$gw150914_bank" "$gw150914_bank" )
    else
        bank_array=( "$extra_bank" "$gw150914_bank" )
    fi
    if test -z "$extra_approx" ; then
        approx_array=( "$gw150914_approx" "$gw150914_approx" )
    else
        approx_array=( "$extra_approx" "$gw150914_approx" )
    fi
fi

n_runs=${#bank_array[@]}

if $silent_build ; then
    # redirect stdout and stderr back to the screen
    exec 1>&-
    exec 2>&-
    exec 1>&3
    exec 2>&4
fi

for (( i=0; i<${n_runs}; i++ ))
do
    rm -f H1-INSPIRAL-OUT.hdf
    echo "\
>> [`date`] pycbc_inspiral using
>>   --bank-file ${bank_array[$i]}
>>   --approximant ${approx_array[$i]}
>>   ROM data from $lal_data_path"
    export \
      CPPFLAGS="$CPPFLAGS `python-config --includes`" \
      LAL_DATA_PATH="$lal_data_path" \
      NO_TMPDIR=1 \
      INITIAL_LOG_LEVEL=10 \
      LEVEL2_CACHE_SIZE=8192 \
      WEAVE_FLAGS='-O3 -march=core2 -w' \
      FIXED_WEAVE_CACHE="$PWD/pycbc_inspiral"
    args="--fixed-weave-cache ${processing_scheme} \
      --sample-rate 2048 \
      --sgchisq-snr-threshold 6.0 \
      --sgchisq-locations mtotal>40:20-30,20-45,20-60,20-75,20-90,20-105,20-120 \
      --segment-end-pad 16 \
      --cluster-method window \
      --low-frequency-cutoff 30 \
      --pad-data 8 \
      --cluster-window 1 \
      --cluster-function symmetric \
      --injection-window 4.5 \
      --segment-start-pad 112 \
      --psd-segment-stride 8 \
      --psd-inverse-length 16 \
      --filter-inj-only \
      --psd-segment-length 16 \
      --snr-threshold 5.5 \
      --segment-length 256 \
      --autogating-width 0.25 \
      --autogating-threshold 100 \
      --autogating-cluster 0.5 \
      --autogating-taper 0.25 \
      --newsnr-threshold 5 \
      --psd-estimation median \
      --strain-high-pass 20 \
      --order -1 \
      --chisq-bins 1.75*(get_freq('fSEOBNRv2Peak',params.mass1,params.mass2,params.spin1z,params.spin2z)-60.)**0.5 \
      --channel-name H1:LOSC-STRAIN \
      --gps-start-time 1126259078 \
      --gps-end-time 1126259846 \
      --output H1-INSPIRAL-OUT.hdf \
      --frame-files $frames \
      --approximant ${approx_array[$i]} \
      --bank-file ${bank_array[$i]} \
      --verbose"
    if $silent_build; then
        "$ENVIRONMENT/dist/pycbc_inspiral/pycbc_inspiral" $args 2>&1|
          awk '{if ((!/Filtering template|points above|power chisq|point chisq|Found chisq|generating SEOBNR|generating SPA/) || (/segment 1/ && NR % 50 == 0) || / 0: generating/ || / 1: generating/ ) print}'
    else
        "$ENVIRONMENT/dist/pycbc_inspiral/pycbc_inspiral" $args
    fi
done

# test for GW150914
echo -e "\\n\\n>> [`date`] test for GW150914"
python $SOURCE/pycbc/tools/einsteinathome/check_GW150914_detection.py H1-INSPIRAL-OUT.hdf

fi # if $run_analysis

# zip weave cache
echo -e "\\n\\n>> [`date`] zipping weave cache"
cache="$ENVIRONMENT/dist/pythoncompiled$appendix.zip"
rm -f "$cache"
# Fetch any extra libraries specified on the command line
if [ ! -z ${extra_libs} ] ; then
  echo -e "\\n\\n>> [`date`] adding extra libraries from ${extra_libs} to bundle" >&3
  tar -C pycbc_inspiral/ -zxvf extra_libs.tar.gz
fi
# addin all ROM files to the cache would blow it up to >300MB, so for now add only the one
# that is actually used in the GW150914 analysis. Use '$roms' instead of
# 'SEOBNRv2ChirpTimeSS.dat' if you want all to be included, currently +280MB
zip -r "$cache" pycbc_inspiral SEOBNRv2ChirpTimeSS.dat

# build additional PyInstaller "onefile" bundles (with unpatched PyInstaller bootloader)
if $build_onefile_bundles; then
    # restore unpatched pyinstaller version
    if $patch_pyinstaller_bootloader; then
        echo -e "\\n\\n>> [`date`] restore unpatched pyinstaller version" >&3
        cd $SOURCE/pyinstaller
        xargs rm -f < installed-files.txt
        tar -xPzf ../pyinstaller-clean-installed.tar.gz
    fi

    export NOW_BUILDING=NULL
    cd "$ENVIRONMENT"

    echo -e "\\n\\n>> [`date`] build pycbc_condition_strain bundle" >&3
    pyinstaller \
        --additional-hooks-dir $hooks/hooks \
        --runtime-hook $hooks/runtime-tkinter.py \
        --runtime-hook $hooks/runtime-scipy.py \
        --hidden-import=pkg_resources $hidden_imports \
        --onefile ./bin/pycbc_condition_strain

    echo -e "\\n\\n>> [`date`] building pycbc_inspiral_osg pyinstaller onefile spec" >&3
    pyi-makespec \
        --additional-hooks-dir $hooks/hooks \
        --runtime-hook $hooks/runtime-tkinter.py \
        --hidden-import=pkg_resources $hidden_imports \
        --onefile ./bin/pycbc_inspiral --name pycbc_inspiral_osg

    echo -e ">> [`date`] patching pyinstaller spec" >&3
    mv pycbc_inspiral_osg.spec pycbc_inspiral_osg.spec1
    ls $SOURCE/test/pycbc_inspiral/ |
        sed "s%.*%a.datas += [('&','$SOURCE/test/pycbc_inspiral/&','DATA')]%" > pycbc_inspiral_osg.spec2
    if test ".$libgomp" != "." -a "$libgomp" != "libgomp.so.1"; then
        echo 'a.datas += [("libgomp.so.1", "'"$libgomp"'", "DATA")]' >> pycbc_inspiral_osg.spec2
    fi
    awk '/exe = EXE/ { while (getline l < "pycbc_inspiral_osg.spec2") {print l} }; {print}' \
        pycbc_inspiral_osg.spec1 > pycbc_inspiral_osg.spec

    echo -e ">> [`date`] running pyinstaller" >&3
    pyinstaller pycbc_inspiral_osg.spec

    if echo ".$pycbc_tag" | egrep '^\.v[0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*$' >/dev/null; then
        mv dist/pycbc_inspiral_osg "dist/pycbc_inspiral_osg_$pycbc_tag"
    fi
fi

# remove lock
rm -f "$PYCBC/lock"

# remove log file
if [ -f $LOG_FILE ] ; then
    rm -f $LOG_FILE
fi

echo -e "\\n\\n>> [`date`] Success $0"
