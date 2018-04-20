#!/bin/bash

set -e

# update apt
grep '^deb http://archive.debian.org/debian/' /etc/apt/sources.list ||
  echo 'deb http://archive.debian.org/debian/ etch main' >>/etc/apt/sources.list
apt-get -y update || true
apt-get -y install debian-keyring debian-archive-keyring
apt-get -y update
# tools for pycbc
apt-get -y install git-core gcc g++ gfortran automake autoconf make libtool pkg-config bzip2
# libraries for pycbc
apt-get -y install libpcre3-dev libfreetype6-dev libjpeg-dev libpng-dev libmysqlclient-dev libpq-dev libssl-dev libsqlite3-dev libdb4.4-dev
# make sure these libraries aren't on the system
apt-get -y remove lapack3-dev atlas3-base atlas3-headers refblas3
# further necessary / useful tools
apt-get -y install openssh-server ntpdate zip gettext curl libcurl3-openssl-dev wget bzip2 libbz2-dev screen emacs

test ".$PREFIX" = "." &&
  PREFIX=/usr/local

export PATH="$PREFIX:$PATH" LD_LIBRARY_PATH="$PREFIX/lib:$LD_LIBRARY_PATH"

test ".$GCC" = "." &&
  GCC=4.8.5

if fgrep -q -x "$PREFIX/lib" /etc/ld.so.conf; then
  echo "$PREFIX/lib already in /etc/ld.so.conf"
else
  echo "$PREFIX/lib" >> /etc/ld.so.conf
fi

compile_tar_gz() {
    echo "==== `date` compiling $1"
    pkg="$1"
    shift
    tar -xzvf "$pkg".tar.gz 
    cd "$pkg"
    ./configure "$@" "--prefix=$PREFIX"
    make
    make install
    cd ..
    rm -r "$pkg"
    ldconfig
    hash -r
}

compile_gcc() {
    echo "==== `date` compiling gcc-$1"
    vers="$1"
    shift
    tar -xjvf gcc-${vers}.tar.bz2 
    mkdir -p gcc-${vers}/build
    cd gcc-${vers}/build
    ../configure "$@" \
        "--prefix=$PREFIX/gcc-$vers" \
        "--with-local-prefix=$PREFIX/gcc-$vers" \
	"--program-suffix=-$vers" \
        --enable-languages=c,c++,objc,fortran \
        --disable-symvers --enable-shared
    make
    make install
    cd ../..
    rm -r gcc-${vers}

    ( cd "$PREFIX/bin" &&
	rm -f gcc-${vers} gfortran-${vers} g++-${vers} &&
	ln -s ../gcc-${vers}/bin/gcc-${vers} &&
	ln -s ../gcc-${vers}/bin/gfortran-${vers} &&
	ln -s ../gcc-${vers}/bin/g++-${vers} )

    ldconfig
    hash -r
}

gnu=ftp://ftp.fu-berlin.de/unix/gnu

for url in `echo "
  $gnu/m4/m4-1.4.9.tar.gz
  $gnu/automake/automake-1.11.tar.gz
  $gnu/autoconf/autoconf-2.63.tar.gz
  $gnu/libtool/libtool-1.5.8.tar.gz
  $gnu/gmp/gmp-4.3.2.tar.gz
  $gnu/mpfr/mpfr-2.4.2.tar.gz
  $gnu/binutils/binutils-2.24.tar.gz
  $gnu/gcc/gcc-$GCC/gcc-$GCC.tar.bz2
  http://www.aei.uni-hannover.de/~bema/zlib-1.2.6.tar.gz
  http://www.aei.uni-hannover.de/~bema/jre-6u45-linux-x64.bin
  http://www.multiprecision.org/mpc/download/mpc-0.9.tar.gz
  http://pkgconfig.freedesktop.org/releases/pkg-config-0.23.tar.gz
  http://www.python.org/ftp/python/2.7.10/Python-2.7.10.tgz
  https://www.kernel.org/pub/software/scm/git/git-1.9.5.tar.gz
  https://curl.haxx.se/download/curl-7.47.1.tar.gz
"` ; do
  test -r `echo $url | sed 's%.*/%%'` ||
  wget --no-check-certificate --passive-ftp "$url" || exit
done

compile_tar_gz m4-1.4.9
compile_tar_gz autoconf-2.63
compile_tar_gz automake-1.11
compile_tar_gz libtool-1.5.8
compile_tar_gz pkg-config-0.23
compile_tar_gz curl-7.47.1 "--with-ssl"

tar -xzvf git-1.9.5.tar.gz
cd git-1.9.5
make configure
./configure --prefix=$PREFIX --with-tcltk
make && make install
cd ..
rm -rf git-1.9.5

compile_tar_gz gmp-4.3.2 "--build=`gcc -dumpmachine`"
compile_tar_gz mpfr-2.4.2 "--with-gmp=$PREFIX"
compile_tar_gz mpc-0.9 "--with-gmp=$PREFIX" "--with-mpfr=$PREFIX"
compile_tar_gz binutils-2.24
compile_gcc $GCC "--with-gmp=$PREFIX" "--with-mpfr=$PREFIX" "--with-mpc=$PREFIX" --disable-multilib
# "--target=`gcc -dumpmachine`" --disable-multilib --with-multilib-list=m64
( cd $PREFIX/bin && for i in gcc g++ gfortran; do rm -f $i && ln -s $i-$GCC $i; done )

rm -f Python-2.7.10.tar.gz
ln -s Python-2.7.10.tgz Python-2.7.10.tar.gz
compile_tar_gz Python-2.7.10 --enable-shared
python -m ensurepip
pip install --upgrade pip
pip install --trusted-host pypi.python.org virtualenv

p=$PWD/jre-6u45-linux-x64.bin
chmod +x $p
mkdir -p /opt
cd /opt
$p
ln -s jre1.6.0_45 java

echo "==== `date` done."
