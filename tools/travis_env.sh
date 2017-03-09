BUILD=${HOME}/build
BUILDDIRNAME="pycbc-build"
PYCBC="$BUILD/$BUILDDIRNAME"
PYTHON_PREFIX="$PYCBC"
ENVIRONMENT="$PYCBC/environment"
PREFIX="$ENVIRONMENT"
PATH="$PREFIX/bin:$PYTHON_PREFIX/bin:$PATH"
export
LD_LIBRARY_PATH="$PREFIX/lib:$PREFIX/bin:$PYTHON_PREFIX/lib:/usr/local/lib:$LD_LIBRARY_PATH"
export
PKG_CONFIG_PATH="$PREFIX/lib/pkgconfig:$PYTHON_PREFIX/lib/pkgconfig:/usr/local/lib/pkgconfig:$PKG_CONFIG_PATH"
source ${BUILD}/pycbc-build/environment/etc/lalsuite-user-env.sh
source ${BUILD}/pycbc-build/environment/bin/activate
export LAL_DATA_PATH=${BUILD}/pycbc-sources/test

