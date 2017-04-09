#!/bin/bash

if [ "x${OS_NAME}" != "xubuntu" ] ; then
  set -e
  echo -e "\\n>> [`date`] Starting ${OS_NAME} ${OS_VERSION} docker container"
  if [ "x${OS_NAME}" == "xcentos" ] ; then DOCKER_IMG="pycbc/ldg-el7" ; fi
  if [ "x${OS_NAME}" == "xscientific" ] ; then DOCKER_IMG="pycbc/sl6-travis" ; fi
  if [ "x${TRAVIS_SECURE_ENV_VARS}" == "xtrue" ] ; then
    cp -R ~/.ssh .
  fi
  mkdir -p build/pycbc-sources
  if [ -f "$HOME/docker-cache/pycbc-build-preinst.tgz" ] ; then
    cp $HOME/docker-cache/pycbc-build-preinst.tgz build/pycbc-sources
  fi
  if [ -f "$HOME/docker-cache/pycbc-build-preinst-lalsuite.tgz" ] ; then
    cp $HOME/docker-cache/pycbc-build-preinst-lalsuite.tgz build/pycbc-sources
  fi
  if [ -f "$HOME/docker-cache/test" ] ; then
    cp -a $HOME/docker-cache/test build/pycbc-sources/test
  fi
  sudo docker run --name buildvm -v `pwd`:/pycbc:rw ${DOCKER_IMG} /bin/bash -c "bash /pycbc/tools/docker_build_dist.sh --os-version=${OS_VERSION} --pull-request=${TRAVIS_PULL_REQUEST} --commit=${TRAVIS_COMMIT} --secure=${TRAVIS_SECURE_ENV_VARS} --tag=${TRAVIS_TAG}"
  echo -e "\\n>> [`date`] CentOS Docker exited"
  if [ "x${OS_NAME}" == "xscientific" ] ; then
    echo -e "\\n>> [`date`] Caching E@H build environment"
    mkdir -p $HOME/docker-cache
    mkdir -p $HOME/docker-cache/test
    sudo docker cp buildvm:/pycbc/build/pycbc-sources/pycbc-build-preinst.tgz $HOME/docker-cache/pycbc-build-preinst.tgz
    sudo docker cp buildvm:/pycbc/build/pycbc-sources/pycbc-build-preinst-lalsuite.tgz $HOME/docker-cache/pycbc-build-preinst-lalsuite.tgz
    sudo docker cp buildvm:/pycbc/build/pycbc-sources/test/H1L1-SBANK_FOR_GW150914ER10.xml.gz $HOME/docker-cache/test/H1L1-SBANK_FOR_GW150914ER10.xml.gz
    sudo docker cp buildvm:/pycbc/build/pycbc-sources/test/H-H1_LOSC_4_V1-1126257414-4096.gwf $HOME/docker-cache/test/H-H1_LOSC_4_V1-1126257414-4096.gwf
  fi
  exit 0
fi

echo -e "\\n>> [`date`] Starting PyCBC test suite"

LOG_FILE=$(mktemp -t pycbc-test-log.XXXXXXXXXX)

BUILD=${HOME}/build
BUILDDIRNAME="pycbc-build"
PYCBC="$BUILD/$BUILDDIRNAME"
PYTHON_PREFIX="$PYCBC"
ENVIRONMENT="$PYCBC/environment"
PREFIX="$ENVIRONMENT"
PATH="$PREFIX/bin:$PYTHON_PREFIX/bin:$PATH"
export LD_LIBRARY_PATH="$PREFIX/lib:$PREFIX/bin:$PYTHON_PREFIX/lib:/usr/local/lib:$LD_LIBRARY_PATH"
export PKG_CONFIG_PATH="$PREFIX/lib/pkgconfig:$PYTHON_PREFIX/lib/pkgconfig:/usr/local/lib/pkgconfig:$PKG_CONFIG_PATH"
source ${BUILD}/pycbc-build/environment/etc/lalsuite-user-env.sh
source ${BUILD}/pycbc-build/environment/bin/activate
export LAL_DATA_PATH=$HOME/build/pycbc-sources/test

RESULT=0

# Using python setup.py test has two issues:
#     Some tests fail for reasons not necessarily related to PyCBC
#     Setup.py seems to returns 0 even when tests fail
# So we rather run specific tests manually
for prog in `find test -name '*.py' -print | egrep -v '(autochisq|bankveto|fft|schemes|long|lalsim|test_waveform)'`
do 
    echo -e ">> [`date`] running unit test for $prog"
    python $prog &> $LOG_FILE
    if test $? -ne 0 ; then
        RESULT=1
        echo -e "    FAILED!"
        echo -e "---------------------------------------------------------"
        cat $LOG_FILE
        echo -e "---------------------------------------------------------"
    else
        echo -e "    Pass."
    fi
done

# check that all executables that do not require
# special environments can return a help message
for prog in `find ${PATH//:/ } -maxdepth 1 -name 'pycbc*' -print 2>/dev/null | egrep -v '(pycbc_live|pycbc_live_nagios_monitor|pycbc_make_grb_summary_page|pycbc_make_offline_grb_workflow|pycbc_mvsc_get_features|pycbc_upload_xml_to_gracedb)'`
do
    echo -e ">> [`date`] running $prog --help"
    $prog --help &> $LOG_FILE
    if test $? -ne 0 ; then
        RESULT=1
        echo -e "    FAILED!"
        echo -e "---------------------------------------------------------"
        cat $LOG_FILE
        echo -e "---------------------------------------------------------"
    else
        echo -e "    Pass."
    fi
done

echo -e "\\n>> [`date`] Building documentation"

python setup.py build_gh_pages &> $LOG_FILE
if test $? -ne 0 ; then
    echo -e "    FAILED!"
    echo -e "---------------------------------------------------------"
    cat $LOG_FILE
    echo -e "---------------------------------------------------------"
    RESULT=1
fi

ls -alrt

ls _gh_pages/latest

ls _gh_pages/latest/_downloads

ls docs

exit ${RESULT}
