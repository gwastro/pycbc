#!/bin/bash -v

set -e

OS_VERSION=${1}
TRAVIS_TAG=${2}
PYCBC_CODE=${3}
LALSUITE_CODE=${4}

echo -e "\\n>> [`date`] Inside CentOS ${OS_VERSION}"
echo -e "\\n>> [`date`] Release tag is ${TRAVIS_TAG}"
echo -e "\\n>> [`date`] Using PyCBC code ${PYCBC_CODE}"
echo -e "\\n>> [`date`] Using lalsuite code ${LALSUITE_CODE}"

cat /etc/issue

if [ "x${OS_VERSION}" == "x6" ] ; then
  echo -e "\\n>> [`date`] Building pycbc_inspiral bundle for CentOS 6"

  # create working dir for build script
  BUILD=/pycbc/build
  mkdir -p ${BUILD}
  export PYTHONUSERBASE=${BUILD}/.local
  export XDG_CACHE_HOME=${BUILD}/.cache

  # run the einstein at home build and test script
  pushd ${BUILD}
  /pycbc/tools/einsteinathome/pycbc_build_eah.sh ${LALSUITE_CODE} ${PYCBC_CODE} --silent-build
  popd
fi

if [ "x${OS_VERSION}" == "x7" ] ; then
  echo -e "\\n>> [`date`] Building pycbc virtual environment for CentOS 7"
fi 

echo -e "\\n>> [`date`] CentOS Docker script exiting"

exit 0
