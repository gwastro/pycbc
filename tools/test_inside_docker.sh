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

cat /etc/redhat-release

if [ "x${OS_VERSION}" == "x6" ] ; then
  echo -e "\\n>> [`date`] Building pycbc_inspiral bundle for CentOS 6"

  # install requirements into docker container
  yum install -y gcc gcc-c++ gcc-gfortran python-devel pcre-devel autoconf automake make zlib-devel libpng-devel libjpeg-devel libsqlite3-dev sqlite-devel
  ln -s /usr/bin/g++ /usr/bin/g++-4.4.7
  ln -s /usr/bin/gfortran /usr/bin/gfortran-4.4.7
  yum -y install epel-release
  yum -y install ant asciidoc fop docbook-style-xsl.noarch R-devel

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

  rpm -ivh http://software.ligo.org/lscsoft/scientific/7.2/x86_64/production/lscsoft-production-config-1.3-1.el7.noarch.rpm
  yum clean all
  yum makecache 
  yum update
  yum -y install lscsoft-backports-config
  yum -y install lscsoft-epel-config
  curl http://download.pegasus.isi.edu/wms/download/rhel/7/pegasus.repo > /etc/yum.repos.d/pegasus.repo
  yum clean all
  yum makecache
  yum update
  yum -y install lscsoft-ius-config
  yum clean all
  yum makecache
  yum -y install git2u-all
  yum -y install lscsoft-all

  TRAVIS_TAG="vX.Y.Z"
  CVMFS_PATH=/cvmfs/oasis.opensciencegrid.org/ligo/sw/pycbc/x86_64_rhel_7/virtualenv
  mkdir -p ${CVMFS_PATH}
  VENV_PATH=${CVMFS_PATH}/pycbc-${TRAVIS_TAG}
  pip install virtualenv
  virtualenv ${VENV_PATH}
  echo 'export PYTHONUSERBASE=${VIRTUAL_ENV}/.local' >> ${VENV_PATH}/bin/activate
  echo 'export XDG_CACHE_HOME=${VIRTUAL_ENV}/.cache' >> ${VENV_PATH}/bin/activate
  source ${VENV_PATH}/bin/activate
  pip install --upgrade pip
  pip install six packaging appdirs
  pip install --upgrade setuptools

  deactivate
fi 

echo -e "\\n>> [`date`] CentOS Docker script exiting"

exit 0
