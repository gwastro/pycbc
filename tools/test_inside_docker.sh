#!/bin/bash -v

set -e

OS_VERSION=${1}
TRAVIS_TAG=${2}
TRAVIS_PULL_REQUEST=${3}
TRAVIS_COMMIT=${4}
TRAVIS_SECURE_ENV_VARS=${5}

# determine the pycbc git branch and origin
yum -q -y install git
pushd /pycbc
git branch -vvv
if test x$TRAVIS_PULL_REQUEST = "xfalse" ; then
    PYCBC_CODE="--pycbc-commit=${TRAVIS_COMMIT}"
else
    PYCBC_CODE="--pycbc-fetch-ref=refs/pull/${TRAVIS_PULL_REQUEST}/merge"
fi
popd

# set the lalsuite checkout to use
LALSUITE_HASH="539c8700af92eb6dd00e0e91b9dbaf5bae51f004"

if [ "x${TRAVIS_SECURE_ENV_VARS}" == "xtrue" ] ; then
  mkdir -p ~/.ssh
  cp /pycbc/.ssh/* ~/.ssh
  chmod 600 ~/.ssh/id_rsa
fi

echo -e "\\n>> [`date`] Inside CentOS ${OS_VERSION}"
echo -e "\\n>> [`date`] Release tag is ${TRAVIS_TAG}"
echo -e "\\n>> [`date`] Using PyCBC code ${PYCBC_CODE}"
echo -e "\\n>> [`date`] Using lalsuite hash ${LALSUITE_HASH}"

if [ "x${OS_VERSION}" == "x6" ] ; then
  echo -e "\\n>> [`date`] Building pycbc_inspiral bundle for CentOS 6"

  # install requirements into docker container
  yum install -q -y gcc gcc-c++ gcc-gfortran python-devel pcre-devel autoconf automake tar zlib-devel libpng-devel libjpeg-devel libsqlite3-dev sqlite-devel wget db4-devel
  ln -s /usr/bin/g++ /usr/bin/g++-4.4.7
  ln -s /usr/bin/gcc /usr/bin/gcc-4.4.7
  ln -s /usr/bin/gfortran /usr/bin/gfortran-4.4.7

  # create working dir for build script
  BUILD=/pycbc/build
  mkdir -p ${BUILD}
  export PYTHONUSERBASE=${BUILD}/.local
  export XDG_CACHE_HOME=${BUILD}/.cache

  # run the einstein at home build and test script
  pushd ${BUILD}
  /pycbc/tools/einsteinathome/pycbc_build_eah.sh --lalsuite-commit=${LALSUITE_HASH} ${PYCBC_CODE} --silent-build
  popd
fi

if [ "x${OS_VERSION}" == "x7" ] ; then
  echo -e "\\n>> [`date`] Building pycbc virtual environment for CentOS 7"

  rpm -ivh http://software.ligo.org/lscsoft/scientific/7.2/x86_64/production/lscsoft-production-config-1.3-1.el7.noarch.rpm
  yum clean all
  yum makecache 
  yum update
  yum -q -y install lscsoft-backports-config
  yum -q -y install lscsoft-epel-config
  curl http://download.pegasus.isi.edu/wms/download/rhel/7/pegasus.repo > /etc/yum.repos.d/pegasus.repo
  yum clean all
  yum makecache
  yum update
  yum -q -y install lscsoft-ius-config
  yum clean all
  yum makecache
  yum --debuglevel=1 -y install lscsoft-all
  yum install -q -y zlib-devel libpng-devel libjpeg-devel libsqlite3-dev sqlite-devel db4-devel

  rpm --nodeps -e `rpm -qa | grep lal`

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

  pip install "numpy>=1.6.4" "h5py>=2.5" unittest2 python-cjson Cython decorator
  SWIG_FEATURES="-cpperraswarn -includeall -I/usr/include/openssl" pip install M2Crypto

  mkdir -p ${VIRTUAL_ENV}/src
  cd ${VIRTUAL_ENV}/src
  git clone https://github.com/lscsoft/lalsuite.git
  cd ${VIRTUAL_ENV}/src/lalsuite
  git checkout ${LALSUITE_HASH}
  ./00boot
  ./configure --prefix=${VIRTUAL_ENV}/opt/lalsuite --enable-swig-python --disable-lalstochastic --disable-lalxml --disable-lalinference --disable-laldetchar --disable-lalapps 2>&1 | grep -v checking
  make -j 2 2>&1 | grep Entering
  make install 2>&1 | grep Entering
  echo 'source ${VIRTUAL_ENV}/opt/lalsuite/etc/lalsuite-user-env.sh' >> ${VIRTUAL_ENV}/bin/activate
  deactivate

  source ${VENV_PATH}/bin/activate
  cd $VIRTUAL_ENV/src/lalsuite/lalapps
  LIBS="-lhdf5_hl -lhdf5 -ldl -lz" ./configure --prefix=${VIRTUAL_ENV}/opt/lalsuite --enable-static-binaries --disable-lalinference --disable-lalburst --disable-lalpulsar --disable-lalstochastic
  cd $VIRTUAL_ENV/src/lalsuite/lalapps/src/lalapps
  make -j 2
  cd $VIRTUAL_ENV/src/lalsuite/lalapps/src/inspiral
  make lalapps_inspinj
  cp lalapps_inspinj $VIRTUAL_ENV/bin
  cd $VIRTUAL_ENV/src/lalsuite/lalapps/src/ring
  make lalapps_coh_PTF_inspiral
  cp lalapps_coh_PTF_inspiral $VIRTUAL_ENV/bin

  pip install http://download.pegasus.isi.edu/pegasus/4.7.4/pegasus-python-source-4.7.4.tar.gz
  pip install git+https://github.com/ligovirgo/dqsegdb@clean_pip_install_1_4_1#egg=dqsegdb

  cd /pycbc
  python setup.py install

  pip install pycbc-pylal

  pip install 'matplotlib==1.5.3'
  pip install ipython
  pip install jupyter

cat << EOF >> $VIRTUAL_ENV/bin/activate

# If a suitable MKL exists, set it up
if [ -f /opt/intel/composer_xe_2015/mkl/bin/mklvars.sh ] ; then
  . /opt/intel/composer_xe_2015/mkl/bin/mklvars.sh intel64
elif [ -f /ldcg/intel/2017u0/compilers_and_libraries_2017.0.098/linux/mkl/bin/mklvars.sh ] ; then
  . /ldcg/intel/2017u0/compilers_and_libraries_2017.0.098/linux/mkl/bin/mklvars.sh intel64
fi

# Use the revison 11 ROM data from CVMFS
export LAL_DATA_PATH=/cvmfs/oasis.opensciencegrid.org/ligo/sw/pycbc/lalsuite-extra/11/share/lalsimulation
EOF

  deactivate

  if [ "x${TRAVIS_SECURE_ENV_VARS}" == "xtrue" ] ; then
    rsync --rsh=ssh -ravz ${CVMFS_PATH}/ pycbc@sugwg-test1.phy.syr.edu:/home/pycbc/ouser.ligo/ligo/deploy/sw/pycbc/x86_64_rhel_7/virtualenv/
    ssh pycbc@sugwg-test1.phy.syr.edu "ls -al /home/pycbc/ouser.ligo/ligo/deploy/sw/pycbc/x86_64_rhel_7/virtualenv"
  fi
fi 

echo -e "\\n>> [`date`] CentOS Docker script exiting"

exit 0
