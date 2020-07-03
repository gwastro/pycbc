#!/bin/bash

set -e

for i in $*; do
  case $i in
    --pycbc-container=*) PYCBC_CONTAINER="`echo $i|sed 's/^--pycbc-container=//'`";;
    --pull-request=*) TRAVIS_PULL_REQUEST="`echo $i|sed 's/^--pull-request=//'`";;
    --commit=*) TRAVIS_COMMIT="`echo $i|sed 's/^--commit=//'`";;
    --secure=*) TRAVIS_SECURE_ENV_VARS="`echo $i|sed 's/^--secure=//'`";;
    --tag=*) TRAVIS_TAG="`echo $i|sed 's/^--tag=//'`";;
    *) echo -e "unknown option '$i', valid are:\n$usage">&2; exit 1;;
  esac
done

# determine the pycbc git branch and origin
if test x$TRAVIS_PULL_REQUEST = "xfalse" ; then
    PYCBC_CODE="--pycbc-commit=${TRAVIS_COMMIT}"
else
    PYCBC_CODE="--pycbc-fetch-ref=refs/pull/${TRAVIS_PULL_REQUEST}/merge"
fi

# set the lalsuite checkout to use

if [ "x$TRAVIS_TAG" == "x" ] ; then
  TRAVIS_TAG="master"
  RSYNC_OPTIONS="--delete"
else
  RSYNC_OPTIONS=""
fi

echo -e "\\n>> [`date`] Inside container ${PYCBC_CONTAINER}"
echo -e "\\n>> [`date`] Release tag is ${TRAVIS_TAG}"
echo -e "\\n>> [`date`] Using PyCBC code ${PYCBC_CODE}"
echo -e "\\n>> [`date`] Travis pull request is ${TRAVIS_PULL_REQUEST}"
echo -e "\\n>> [`date`] Travis commit is ${TRAVIS_COMMIT}"
echo -e "\\n>> [`date`] Travis secure env is ${TRAVIS_SECURE_ENV_VARS}"
echo -e "\\n>> [`date`] Travis tag is ${TRAVIS_TAG}"

if [ "x${TRAVIS_SECURE_ENV_VARS}" == "xtrue" ] ; then
  mkdir -p ~/.ssh
  cp /pycbc/.ssh/* ~/.ssh
  chmod 600 ~/.ssh/id_rsa
fi

if [ "x${PYCBC_CONTAINER}" == "xpycbc_rhel_virtualenv" ]; then

  ENV_OS="x86_64_rhel_7"
  yum -y install python2-pip python-setuptools which
  yum -y install curl
  curl http://download.pegasus.isi.edu/wms/download/rhel/7/pegasus.repo > /etc/yum.repos.d/pegasus.repo
  yum clean all
  yum makecache
  yum -y install openssl-devel openssl-static
  yum -y install pegasus
  yum -y install ligo-proxy-utils
  yum -y install ecp-cookie-init
  yum -y install python-virtualenv
  yum -y install hdf5-static libxml2-static zlib-static libstdc++-static cfitsio-static glibc-static fftw-static gsl-static --skip-broken

  CVMFS_PATH=/cvmfs/oasis.opensciencegrid.org/ligo/sw/pycbc/${ENV_OS}/virtualenv
  mkdir -p ${CVMFS_PATH}

  VENV_PATH=${CVMFS_PATH}/pycbc-${TRAVIS_TAG}
  virtualenv ${VENV_PATH}
  echo 'export PYTHONUSERBASE=${VIRTUAL_ENV}/.local' >> ${VENV_PATH}/bin/activate
  echo "export XDG_CACHE_HOME=\${HOME}/cvmfs-pycbc-${TRAVIS_TAG}/.cache" >> ${VENV_PATH}/bin/activate
  source ${VENV_PATH}/bin/activate
  mkdir -p ${VIRTUAL_ENV}/.local
  echo -e "[easy_install]\\nzip_ok = false\\n" > ~/.pydistutils.cfg
  echo -e "[easy_install]\\nzip_ok = false\\n" > ${VIRTUAL_ENV}/.local/.pydistutils.cfg

  echo -e "\\n>> [`date`] Upgrading pip and setuptools"
  pip install --upgrade pip setuptools
  pip install six packaging appdirs

  echo -e "\\n>> [`date`] Installing PyCBC dependencies from requirements.txt"
  cd /pycbc
  pip install -r requirements.txt
  pip install -r companion.txt

  echo -e "\\n>> [`date`] Installing PyCBC from source"
  pip install .

  echo -e "\\n>> [`date`] Installing ipython and jupyter"
  pip install jupyter

  cat << EOF >> $VIRTUAL_ENV/bin/activate

# if a suitable MKL exists, set it up
if [ -f /opt/intel/composer_xe_2015/mkl/bin/mklvars.sh ] ; then
  # location on syracuse cluster
  . /opt/intel/composer_xe_2015/mkl/bin/mklvars.sh intel64
elif [ -f /opt/intel/2015/composer_xe_2015/mkl/bin/mklvars.sh ] ; then
  # location on atlas cluster
  . /opt/intel/2015/composer_xe_2015/mkl/bin/mklvars.sh intel64
elif [ -f /ldcg/intel/2017u0/compilers_and_libraries_2017.0.098/linux/mkl/bin/mklvars.sh ] ; then
  # location on cit cluster
  . /ldcg/intel/2017u0/compilers_and_libraries_2017.0.098/linux/mkl/bin/mklvars.sh intel64
elif [ -f /apps/compilers/intel/2019.3/compilers_and_libraries/linux/mkl/bin/mklvars.sh ] ; then
  # location on ARCCA Hawk cluster
  . /apps/compilers/intel/2019.3/compilers_and_libraries/linux/mkl/bin/mklvars.sh intel64
fi

# Use the ROM data from CVMFS
export LAL_DATA_PATH=/cvmfs/oasis.opensciencegrid.org/ligo/sw/pycbc/lalsuite-extra/e02dab8c/share/lalsimulation
EOF

  deactivate

  echo -e "\\n>> [`date`] Running test_coinc_search_workflow.sh"
  mkdir -p /pycbc/workflow-test
  pushd /pycbc/workflow-test
  /pycbc/tools/test_coinc_search_workflow.sh ${VENV_PATH} ${TRAVIS_TAG}
  popd

  if [ "x${TRAVIS_SECURE_ENV_VARS}" == "xtrue" ] ; then
    echo -e "\\n>> [`date`] Setting virtual environment permissions for deployment"
    find ${VENV_PATH} -type d -exec chmod go+rx {} \;
    chmod -R go+r ${VENV_PATH}

    echo -e "\\n>> [`date`] Deploying virtual environment ${VENV_PATH}"
    echo -e "\\n>> [`date`] Deploying virtual environment to sugwg-condor.phy.syr.edu"
    ssh pycbc@sugwg-condor.phy.syr.edu "mkdir -p /home/pycbc/ouser.ligo/ligo/deploy/sw/pycbc/${ENV_OS}/virtualenv/pycbc-${TRAVIS_TAG}"
    rsync --rsh=ssh $RSYNC_OPTIONS -qraz ${VENV_PATH}/ pycbc@sugwg-condor.phy.syr.edu:/home/pycbc/ouser.ligo/ligo/deploy/sw/pycbc/${ENV_OS}/virtualenv/pycbc-${TRAVIS_TAG}/
    if [ "x${TRAVIS_TAG}" != "xmaster" ] ; then
      echo -e "\\n>> [`date`] Deploying release ${TRAVIS_TAG} to CVMFS"
      # remove lalsuite source and deploy on cvmfs
      rm -rf ${VENV_PATH}/src/lalsuite
      ssh ouser.ligo@oasis-login.opensciencegrid.org "mkdir -p /home/login/ouser.ligo/ligo/deploy/sw/pycbc/${ENV_OS}/virtualenv/pycbc-${TRAVIS_TAG}"
      rsync --rsh=ssh $RSYNC_OPTIONS -qraz ${VENV_PATH}/ ouser.ligo@oasis-login.opensciencegrid.org:/home/login/ouser.ligo/ligo/deploy/sw/pycbc/${ENV_OS}/virtualenv/pycbc-${TRAVIS_TAG}/
      ssh ouser.ligo@oasis-login.opensciencegrid.org osg-oasis-update
    fi
    echo -e "\\n>> [`date`] virtualenv deployment complete"
  fi
fi

echo -e "\\n>> [`date`] Docker script exiting"

exit 0
