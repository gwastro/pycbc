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
LALSUITE_HASH="539c8700af92eb6dd00e0e91b9dbaf5bae51f004"

if [ "x$TRAVIS_TAG" == "x" ] ; then
  TRAVIS_TAG="latest"
  RSYNC_OPTIONS="--delete"
else
  RSYNC_OPTIONS=""
fi

echo -e "\\n>> [`date`] Inside container ${PYCBC_CONTAINER}"
echo -e "\\n>> [`date`] Release tag is ${TRAVIS_TAG}"
echo -e "\\n>> [`date`] Using PyCBC code ${PYCBC_CODE}"
echo -e "\\n>> [`date`] Using lalsuite hash ${LALSUITE_HASH}"
echo -e "\\n>> [`date`] Travis pull request is ${TRAVIS_PULL_REQUEST}"
echo -e "\\n>> [`date`] Travis commit is ${TRAVIS_COMMIT}"
echo -e "\\n>> [`date`] Travis secure env is ${TRAVIS_SECURE_ENV_VARS}"
echo -e "\\n>> [`date`] Travis tag is ${TRAVIS_TAG}"

if [ "x${TRAVIS_SECURE_ENV_VARS}" == "xtrue" ] ; then
  mkdir -p ~/.ssh
  cp /pycbc/.ssh/* ~/.ssh
  chmod 600 ~/.ssh/id_rsa
fi

if [ "x${PYCBC_CONTAINER}" == "xpycbc_inspiral_bundle" ] ; then
  echo -e "\\n>> [`date`] Building pycbc_inspiral bundle for CentOS 6"

  # create working dir for build script
  BUILD=/pycbc/build
  mkdir -p ${BUILD}
  export PYTHONUSERBASE=${BUILD}/.local
  export XDG_CACHE_HOME=${BUILD}/.cache

  # get library to build optimized pycbc_inspiral bundle
  wget_opts="-c --passive-ftp --no-check-certificate --tries=5 --timeout=30 --no-verbose"
  primary_url="https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-software/download/03f2048c770492f66f80528493fd6cecded63769"
  secondary_url="https://www.atlas.aei.uni-hannover.de/~dbrown"
  pushd /pycbc
  for p in "x86_64/composer_xe_2015.0.090/composer_xe_2015.0.090.tar.gz" "bank-files/testbank_TF2v4ROM.hdf" ; do
    set +e
    test -r `basename $p` || wget $wget_opts ${primary_url}/${p}
    set -e
    test -r `basename $p` || wget $wget_opts ${secondary_url}/${p}
  done
  popd

  # run the einstein at home build and test script
  echo -e "\\n>> [`date`] Running pycbc_build_eah.sh"
  pushd ${BUILD}
  /pycbc/tools/einsteinathome/pycbc_build_eah.sh --lalsuite-commit=${LALSUITE_HASH} ${PYCBC_CODE} --clean-pycbc --silent-build --with-extra-approximant='SPAtmplt:mtotal<4' --with-extra-approximant='SEOBNRv4_ROM:else'  --with-extra-approximant=--use-compressed-waveforms --with-extra-libs=file:///pycbc/composer_xe_2015.0.090.tar.gz --processing-scheme=mkl

  if [ "x${TRAVIS_SECURE_ENV_VARS}" == "xtrue" ] ; then
    echo -e "\\n>> [`date`] Deploying pycbc_inspiral bundle"
    BUNDLE_DEST=/home/pycbc/ouser.ligo/ligo/deploy/sw/pycbc/x86_64_rhel_6/bundle/${TRAVIS_TAG}
    echo -e "\\n>> [`date`] Deploying pycbc_inspiral bundle to sugwg-condor.phy.syr.edu"
    ssh pycbc@sugwg-condor.phy.syr.edu "mkdir -p ${BUNDLE_DEST}"
    scp ${BUILD}/pycbc-build/environment/dist/pycbc_inspiral_osg* pycbc@sugwg-condor.phy.syr.edu:${BUNDLE_DEST}/pycbc_inspiral
    if [ "x${TRAVIS_TAG}" == "xlatest" ] ; then
      PYCBC_INSPIRAL_SUFFIX="_osg_${TRAVIS_TAG}"
      BUNDLE_DEST=/home/login/ouser.ligo/ligo/deploy/sw/pycbc/x86_64_rhel_6/bundle/${TRAVIS_TAG}
      echo -e "\\n>> [`date`] Deploying pycbc_inspiral${PYCBC_INSPIRAL_SUFFIX} to CVMFS"
      ssh ouser.ligo@oasis-login.opensciencegrid.org "mkdir -p ${BUNDLE_DEST}"
      scp ${BUILD}/pycbc-build/environment/dist/pycbc_inspiral${PYCBC_INSPIRAL_SUFFIX} ouser.ligo@oasis-login.opensciencegrid.org:${BUNDLE_DEST}/pycbc_inspiral
      ssh ouser.ligo@oasis-login.opensciencegrid.org osg-oasis-update
    fi
    echo -e "\\n>> [`date`] pycbc_inspiral deployment complete"
  fi
  popd
fi

if [ "x${PYCBC_CONTAINER}" == "xpycbc_virtualenv" ] ; then
  echo -e "\\n>> [`date`] Building pycbc virtual environment for CentOS 7"

  echo -e "\\n>> [`date`] Removing LAL RPMs"
  yum -y -q remove "*lal*"

  CVMFS_PATH=/cvmfs/oasis.opensciencegrid.org/ligo/sw/pycbc/x86_64_rhel_7/virtualenv
  mkdir -p ${CVMFS_PATH}

  VENV_PATH=${CVMFS_PATH}/pycbc-${TRAVIS_TAG}
  pip install virtualenv
  virtualenv ${VENV_PATH}
  echo 'export PYTHONUSERBASE=${VIRTUAL_ENV}/.local' >> ${VENV_PATH}/bin/activate
  echo "export XDG_CACHE_HOME=\${HOME}/cvmfs-pycbc-${TRAVIS_TAG}/.cache" >> ${VENV_PATH}/bin/activate
  source ${VENV_PATH}/bin/activate
  mkdir -p ${VIRTUAL_ENV}/.local
  echo -e "[easy_install]\\nzip_ok = false\\n" > ~/.pydistutils.cfg
  echo -e "[easy_install]\\nzip_ok = false\\n" > ${VIRTUAL_ENV}/.local/.pydistutils.cfg
  
  echo -e "\\n>> [`date`] Upgrading pip and setuptools"
  pip install --upgrade pip
  pip install six packaging appdirs
  pip install --upgrade setuptools

  echo -e "\\n>> [`date`] Installing base python packages required to build lalsuite"
  pip install "numpy>=1.6.4" "h5py>=2.5" unittest2 python-cjson Cython decorator
  echo -e "\\n>> [`date`] Installing scipy"
  pip install "scipy>=0.13.0" &>/dev/null
  echo -e "\\n>> [`date`] Installing M2Crypto"
  SWIG_FEATURES="-cpperraswarn -includeall -I/usr/include/openssl" pip install M2Crypto

  echo -e "\\n>> [`date`] Installing LAL"
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

  echo -e "\\n>> [`date`] Installing LALApps"
  source ${VENV_PATH}/bin/activate
  cd $VIRTUAL_ENV/src/lalsuite/lalapps
  LIBS="-lhdf5_hl -lhdf5 -ldl -lz" ./configure --prefix=${VIRTUAL_ENV}/opt/lalsuite --enable-static-binaries --disable-lalinference --disable-lalburst --disable-lalpulsar --disable-lalstochastic 2>&1 | grep -v checking
  cd $VIRTUAL_ENV/src/lalsuite/lalapps/src/lalapps
  make -j 2 2>&1 | grep Entering
  cd $VIRTUAL_ENV/src/lalsuite/lalapps/src/inspiral
  make lalapps_inspinj
  cp lalapps_inspinj $VIRTUAL_ENV/bin
  cd $VIRTUAL_ENV/src/lalsuite/lalapps/src/ring
  make lalapps_coh_PTF_inspiral
  cp lalapps_coh_PTF_inspiral $VIRTUAL_ENV/bin

  echo -e "\\n>> [`date`] Installing Pegasus and DQSegDB"
  pip install http://download.pegasus.isi.edu/pegasus/4.7.4/pegasus-python-source-4.7.4.tar.gz
  pip install dqsegdb

  echo -e "\\n>> [`date`] Install matplotlib 1.5.3"
  pip install 'matplotlib==1.5.3'

  echo -e "\\n>> [`date`] Installing PyCBC and dependencies"
  cd /pycbc
  python setup.py install

  echo -e "\\n>> [`date`] Installing PyCBC PyLAL"
  pip install pycbc-pylal

  echo -e "\\n>> [`date`] Installing modules needed to build documentation"
  pip install "Sphinx>=1.5.0"
  pip install sphinx-rtd-theme
  pip install git+https://github.com/ligo-cbc/sphinxcontrib-programoutput.git#egg=sphinxcontrib-programoutput

  echo -e "\\n>> [`date`] Installing ipython and jupyter"
  pip install ipython
  pip install jupyter

cat << EOF >> $VIRTUAL_ENV/bin/activate

# if a suitable MKL exists, set it up
if [ -f /opt/intel/composer_xe_2015/mkl/bin/mklvars.sh ] ; then
  # location on syracuse cluster
  . /opt/intel/composer_xe_2015/mkl/bin/mklvars.sh intel64
elif [ -f /opt/intel/2015/composer_xe_2015/mkl/bin/mklvars.sh ] ; then
  # location on atlas cluster
  . /opt/intel/2015/composer_xe_2015/mkl/bin/mklvars.sh
elif [ -f /ldcg/intel/2017u0/compilers_and_libraries_2017.0.098/linux/mkl/bin/mklvars.sh ] ; then
  # location on cit cluster
  . /ldcg/intel/2017u0/compilers_and_libraries_2017.0.098/linux/mkl/bin/mklvars.sh intel64
fi

# Use the revison 11 ROM data from CVMFS
export LAL_DATA_PATH=/cvmfs/oasis.opensciencegrid.org/ligo/sw/pycbc/lalsuite-extra/11/share/lalsimulation
EOF

  deactivate

  if [ "x${TRAVIS_SECURE_ENV_VARS}" == "xtrue" ] ; then
    echo -e "\\n>> [`date`] Running test_coinc_search_workflow.sh"
    mkdir -p /pycbc/workflow-test
    pushd /pycbc/workflow-test
    /pycbc/tools/test_coinc_search_workflow.sh ${TRAVIS_TAG}
    popd
  fi

  if [ "x${TRAVIS_SECURE_ENV_VARS}" == "xtrue" ] ; then
    echo -e "\\n>> [`date`] Setting virtual environment permissions for deployment"
    find ${VENV_PATH} -type d -exec chmod go+rx {} \;
    chmod -R go+r ${VENV_PATH}
  
    echo -e "\\n>> [`date`] Deploying virtual environment ${VENV_PATH}"
    echo -e "\\n>> [`date`] Deploying virtual environment to sugwg-condor.phy.syr.edu"
    ssh pycbc@sugwg-condor.phy.syr.edu "mkdir -p /home/pycbc/ouser.ligo/ligo/deploy/sw/pycbc/x86_64_rhel_7/virtualenv/pycbc-${TRAVIS_TAG}"
    rsync --rsh=ssh $RSYNC_OPTIONS -qraz ${VENV_PATH}/ pycbc@sugwg-condor.phy.syr.edu:/home/pycbc/ouser.ligo/ligo/deploy/sw/pycbc/x86_64_rhel_7/virtualenv/pycbc-${TRAVIS_TAG}/
    if [ "x${TRAVIS_TAG}" == "xlatest" ] ; then
      echo -e "\\n>> [`date`] Deploying release ${TRAVIS_TAG} to CVMFS"
      # remove lalsuite source and deploy on cvmfs
      rm -rf ${VENV_PATH}/src/lalsuite
      ssh ouser.ligo@oasis-login.opensciencegrid.org "mkdir -p /home/login/ouser.ligo/ligo/deploy/sw/pycbc/x86_64_rhel_7/virtualenv/pycbc-${TRAVIS_TAG}"
      rsync --rsh=ssh $RSYNC_OPTIONS -qraz ${VENV_PATH}/ ouser.ligo@oasis-login.opensciencegrid.org:/home/login/ouser.ligo/ligo/deploy/sw/pycbc/x86_64_rhel_7/virtualenv/pycbc-${TRAVIS_TAG}/
      ssh ouser.ligo@oasis-login.opensciencegrid.org osg-oasis-update
    fi
    echo -e "\\n>> [`date`] virtualenv deployment complete"
  fi
fi 

echo -e "\\n>> [`date`] CentOS Docker script exiting"

exit 0
