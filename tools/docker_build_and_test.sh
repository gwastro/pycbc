#!/bin/bash

set -e

echo -e "\\n>> [`date`] Starting docker container ${DOCKER_IMG}"

if [ "x${TRAVIS_SECURE_ENV_VARS}" == "xtrue" ] ; then
  cp -R ~/.ssh .
fi

if [ "x${PYCBC_CONTAINER}" == "xpycbc_inspiral_bundle" ] ; then
  echo -e "\\n>> [`date`] Looking for cached E@H build environment"
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
fi

sudo docker run --name buildvm -v `pwd`:/pycbc:rw ${DOCKER_IMG} /bin/bash -c "bash /pycbc/tools/docker_build_dist.sh --pycbc-container=${PYCBC_CONTAINER} --pull-request=${TRAVIS_PULL_REQUEST} --commit=${TRAVIS_COMMIT} --secure=${TRAVIS_SECURE_ENV_VARS} --tag=${TRAVIS_TAG}"

echo -e "\\n>> [`date`] Docker exited"

if [ "x${PYCBC_CONTAINER}" == "xpycbc_inspiral_bundle" ] ; then
  echo -e "\\n>> [`date`] Caching E@H build environment"
  mkdir -p $HOME/docker-cache
  mkdir -p $HOME/docker-cache/test
  sudo docker cp buildvm:/pycbc/build/pycbc-sources/pycbc-build-preinst.tgz $HOME/docker-cache/pycbc-build-preinst.tgz
  sudo docker cp buildvm:/pycbc/build/pycbc-sources/pycbc-build-preinst-lalsuite.tgz $HOME/docker-cache/pycbc-build-preinst-lalsuite.tgz
  sudo docker cp buildvm:/pycbc/build/pycbc-sources/test/H1L1-SBANK_FOR_GW150914ER10.xml.gz $HOME/docker-cache/test/H1L1-SBANK_FOR_GW150914ER10.xml.gz
  sudo docker cp buildvm:/pycbc/build/pycbc-sources/test/H-H1_LOSC_4_V1-1126257414-4096.gwf $HOME/docker-cache/test/H-H1_LOSC_4_V1-1126257414-4096.gwf
fi

exit 0
