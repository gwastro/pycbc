#!/bin/bash

set -e

echo -e "\\n>> [`date`] Starting docker container ${DOCKER_IMG}"

if [ "x${TRAVIS_SECURE_ENV_VARS}" == "xtrue" ] ; then
  cp -R ~/.ssh .
fi

sudo docker run --name buildvm -v `pwd`:/pycbc:rw ${DOCKER_IMG} /bin/bash -c "bash /pycbc/tools/docker_build_dist.sh --pycbc-container=${PYCBC_CONTAINER} --pull-request=${TRAVIS_PULL_REQUEST} --commit=${TRAVIS_COMMIT} --secure=${TRAVIS_SECURE_ENV_VARS} --tag=${TRAVIS_TAG}"

echo -e "\\n>> [`date`] Docker exited"

exit 0
