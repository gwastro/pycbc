#!/bin/bash

set -e

echo -e "\\n>> [`date`] Starting docker container ${DOCKER_IMG}"

if [ "x${DOCKER_SECURE_ENV_VARS}" == "xtrue" ] ; then
  cp -R ~/.ssh .
fi

SOURCE_TAG=`git show-ref | grep ${GITHUB_SHA} | grep -E -o "refs/tags.{0,100}" | cut -c11-`

sudo docker run --name buildvm -v `pwd`:/pycbc:rw ${DOCKER_IMG} /bin/bash -c "bash /pycbc/tools/docker_build_dist.sh --pycbc-container=${PYCBC_CONTAINER} --pycbc-code=${GITHUB_REF} --secure=${DOCKER_SECURE_ENV_VARS} --tag=${SOURCE_TAG}"

echo -e "\\n>> [`date`] Docker exited"

exit 0
