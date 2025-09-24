SOURCE_TAG=`git show-ref | grep ${GITHUB_SHA} | grep -E -o "refs/tags.{0,100}" | cut -c11-`
MASTER_HASH=`git rev-parse origin/master`

if [ "x${SOURCE_TAG}" == "x" ] ; then
    if [ "x${MASTER_HASH}" != "x${GITHUB_SHA}" ] ; then
        echo "Not on master or a tag, so not pushing the docker image."
        exit
    fi
    SOURCE_TAG=latest
fi

docker login --username "$DOCKER_USERNAME" --password "$DOCKER_PASSWORD"
docker tag ${DOCKER_IMG} ${DOCKER_IMG}:${SOURCE_TAG} || exit 1 ;
docker push ${DOCKER_IMG}:${SOURCE_TAG} || exit 1 ;
