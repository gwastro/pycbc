# Find tag (if any) corresponding to the commit being tested
SOURCE_TAG=`git tag --points-at ${GITHUB_SHA}`
echo "Source tag: ${SOURCE_TAG}"

# Find latest hash of upstream master branch
MASTER_HASH=`git ls-remote origin master | cut -f 1`
echo "Master hash: ${MASTER_HASH}"

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
