echo "$DOCKER_PASSWORD" | docker login --username "$DOCKER_USERNAME" --password-stdin
SOURCE_TAG=`git show-ref | grep "8d86765d0b3d059dd6e5c6195f00f2cf252e661f" | grep -E -o "refs/tags.{0,100}" | cut -c11-`
if [ "x${SOURCE_TAG}"  == "x" ] ; then
    docker tag ${DOCKER_IMG} ${DOCKER_IMG}:latest || exit 1 ;
    docker push ${DOCKER_IMG}:latest || exit 1 ;
else
    docker tag ${DOCKER_IMG} ${DOCKER_IMG}:${SOURCE_TAG} || exit 1 ;
    docker push ${DOCKER_IMG}:${SOURCE_TAG} || exit 1 ;
fi ;
