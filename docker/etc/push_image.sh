docker login -u="$DOCKER_USERNAME" -p="$DOCKER_PASSWORD"
if [ "x${SOURCE_TAG}"  == "x" ] ; then
    docker tag ${DOCKER_IMG} ${DOCKER_IMG}:latest || exit 1 ;
    docker push ${DOCKER_IMG}:latest || exit 1 ;
else
    docker tag ${DOCKER_IMG} ${DOCKER_IMG}:${SOURCE_TAG} || exit 1 ;
    docker push ${DOCKER_IMG}:${SOURCE_TAG} || exit 1 ;
fi ;
