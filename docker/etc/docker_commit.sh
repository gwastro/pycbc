docker commit --change='CMD ["/bin/su", "-l", "pycbc"]' pycbc_inst ${DOCKER_IMG}
