#!/bin/bash

set -e

# Please don't try and run this script directly!
ssh -o StrictHostKeyChecking=no cvmfs.pycbc@cvmfs-software.ligo.caltech.edu <<EOF
sudo -u repo.software cvmfs_server transaction -t 300 software.igwn.org
mkdir -p /cvmfs/software.igwn.org/pycbc/${ENV_OS}/virtualenv/pycbc-${SOURCE_TAG}
EOF
rsync --rsh="ssh -o StrictHostKeyChecking=no" $RSYNC_OPTIONS -qraz ${VENV_PATH}/ cvmfs.pycbc@cvmfs-software.ligo.caltech.edu:/cvmfs/software.igwn.org/pycbc/${ENV_OS}/virtualenv/pycbc-${SOURCE_TAG}/
ssh -o StrictHostKeyChecking=no cvmfs.pycbc@cvmfs-software.ligo.caltech.edu "sudo -u repo.software cvmfs_server publish"
