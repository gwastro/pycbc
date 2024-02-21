#!/bin/bash

set -e

# Please don't try and run this script directly!
ssh cvmfs.pycbc@cvmfs-software.ligo.caltech.edu <<EOF
sudo -u repo.software cvmfs_server transaction -t 300 software.igwn.org
sudo -u repo.software chmod go+rx /cvmfs/software.igwn.org/
mkdir -p /cvmfs/software.igwn.org/pycbc/${ENV_OS}/virtualenv/pycbc-${SOURCE_TAG}
touch /cvmfs/software.igwn.org/pycbc/${ENV_OS}/virtualenv/pycbc-${SOURCE_TAG}/.cvmfscatalog
EOF
rsync --rsh="ssh" $RSYNC_OPTIONS --filter='Pp .cvmfscatalog' --filter='Pp .cvmfsautocatalog' -qraz ${VENV_PATH}/ cvmfs.pycbc@cvmfs-software.ligo.caltech.edu:/cvmfs/software.igwn.org/pycbc/${ENV_OS}/virtualenv/pycbc-${SOURCE_TAG}/
ssh cvmfs.pycbc@cvmfs-software.ligo.caltech.edu "sudo -u repo.software cvmfs_server publish"
