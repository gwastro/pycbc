#!/bin/bash -v
set -x
id
modprobe fuse
cd /scratch
pip install -r requirements.txt
pip install .
cd /
mkdir -p /opt/pycbc/pycbc-software/share/lal-data
mount /cvmfs/config-osg.opensciencegrid.org
mount /cvmfs/oasis.opensciencegrid.org
rsync --exclude='SEOBNRv1ROM*' --exclude='SEOBNRv2ROM_DS_HI_v1.0.hdf5' --exclude='NRSur4d2s_FDROM.hdf5' -ravz /cvmfs/oasis.opensciencegrid.org/ligo/sw/pycbc/lalsuite-extra/current/share/lalsimulation/ /opt/pycbc/pycbc-software/share/lal-data/
umount /cvmfs/config-osg.opensciencegrid.org
umount /cvmfs/oasis.opensciencegrid.org
exit 0
