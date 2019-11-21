#!/bin/bash -v
set -e
cd /scratch
pip install -r requirements.txt
pip install .
cd /
mkdir -p /opt/pycbc/src
cp -a /scratch /opt/pycbc/src/pycbc
mount /cvmfs/config-osg.opensciencegrid.org
mount /cvmfs/oasis.opensciencegrid.org
mkdir -p /opt/pycbc/pycbc-software/share/lal-data
rsync --exclude='SEOBNRv1ROM*' --exclude='SEOBNRv2ROM_DS_HI_v1.0.hdf5' --exclude='NRSur4d2s_FDROM.hdf5' -ravz /cvmfs/oasis.opensciencegrid.org/ligo/sw/pycbc/lalsuite-extra/current/share/lalsimulation/ /opt/pycbc/pycbc-software/share/lal-data/
umount /cvmfs/config-osg.opensciencegrid.org
umount /cvmfs/oasis.opensciencegrid.org
chown -R 1000:1000 /opt/pycbc/src /opt/pycbc/pycbc-software
chmod -R u=rwX,g=rX,o=rX /opt/pycbc
exit 0
