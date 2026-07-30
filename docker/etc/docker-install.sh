#!/bin/bash -v

set -e

echo //// DEBUG ///////////////////////////////////////
python -c "import sys; print(sys.executable); print(sys.path); import traceback; import importlib; print('stdlib OK')"
echo //////////////////////////////////////////////////

# Install PyCBC
cd /scratch
python -m pip install --verbose --upgrade pip
python -m pip install --verbose -r requirements.txt
echo //// DEBUG ///////////////////////////////////////
python -c "import sys; print(sys.executable); print(sys.path); import traceback; import importlib; print('stdlib OK')"
echo //////////////////////////////////////////////////
export PYTHONVERBOSE=true
python -m pip install --verbose -r requirements-igwn.txt
python -m pip install -r companion.txt
python -m pip install .
cd /

# Copy PyCBC source repository into the image
mkdir -p /opt/pycbc/src
cp -a /scratch /opt/pycbc/src/pycbc

# Copy the lalsuite extra data into the image
chmod a+rw /dev/fuse
mount /cvmfs/config-osg.opensciencegrid.org
mount /cvmfs/software.igwn.org
mkdir -p /opt/pycbc/pycbc-software/share/lal-data
rsync \
	--exclude='SEOBNRv1ROM*' \
	--exclude='SEOBNRv2ROM_DS_HI_v1.0.hdf5' \
	--exclude='NRSur4d2s_FDROM.hdf5' \
	-ravz \
	/cvmfs/software.igwn.org/pycbc/lalsuite-extra/current/share/lalsimulation/ \
	/opt/pycbc/pycbc-software/share/lal-data/
umount /cvmfs/config-osg.opensciencegrid.org
umount /cvmfs/software.igwn.org

# Set ownership and permissions correctly
chown -R 1000:1000 /opt/pycbc/src /opt/pycbc/pycbc-software
chmod -R u=rwX,g=rX,o=rX /opt/pycbc

# Cleanup any leftover stuff cached by dnf and pip
dnf clean all
python -m pip cache purge

exit 0
