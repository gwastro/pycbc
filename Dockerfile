FROM igwn/base:el8

ADD docker/etc/profile.d/pycbc.sh /etc/profile.d/pycbc.sh
ADD docker/etc/profile.d/pycbc.csh /etc/profile.d/pycbc.csh
ADD docker/etc/cvmfs/default.local /etc/cvmfs/default.local
ADD docker/etc/cvmfs/60-osg.conf /etc/cvmfs/60-osg.conf
ADD docker/etc/cvmfs/config-osg.opensciencegrid.org.conf /etc/cvmfs/config-osg.opensciencegrid.org.conf

# Set up extra repositories
RUN <<EOF
dnf -y install https://ecsft.cern.ch/dist/cvmfs/cvmfs-release/cvmfs-release-latest.noarch.rpm
dnf -y install cvmfs cvmfs-config-default
dnf makecache
dnf -y groupinstall "Development Tools" "Scientific Support"
rpm -e --nodeps git perl-Git
dnf -y install \
    fftw \
    fftw-devel \
    fftw-libs \
    fftw-libs-double \
    fftw-libs-long \
    fftw-libs-single \
    gsl \
    gsl-devel \
    hdf5 \
    hdf5-devel \
    libjpeg-devel \
    libpng-devel \
    openssl-devel \
    osg-ca-certs \
    python3.11 \
    python3.11-devel \
    python3.11-pip \
    rsync \
    sqlite-devel \
    swig \
    which \
    zlib-devel
alternatives --set python /usr/bin/python3.11
python -m pip install --upgrade pip setuptools wheel cython
python -m pip install mkl ipython jupyter jupyterhub jupyterlab lalsuite
dnf -y install https://repo.opensciencegrid.org/osg/3.5/el8/testing/x86_64/osg-wn-client-3.5-5.osg35.el8.noarch.rpm
EOF

# set up environment
RUN <<EOF
mkdir -p /cvmfs/config-osg.opensciencegrid.org /cvmfs/software.igwn.org /cvmfs/gwosc.osgstorage.org
echo "config-osg.opensciencegrid.org /cvmfs/config-osg.opensciencegrid.org cvmfs ro,noauto 0 0" >> /etc/fstab
echo "software.igwn.org /cvmfs/software.igwn.org cvmfs ro,noauto 0 0" >> /etc/fstab
echo "gwosc.osgstorage.org /cvmfs/gwosc.osgstorage.org cvmfs ro,noauto 0 0" >> /etc/fstab
mkdir -p /oasis /scratch /projects /usr/lib64/slurm /var/run/munge
groupadd -g 1000 pycbc
useradd -u 1000 -g 1000 -d /opt/pycbc -k /etc/skel -m -s /bin/bash pycbc
EOF

# Install MPI software needed for pycbc_inference at the end.
RUN <<EOF
dnf -y install \
    libibmad \
    libibmad-devel \
    libibumad \
    libibumad-devel \
    libibverbs \
    libibverbs-devel \
    libmlx4 \
    libmlx5 \
    librdmacm \
    librdmacm-devel \
    openmpi \
    openmpi-devel
python -m pip install schwimmbad
MPICC=/lib64/openmpi/bin/mpicc CFLAGS='-I /usr/include/openmpi-x86_64/ -L /usr/lib64/openmpi/lib/ -lmpi' python -m pip install --no-cache-dir mpi4py
echo "/usr/lib64/openmpi/lib/" > /etc/ld.so.conf.d/openmpi.conf
EOF

# Now update all of our library installations
RUN rm -f /etc/ld.so.cache && /sbin/ldconfig

# Explicitly set the path so that it is not inherited from build the environment
ENV PATH "/usr/local/bin:/usr/bin:/bin:/lib64/openmpi/bin/bin"

# Set the default LAL_DATA_PATH to point at CVMFS first, then the container.
# Users wanting it to point elsewhere should start docker using:
#   docker <cmd> -e LAL_DATA_PATH="/my/new/path"
ENV LAL_DATA_PATH "/cvmfs/software.igwn.org/pycbc/lalsuite-extra/current/share/lalsimulation:/opt/pycbc/pycbc-software/share/lal-data"

# When the container is started with
#   docker run -it pycbc/pycbc-el8:latest
# the default is to start a login shell as the pycbc user.
# This can be overridden to log in as root with
#   docker run -it pycbc/pycbc-el8:latest /bin/bash -l
CMD ["/bin/su", "-l", "pycbc"]
