FROM igwn/base:el8

ADD docker/etc/profile.d/pycbc.sh /etc/profile.d/pycbc.sh
ADD docker/etc/profile.d/pycbc.csh /etc/profile.d/pycbc.csh

# We are going to use pip as root a lot; tell it not to complain
ENV PIP_ROOT_USER_ACTION=ignore

# Set up extra repositories
RUN <<EOF
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
python -m pip install mkl lalsuite
dnf clean all
python -m pip cache purge
EOF

# set up environment
RUN <<EOF
mkdir -p /oasis /scratch /projects /usr/lib64/slurm /var/run/munge
groupadd -g 1000 pycbc
useradd -u 1000 -g 1000 -d /opt/pycbc -k /etc/skel -m -s /bin/bash pycbc
EOF

# Now update all of our library installations
RUN rm -f /etc/ld.so.cache && /sbin/ldconfig

# Explicitly set the path so that it is not inherited from build the environment
ENV PATH="/usr/local/bin:/usr/bin:/bin:/lib64/openmpi/bin/bin"

# When the container is started with
#   docker run -it pycbc/pycbc-el8:latest
# the default is to start a login shell as the pycbc user.
# This can be overridden to log in as root with
#   docker run -it pycbc/pycbc-el8:latest /bin/bash -l
CMD ["/bin/su", "-l", "pycbc"]
