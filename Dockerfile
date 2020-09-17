FROM centos:centos7

ADD docker/etc/profile.d/pycbc.sh /etc/profile.d/pycbc.sh
ADD docker/etc/profile.d/pycbc.csh /etc/profile.d/pycbc.csh
ADD docker/etc/cvmfs/default.local /etc/cvmfs/default.local
ADD docker/etc/cvmfs/60-osg.conf /etc/cvmfs/60-osg.conf
ADD docker/etc/cvmfs/config-osg.opensciencegrid.org.conf /etc/cvmfs/config-osg.opensciencegrid.org.conf

# Set up extra repositories
RUN rpm -ivh http://software.ligo.org/lscsoft/scientific/7/x86_64/production/l/lscsoft-production-config-1.3-1.el7.noarch.rpm && yum install -y https://repo.ius.io/ius-release-el7.rpm && yum install -y https://ecsft.cern.ch/dist/cvmfs/cvmfs-release/cvmfs-release-latest.noarch.rpm && yum install -y cvmfs cvmfs-config-default && yum -y install http://repo.opensciencegrid.org/osg/3.4/osg-3.4-el7-release-latest.rpm && yum clean all && yum makecache && \
    yum -y groupinstall "Compatibility Libraries" \
                        "Development Tools" \
                        "Scientific Support" && \
    rpm -e --nodeps git perl-Git && yum -y install python2-pip python-setuptools zlib-devel libpng-devel libjpeg-devel libsqlite3-dev sqlite-devel db4-devel openssl-devel git2u-all fftw-libs-single fftw-devel fftw fftw-libs-long fftw-libs fftw-libs-double gsl gsl-devel libframe-utils libframe-devel libframe libmetaio libmetaio-devel libmetaio-utils hdf5 hdf5-devel python-devel which osg-wn-client osg-ca-certs && pip install --upgrade pip==19.3.1 setuptools==44.0.0 && pip install mkl ipython jupyter

# set up environment
RUN cd / && \
    mkdir -p /cvmfs/config-osg.opensciencegrid.org /cvmfs/oasis.opensciencegrid.org /cvmfs/gwosc.osgstorage.org && echo "config-osg.opensciencegrid.org /cvmfs/config-osg.opensciencegrid.org cvmfs ro,noauto 0 0" >> /etc/fstab && echo "oasis.opensciencegrid.org /cvmfs/oasis.opensciencegrid.org cvmfs ro,noauto 0 0" >> /etc/fstab && echo "gwosc.osgstorage.org /cvmfs/gwosc.osgstorage.org cvmfs ro,noauto 0 0" >> /etc/fstab && mkdir -p /oasis /scratch /projects /usr/lib64/slurm /var/run/munge && \
    groupadd -g 1000 pycbc && useradd -u 1000 -g 1000 -d /opt/pycbc -k /etc/skel -m -s /bin/bash pycbc

# Install MPI software needed for pycbc_inference
RUN yum install -y libibverbs libibverbs-devel libibmad libibmad-devel libibumad libibumad-devel librdmacm librdmacm-devel libmlx5 libmlx4 && curl http://mvapich.cse.ohio-state.edu/download/mvapich/mv2/mvapich2-2.1.tar.gz | tar -C /tmp -zxf - && \
    cd /tmp/mvapich2-2.1 && ./configure --prefix=/opt/mvapich2-2.1 && make install && \
    cd / && rm -rf /tmp/mvapich2-2.1 && \
    pip install schwimmbad && \
    MPICC=/opt/mvapich2-2.1/bin CFLAGS='-I /opt/mvapich2-2.1/include -L /opt/mvapich2-2.1/lib -lmpi' pip install --no-cache-dir mpi4py
RUN echo "/opt/mvapich2-2.1/lib" > /etc/ld.so.conf.d/mvapich2-2.1.conf

# Now update all of our library installations
RUN rm -f /etc/ld.so.cache && /sbin/ldconfig

# Explicitly set the path so that it is not inherited from build the environment
ENV PATH "/usr/local/bin:/usr/bin:/bin:/opt/mvapich2-2.1/bin"

# Set the default LAL_DATA_PATH to point at CVMFS first, then the container.
# Users wanting it to point elsewhere should start docker using:
#   docker <cmd> -e LAL_DATA_PATH="/my/new/path"
ENV LAL_DATA_PATH "/cvmfs/oasis.opensciencegrid.org/ligo/sw/pycbc/lalsuite-extra/current/share/lalsimulation:/opt/pycbc/pycbc-software/share/lal-data"

# When the container is started with
#   docker run -it pycbc/pycbc-el7:latest
# the default is to start a loging shell as the pycbc user.
# This can be overridden to log in as root with
#   docker run -it pycbc/pycbc-el7:latest /bin/bash -l
CMD ["/bin/su", "-l", "pycbc"]
