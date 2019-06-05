FROM centos:centos7

# Set up extra repositories
RUN rpm -ivh http://software.ligo.org/lscsoft/scientific/7/x86_64/production/lscsoft-production-config-1.3-1.el7.noarch.rpm
RUN yum install -y https://centos7.iuscommunity.org/ius-release.rpm
RUN yum clean all
RUN yum makecache

# Install tool sets
RUN yum -y groupinstall "Compatibility Libraries" \
                        "Development Tools" \
                        "Scientific Support"

# Update pip and setuptools
RUN yum -y install python2-pip python-setuptools

# Add extra development libraries
RUN yum -y install zlib-devel libpng-devel libjpeg-devel libsqlite3-dev sqlite-devel db4-devel openssl-devel

# Install full git
RUN rpm -e --nodeps git perl-Git
RUN yum -y install git2u-all

# Add numerical libraries needed for lal
RUN yum install -y fftw-libs-single fftw-devel fftw fftw-libs-long fftw-libs fftw-libs-double
RUN yum install -y gsl gsl-devel

# Add LIGO libraries needed for lal
RUN yum install -y libframe-utils libframe-devel libframe
RUN yum install -y libmetaio libmetaio-devel libmetaio-utils

# Add HDF5 development tools for PyCBC
RUN yum install -y hdf5 hdf5-devel

# Install MKL
RUN pip install mkl

# Install Python development tools
RUN yum install -y python-devel which

# Upgrade to latest python installers
RUN pip install --upgrade pip setuptools

# Install PyCBC dependencies
RUN pip install ipython jupyter

# set up environment
ADD docker/etc/profile.d/pycbc.sh /etc/profile.d/pycbc.sh
ADD docker/etc/profile.d/pycbc.csh /etc/profile.d/pycbc.csh

# add singularity profiles
COPY docker/.singularity.d /.singularity.d
RUN cd / && \
    ln -s .singularity.d/actions/exec .exec && \
    ln -s .singularity.d/actions/run .run && \
    ln -s .singularity.d/actions/test .shell && \
    ln -s .singularity.d/runscript singularity

# Install MPI software needed for pycbc_inference
RUN yum install -y libibverbs libibverbs-devel libibmad libibmad-devel libibumad libibumad-devel librdmacm librdmacm-devel libmlx5 libmlx4
RUN curl http://mvapich.cse.ohio-state.edu/download/mvapich/mv2/mvapich2-2.1.tar.gz | tar -C /tmp -zxf - && \
    cd /tmp/mvapich2-2.1 && ./configure --prefix=/opt/mvapich2-2.1 && make install && \
    rm -rf /tmp/mvapich2-2.1
RUN pip install schwimmbad
RUN MPICC=/opt/mvapich2-2.1/bin CFLAGS='-I /opt/mvapich2-2.1/include -L /opt/mvapich2-2.1/lib -lmpi' pip install --no-cache-dir mpi4py

# Install and set up cvmfs using static mounts
RUN yum install -y https://ecsft.cern.ch/dist/cvmfs/cvmfs-release/cvmfs-release-latest.noarch.rpm && yum install -y cvmfs cvmfs-config-default
ADD docker/etc/cvmfs/default.local /etc/cvmfs/default.local
ADD docker/etc/cvmfs/60-osg.conf /etc/cvmfs/60-osg.conf
ADD docker/etc/cvmfs/config-osg.opensciencegrid.org.conf /etc/cvmfs/config-osg.opensciencegrid.org.conf
RUN mkdir -p /cvmfs/config-osg.opensciencegrid.org /cvmfs/oasis.opensciencegrid.org /cvmfs/gwosc.osgstorage.org
RUN echo "config-osg.opensciencegrid.org /cvmfs/config-osg.opensciencegrid.org cvmfs ro,noauto 0 0" >> /etc/fstab && echo "oasis.opensciencegrid.org /cvmfs/oasis.opensciencegrid.org cvmfs ro,noauto 0 0" >> /etc/fstab && echo "gwosc.osgstorage.org /cvmfs/gwosc.osgstorage.org cvmfs ro,noauto 0 0" >> /etc/fstab

# Mount points for the XSEDE Comet machine
RUN mkdir -p /oasis /scratch /projects /usr/lib64/slurm /var/run/munge

# Create a PyCBC user
RUN groupadd -g 1000 pycbc && useradd -u 1000 -g 1000 -d /opt/pycbc -k /etc/skel -m -s /bin/bash pycbc

# Add the script that installs PyCBC inside the container
ADD docker/etc/docker-install.sh /etc/docker-install.sh

# When the container is started with 
#   docker run -it pycbc/pycbc-el7:latest 
# the default is to start a loging shell as the pycbc user. This can be overridden by e.g.
#   docker run -it pycbc/pycbc-el7:latest root /bin/bash -l 
# which will start a container with a root login shell
ENTRYPOINT ["/bin/su", "-l"]
CMD ["pycbc", "/bin/bash", "-l"]
