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
RUN mkdir -p /opt/intel/composer_xe_2015.0.090/mkl/lib/intel64
RUN curl https://git.ligo.org/ligo-cbc/pycbc-software/raw/efd37637fbb568936dfb92bc7aa8a77359c9aa36/x86_64/composer_xe_2015.0.090/composer_xe_2015.0.090.tar.gz | tar -C /opt/intel/composer_xe_2015.0.090/mkl/lib/intel64 -zxvf -
RUN curl https://software.intel.com/en-us/license/intel-simplified-software-license > /opt/intel/composer_xe_2015.0.090/mkl/lib/intel64/intel-simplified-software-license.html
RUN chmod go+rx /opt/intel/composer_xe_2015.0.090/mkl/lib/intel64/*.so
RUN chmod go+r /opt/intel/composer_xe_2015.0.090/mkl/lib/intel64/intel-simplified-software-license.html

# Install Python development tools
RUN yum install -y python-devel which

# Upgrade to latest python installers
RUN pip install --upgrade pip setuptools

# Install PyCBC dependencies
RUN pip install ipython jupyter

# set up environment
ADD etc/profile.d/pycbc.sh /etc/profile.d/pycbc.sh
ADD etc/profile.d/pycbc.csh /etc/profile.d/pycbc.csh

# add singularity profiles
COPY .singularity.d /.singularity.d
RUN cd / && \
    ln -s .singularity.d/actions/exec .exec && \
    ln -s .singularity.d/actions/run .run && \
    ln -s .singularity.d/actions/test .shell && \
    ln -s .singularity.d/runscript singularity

# Install MPI software needed for pycbc_inference
RUN yum install -y libibverbs libibverbs-devel libibmad libibmad-devel libibumad libibumad-devel librdmacm librdmacm-devel libmlx5 libmlx4
RUN curl http://mvapich.cse.ohio-state.edu/download/mvapich/mv2/mvapich2-2.1.tar.gz | tar -C /tmp -zxvf - && \
    cd /tmp/mvapich2-2.1 && ./configure --prefix=/opt/mvapich2-2.1 && make -j 8 install && \
    rm -rf /tmp/mvapich2-2.1
RUN pip install schwimmbad
RUN MPICC=/opt/mvapich2-2.1/bin CFLAGS='-I /opt/mvapich2-2.1/include -L /opt/mvapich2-2.1/lib -lmpi' pip install --no-cache-dir mpi4py

# Mount points for the XSEDE Comet machine
RUN mkdir -p /oasis /scratch /projects /usr/lib64/slurm /var/run/munge

# Create a PyCBC user
RUN groupadd -g 1000 pycbc && useradd -u 1000 -g 1000 -d /opt/pycbc -k /etc/skel -m -s /bin/bash pycbc

# When the container is started with 
#   docker run -it pycbc/pycbc-el7:latest 
# the default is to start a loging shell as the pycbc user. This can be overridden by e.g.
#   docker run -it pycbc/pycbc-el7:latest root /bin/bash -l 
# which will start a container with a root login shell
ENTRYPOINT ["/bin/su", "-l"]
CMD ["pycbc", "/bin/bash", "-l"]
