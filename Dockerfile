FROM ligo/lalsuite-dev:el7

USER root

RUN curl http://download.pegasus.isi.edu/wms/download/rhel/7/pegasus.repo > /etc/yum.repos.d/pegasus.repo
RUN yum clean all
RUN yum makecache
RUN yum -q -y update
RUN yum -q -y install lscsoft-backports-config
RUN yum -q -y install lscsoft-epel-config
RUN yum -q -y install lscsoft-ius-config
RUN yum clean all
RUN yum makecache
RUN yum -q -y install git2u-all lscsoft-all
RUN yum -q -y install zlib-devel libpng-devel libjpeg-devel libsqlite3-dev sqlite-devel db4-devel
RUN yum -q -y install tkinter libpng-devel lynx telnet
RUN yum -q -y install compat-glibc compat-glibc-headers
RUN yum -q -y install gd-devel audit-libs-devel libcap-devel nss-devel
RUN yum -q -y install xmlto asciidoc hmaccalc newt-devel 'perl(ExtUtils::Embed)' pesign elfutils-devel binutils-devel numactl-devel pciutils-devel
RUN yum -q -y install dejagnu sharutils gcc-gnat libgnat dblatex gmp-devel mpfr-devel libmpc-devel
RUN yum -q -y install libuuid-devel netpbm-progs nasm
RUN yum -q -y install libstdc++-static
RUN yum -q -y install gettext-devel avahi-devel dyninst-devel crash-devel latex2html emacs libvirt-devel
RUN yum -q -y install xmlto-tex patch
RUN yum -q -y install ant asciidoc xsltproc fop docbook-style-xsl.noarch
RUN yum -q -y install vim-enhanced
RUN yum -q -y install openssh-server
RUN yum -q -y install globus-gsi-cert-utils-progs gsi-openssh-clients osg-ca-certs ligo-proxy-utils
RUN yum -y install wget
RUN wget http://htcondor.org/yum/RPM-GPG-KEY-HTCondor
RUN rpm --import RPM-GPG-KEY-HTCondor
RUN wget http://htcondor.org/yum/repo.d/htcondor-stable-rhel7.repo
RUN mv htcondor-stable-rhel7.repo /etc/yum.repos.d/htcondor-stable-rhel7.repo
RUN yum -q -y install condor condor-classads condor-python condor-procd condor-external-libs
RUN yum -q -y install ecp-cookie-init
RUN yum -q -y install man-db

# set up cvmfs
RUN yum -q -y install osg-oasis
RUN echo "/cvmfs /etc/auto.cvmfs" > /etc/auto.master
ADD tools/cvmfs.default.local /etc/cvmfs/default.local

# uninstall lalsuite
RUN yum -q -y remove "*lal*"

# enable ssh
EXPOSE 22
ADD tools/pycbc-sshd /usr/bin/pycbc-sshd
RUN chmod +x /usr/bin/pycbc-sshd
RUN mkdir -p /var/run/sshd

# create a regular user account and switch to it
RUN useradd -ms /bin/bash pycbc
USER pycbc
WORKDIR /home/pycbc
RUN cp -R /etc/skel/.??* ~

RUN pip install virtualenv
RUN virtualenv pycbc-software ; \
      source ~/pycbc-software/bin/activate ; \
      pip install --upgrade pip ; \
      pip install six packaging appdirs ; \
      pip install --upgrade setuptools ; \
      pip install "numpy>=1.6.4" "h5py>=2.5" unittest2 python-cjson Cython decorator ; \
      pip install "scipy>=0.13.0" ; \
      SWIG_FEATURES="-cpperraswarn -includeall -I/usr/include/openssl" pip install M2Crypto ; \
      deactivate

RUN source ~/pycbc-software/bin/activate ; \
      mkdir -p ~/src ; \
      cd ~/src ; \
      git clone https://github.com/lscsoft/lalsuite.git ; \
      cd ~/src/lalsuite ; \
      git checkout 539c8700af92eb6dd00e0e91b9dbaf5bae51f004 ; \
      ./00boot ; \
      ./configure --prefix=${VIRTUAL_ENV}/opt/lalsuite --enable-swig-python \
        --disable-lalstochastic --disable-lalxml --disable-lalinference \
        --disable-laldetchar --disable-lalapps 2>&1 | grep -v checking ; \
      make install ; \
      echo 'source ${VIRTUAL_ENV}/opt/lalsuite/etc/lalsuite-user-env.sh' >> ${VIRTUAL_ENV}/bin/activate ; \
      deactivate

RUN source ~/pycbc-software/bin/activate ; \
      cd ~/src/lalsuite/lalapps ; \
      LIBS="-lhdf5_hl -lhdf5 -ldl -lz" ./configure --prefix=${VIRTUAL_ENV}/opt/lalsuite \
         --enable-static-binaries --disable-lalinference \
         --disable-lalburst --disable-lalpulsar \
         --disable-lalstochastic ; \
      cd ~/src/lalsuite/lalapps/src/lalapps ; \
      make ; \
      cd ~/src/lalsuite/lalapps/src/inspiral ;\
      make lalapps_inspinj ; \
      cp lalapps_inspinj $VIRTUAL_ENV/bin ; \
      cd ~/src/lalsuite/lalapps/src/ring ; \
      make lalapps_coh_PTF_inspiral ; \
      cp lalapps_coh_PTF_inspiral $VIRTUAL_ENV/bin ;\
      deactivate

RUN source ~/pycbc-software/bin/activate ; \
        pip install http://download.pegasus.isi.edu/pegasus/4.7.4/pegasus-python-source-4.7.4.tar.gz ; \
        pip install dqsegdb ; \
        pip install 'matplotlib==1.5.3' ; \
        pip install "Sphinx>=1.4.2" numpydoc sphinx-rtd-theme ; \
        pip install "git+https://github.com/ligo-cbc/sphinxcontrib-programoutput.git#egg=sphinxcontrib-programoutput" ; \
        pip install ipython jupyter ; \
        deactivate

RUN echo 'source ${HOME}/pycbc-software/bin/activate' >> ~/.bash_profile
