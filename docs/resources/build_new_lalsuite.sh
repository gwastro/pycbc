#!/bin/bash -v
## -------------------- Clone the LALSuite Git Repository -------------------- ## 

# LSCSOFT_SRCDIR=/${HOME}/src/lscsoft
# mkdir -p ${LSCSOFT_SRCDIR}
#
# git config --global user.name "Albert Einstein"
# git config --global user.email albert.einstein@ligo.org
# cd ${LSCSOFT_SRCDIR}
# git clone albert.einstein@ligo-vcs.phys.uwm.edu:/usr/local/git/lalsuite.git


## -------------------- Build a new branch of LALSuite -------------------- ##

# Command line arguements
# 1 - source branch name
# 2 - source directory
# 3 - install directory
# 4 - number of parallel tasks during install

export BRANCH_NAME=${1}
export LSCSOFT_SRCDIR=${2}
export INSTALL_DIR=${3}
MAKE_CORES=${4}

set -e

if [ ${#} -ne 4 ] ; then
    echo "Usage: ${0} branch subdirectory ncores"
    echo "   build branch and install lalsuite into ${INSTALL_DIR} using ncores"
    exit 1
fi


#read -p "Install new lalsuite into ${INSTALL_DIR} (y/n)? "
#if [ "$REPLY" != "y" ] ; then
#  echo "Exiting without installing."
#  exit 1
#fi

cd ${LSCSOFT_SRCDIR}/lalsuite
git checkout ${BRANCH_NAME}
git clean -dxf

mkdir -p ${INSTALL_DIR}

SOURCE_SCRIPTDIR=${INSTALL_DIR}/etc
rm -f ${INSTALL_DIR}/lscsoftrc

# Build and Install all the lal packages expect LALApps
for d in lal lalframe lalmetaio lalxml lalsimulation lalburst lalinspiral lalstochastic lalpulsar lalinference;
    do pushd $d;
    if [ -f config.status ] ; then
        make distclean || echo "Warning: could not make distclean. Carrying on anyway"
    fi
    ./00boot
    ./configure --prefix=${INSTALL_DIR} --enable-swig-python --enable-pthread-lock
     make -j ${MAKE_CORES} install;
    echo "source ${INSTALL_DIR}/etc/${d}-user-env.sh" >> ${SOURCE_SCRIPTDIR}/lscsoftrc;
    popd;
    source ${SOURCE_SCRIPTDIR}/lscsoftrc;
done

# Build and Install LALApps
pushd lalapps
if [ -f config.status ] ; then
    make distclean || echo "Warning: could not make distclean. Carrying on anyway"
fi
./00boot
./configure --prefix=${INSTALL_DIR}
make -j ${MAKE_CORES} install
popd
echo "source ${INSTALL_DIR}/etc/lalapps-user-env.sh" >> ${SOURCE_SCRIPTDIR}/lscsoftrc
source ${SOURCE_SCRIPTDIR}/lscsoftrc

# Build and Install pylal & glue
for d in glue pylal;
    do pushd $d;
    python setup.py clean --all
    python setup.py install --prefix=${INSTALL_DIR}
    echo "source ${INSTALL_DIR}/etc/${d}-user-env.sh" >> ${SOURCE_SCRIPTDIR}/lscsoftrc;
    popd;
    source ${SOURCE_SCRIPTDIR}/lscsoftrc;
done

cd ~

lal-version
lalapps_version
