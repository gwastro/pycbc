#!/bin/bash

NAME=${HOME}/bundle_env2.6

if [ ! -e ${HOME}/src/lalsuite ] && [ ! -e ${HOME}/lalsuite ]
then
	echo "Please checkout lalsuite under ${HOME}/src/lalsuite"
	echo "This can be done with"
	echo "mkdir -p ${HOME}/src"
	echo "cd ${HOME}/src"
	echo "git clone git://versions.ligo.org/lalsuite.git"
	echo "cd lalsuite"
	echo "git checkout o1_lalinference_20151210_SpinFix"
	echo "git cherry-pick d378c1f37597447ab53e1114daba2805477bfd0f"

	exit 1
fi

if [ `hostname` == 'sugar-dev2.phy.syr.edu' ]
then
	source /cvmfs/oasis.opensciencegrid.org/osg/modules/lmod/current/init/bash
	module load python/2.7
	NAME=${HOME}/bundle_env2.7
fi

cd ~
rm -rf local .local

VERSION=`python -c 'import sys; print "%d.%d" % (sys.version_info[0], sys.version_info[1])'`
mkdir -p ${HOME}/local/pip-7.1.0/lib/python${VERSION}/site-packages
export PYTHONPATH=${HOME}/local/pip-7.1.0/lib/python${VERSION}/site-packages
export PATH=${HOME}/local/pip-7.1.0/bin:${PATH}
easy_install --prefix=${HOME}/local/pip-7.1.0 https://pypi.python.org/packages/source/p/pip/pip-7.1.0.tar.gz#md5=d935ee9146074b1d3f26c5f0acfd120e
rm -rf .cache/pip/

pip install virtualenv --upgrade --user

virtualenv $NAME
source $NAME/bin/activate

mkdir -p $VIRTUAL_ENV/src
cd $VIRTUAL_ENV/src
curl https://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.12/src/hdf5-1.8.12.tar.gz > hdf5-1.8.12.tar.gz
tar -zxvf hdf5-1.8.12.tar.gz
cd hdf5-1.8.12
./configure --prefix=$VIRTUAL_ENV/opt/hdf5-1.8.12
make install

pip install "numpy==1.9.3" unittest2 "python-cjson==1.1.0" "Cython==0.23.2"
HDF5_DIR=$VIRTUAL_ENV/opt/hdf5-1.8.12 pip install "h5py==2.5.0"
SWIG_FEATURES="-cpperraswarn -includeall -I/usr/include/openssl" pip install M2Crypto

pushd ~/src/lalsuite
git pull
git clean -dxf
./00boot
./configure --prefix=$NAME/opt/lalsuite --enable-swig-python --disable-lalstochastic --disable-lalxml --disable-lalinference --disable-laldetchar --disable-lalapps
make -j
make install
echo 'source ${VIRTUAL_ENV}/opt/lalsuite/etc/lalsuiterc' >> $NAME/bin/activate
deactivate
source $NAME/bin/activate


cd ~
rm -f requirements.txt
wget https://raw.githubusercontent.com/ligo-cbc/pycbc/master/requirements.txt
sed -i 's+pegasus-wms==4.5.2+http://download.pegasus.isi.edu/pegasus/4.5.3/pegasus-python-source-4.5.3.tar.gz+' requirements.txt 
pip install -r requirements.txt 

cd ${NAME}
git clone git@github.com:ligo-cbc/pycbc.git
cd pycbc/
git remote add upstream git@github.com:ligo-cbc/pycbc.git
python setup.py install

echo
echo
echo " ============ DONE! ============ "
echo "Everything is set up, but you may want to"
echo "cd ${NAME}/src"
echo "rm -rf pycbc"
echo "checkout from your own github account"
echo "cd pycbc"
echo "git remote add upstream git@github.com:ligo-cbc/pycbc.git"


