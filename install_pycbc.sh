#!/bin/bash
# Copyright 2015 Amber Lenon, Steve Reyes, Duncan Brown.

# Exit if any command fails
set -e 

# this makes sure everything is logged to a file
LOGPATH=${PWD}/install_pycbc_`date +%Y%m%d%H%M%S`.log
if [ "$1" != "noscript" ] ; then
    # just in case the user is calling us via bash or sh
    chmod +x $0
    exec script -q -c "$0 noscript" ${LOGPATH}
    exit 1;
fi


# Set up Python
while true; do
  echo "Do you want to use the OSG build of Python 2.7?"
  echo
  echo "This is used for building bundled executables for use on OSG."
  echo "If you are not sure, say no to use the standard system Python."
  echo
  read -p "Enter yes or no: " OSG_PYTHON
  echo

  if [[ $OSG_PYTHON == "yes" ]] ; then
    if [ ! -e /cvmfs/oasis.opensciencegrid.org/osg/modules/lmod/current/init/bash ] ; then
      echo "Error: Could not find sctipt to set up OSG modules."
      echo "Check that /cvmfs/oasis.opensciencegrid.org is mounted on this machine."
      exit 1
    fi
    unset PYTHONPATH
    unset LD_LIBRARY_PATH
    set +e
    source /cvmfs/oasis.opensciencegrid.org/osg/modules/lmod/current/init/bash
    set -e
    module load python/2.7
    break
  elif [[ $OSG_PYTHON == "no" ]] ; then
    break
  else
    echo "Please enter yes or no."
  fi
done

py_path=`which python`
py_ver=`python --version 2>&1`
echo "Using Python from ${py_path} which is ${py_ver}"
echo

#Check pip and virtualenv versions
while true ; do
  pip_ver=`pip --version`
  echo "You are using pip version: ${pip_ver}"
  read -rp "Is your version of pip greater than or equal to 7.1.0? (Enter yes or no) " pip_version
  
  if [[ $pip_version == "yes" ]] ; then
    break
  
  elif [[ $pip_version == "no" ]] ; then
    echo
    echo "You must have at least version 7.1.0 of pip to install virtualenv."
    echo "To set up pip follow the instructions at:"
    echo "http://ligo-cbc.github.io/pycbc/latest/html/install_virtualenv.html"
    echo
    exit 1
  
  else
    echo "Please enter yes or no"
  fi
done

pip install virtualenv --upgrade --user
export PATH=${HOME}/.local/bin:${PATH}

#Installing pyCBC
while true ; do
  echo
  echo "Are you installing a virtualenv to build bundled exectuables for distribution?"
  read -p "Enter yes or no (if you are not sure, say no): " IS_BUNDLE_ENV

  if [[ $IS_BUNDLE_ENV == "yes" ]] ; then
    UNIQUE_ID=pycbc-`uuidgen`
    NAME=${TMPDIR}/${UNIQUE_ID}
    break
  elif [[ $IS_BUNDLE_ENV == "no" ]] ; then
    while true; do
    # Ask the user where they want pycbc installed
    read -p "Enter the location where you want the virtual env created: " NAME

    if [[ $NAME == ~* ]] ; then
      if [[ ! "$NAME" =~ "/" ]] ; then
        NAME=${HOME}
      else
        # chomp the ~
        NAME=${NAME/##~}
        # chomp anything else before the first slash to catch e.g. ~alenon/
        NAME=${NAME#*\/}
        # expand to the user's home directory
        NAME=${HOME}/${NAME}
      fi
    fi

    if [[ -z $NAME ]] ; then
      echo "ERROR: you must specify a path for your virtual environment."
      continue
    fi

    if [[ -d $NAME ]] ; then
       echo "ERROR: the directory $NAME already exists."
       echo "If you want to use this path, remove the directory and try again."
    else
      break
    fi
    done
    break
  else
    echo "Please enter yes or no."
  fi
done

#Virualenv check
echo
echo "You are installing PyCBC in $NAME."
echo


#Number of Processors for the installation
read -p "Enter the number of processors that you want to use for builds: " nproc

if [[ -z $nproc ]] ; then
  nproc=1
fi

if [[ $nproc -lt 1 ]] ; then
   echo "ERROR: invalid number of processors specified: $nproc"
   continue
fi

if [[ $nproc -gt 24 ]] ; then
   echo "ERROR: please do not use more than 24 CPUs for the build. You asked for $nproc"
   continue
fi

echo "Using $nproc processors for parallel build."
echo


while true ; do
#Inputs for code

#LIGO.ORG username
read -p "Enter your LIGO.ORG username in (e.g. albert.einstein): " directory

#Lalsuite
echo "--- select lalsuite branch or tag -------------------------------"
echo
echo "Enter the name of the lalsuite branch or tag that you want to use."
echo 
echo "A list of branches can be found at:"
echo "   https://versions.ligo.org/cgit/lalsuite/refs/heads"
echo "A list of tags can be found at:"
echo "   https://versions.ligo.org/cgit/lalsuite/refs/tags"
echo
read -rp "Enter a valid tag or branch name (e.g. master or lalsuite_o1_branch): " lalbranch

#github
echo "Would you like to install a released version of PyCBC or a development copy?"
  echo "Please choose either"
  echo "  1. Release version"
  echo "  2. Development version"
  echo
  read -rp "Enter either 1 or 2: " dev_or_rel

  #released version
  if [[ $dev_or_rel -eq 1 ]] ; then
    #Installing a released version of pyCBC
    #Choose release version
    echo "Please enter the name of a valid release tag. These can be found at"
    echo "   https://github.com/ligo-cbc/pycbc/releases"
    read -rp "Enter the tag name (e.g. v1.1.5): " reltag

  #development
  elif [[ $dev_or_rel -eq 2 ]] ; then

    #Fork pyCBC to your Github account
    echo
    echo To install a development version, you will need a GitHub account.
    echo If you do not have a github account you will need to set one up.
    echo
    echo Once you have set up your GitHub account, follow the instruction at
    echo
    echo     https://help.github.com/articles/fork-a-repo/
    echo
    echo to fork the respository
    echo
    echo     https://github.com/ligo-cbc/pycbc
    echo
    echo into your own GitHub account. You will then need to enable ssh 
    echo keys on your GitHub account. You can follow these instructions:
    echo
    echo     https://help.github.com/articles/generating-ssh-keys/
    echo 
    read -rsp $'Once you have completed these steps, hit [Enter] to continue...\n' -n1 key
    echo
    #Create an ssh agent to connect to GitHub
    echo "Do you already have an ssh agent running with the key connected to GitHub?"
    while true ; do
      read -p "Enter yes or no (if you are not sure, enter no): " ssh_key
      if [[ $ssh_key == "yes" ]] ; then
        created_socket=""
        echo "Using $SSH_AUTH_SOCK"
        break
      elif [[ $ssh_key == "no" ]] ; then
        created_socket="yes"
        echo "Creating ssh agent to connect to GitHub:"
        eval `ssh-agent`
        echo "Please enter your ssh key passphrase."
        ssh-add
        break
      else
        echo "ERROR: please enter yes or no."
      fi
    done

    #Input Username
    read -rp "GitHub Username: " github

  else
    echo "You must enter 1 or 2."
 fi

#MKL Optimized Code
echo "To run MKL optimized code, you need to enter the path and architecture"
echo "for the Intel optimized toolkit on your cluster. For example:"
echo "on sugar, enter"
echo "     /opt/intel/bin/compilervars.sh intel64"
echo "on atlas, enter"
echo "     /opt/intel/2015/bin/compilervars.sh intel64"
echo "on ldas-grid, enter"
echo "     /opt/intel/composerxe/bin/compilervars.sh intel64"
echo "If you do not have these tools installed, just press return."
echo 
read -p "Enter path and architecture for Intel compilervars.sh or press return: " intel_path

#ROM Data Path
LAL_DATA_PATH=""
while true; do
echo
read -rp "Is the LAL Reduce Order Model (ROM) data installed on your cluster (Enter yes or no, if unsure type no)? " install_rom

if [[ ${install_rom} == "yes" ]] ; then
 echo
 read -rp "Please enter the path to the ROM data: " rom_path
     if [[ ${rom_path} == ~* ]] ; then
      if [[ ! "${rom_path}" =~ "/" ]] ; then
        rom_path=${HOME}
      else
        # chomp the ~
        rom_path=${rom_path/##~}
        # chomp anything else before the first slash to catch e.g. ~alenon/
        rom_path=${rom_path#*\/}
        # expand to the user's home directory
        rom_path=${HOME}/${rom_path}
      fi
    fi
 LAL_DATA_PATH="${rom_path}"
 break
#exists, and print message if not (and go back to top of loop).

elif [[ ${install_rom} == "no" ]] ; then
 echo
 read -rp "Do you want to download the ROM data now (Enter yes or no): " rom_download

  if [[ ${rom_download} == "no" ]] ; then
    echo "Please determine the location of the ROM data, or chose yes to download it. Exiting install."
    exit 1

  elif [[ ${rom_download} == "yes" ]] ; then
    echo
    read -rp "Enter the path to store the ROM data: " rom_store
    if [[ ${rom_store} == ~* ]] ; then
      if [[ ! "${rom_store}" =~ "/" ]] ; then
        rom_store=${HOME}
        break
      else
        # chomp the ~
        rom_store=${rom_store/##~}
        # chomp anything else before the first slash to catch e.g. ~alenon/
        rom_store=${rom_store#*\/}
        # expand to the user's home directory
        rom_store=${HOME}/${rom_store}
    
       if [ ! -d ${rom_store} ] ; then
        mkdir -p ${rom_store}
       fi
    
    LAL_DATA_PATH="${rom_store}"
    echo 
    break
     fi
   fi

  else
   continue

  fi

else
 continue

fi
done

echo "------------PLease check the inputs carefully-----------------"
echo
echo "LIGO.ORG username: " $directory
echo "Lalsuite Branch or Tag: " $lalbranch
echo "Development or released: " $dev_or_rel
echo "   (1 is release, 2 is development)"

if [[ $dev_or_rel -eq 1 ]] ; then
echo "Released Tag: " $reltag
fi

if [[ $dev_or_rel -eq 2 ]] ; then
echo "Github Username: " $github
fi

echo "Path and architecture for Intel compilervars.sh:" $intel_path
if [[ ${install_rom} == "yes" ]] ; then
echo "ROM data is located in ${LAL_DATA_PATH}"
fi

if [[ ${install_rom} == "no" ]] ; then
echo "ROM data will be install in ${LAL_DATA_PATH}"
fi
echo
echo
read -rp "Are these correct? (Enter yes or no) " questions 

 if [[ $questions == "yes" ]] ; then
 break

 elif [[ $questions == "no" ]] ; then
 continue

 else
  echo "You must enter yes or no or ctrl-c to exit."
  continue

 fi
done

while true ; do

while true ; do
echo "What would you like to do with your pip cache?"
echo "1. Use the existing pip cache.(Fastest)"
echo "2. Ignore the pip cache."
echo "3. Remove the pip cache. (Safest and slowest)"
echo
read -rp "Enter 1, 2 or 3: " pip_cache


if [[ $pip_cache != 1 ]] && [[ $pip_cache != 2 ]] && [[ $pip_cache != 3 ]] ; then
 echo "You must enter 1, 2, or 3." 
 echo
 echo
 continue

else
 break

fi
done

#Pip Cache Check
read -rp "You entered [ $pip_cache ]. Are you sure? (Enter yes or no) " check

if [[ $check == "yes" ]] ; then 

 if [[ $pip_cache -eq 1 ]] ; then 
 cache=""
 break
 fi

 if [[ $pip_cache -eq 2 ]] ; then
 cache="--no-cache-dir"
 break
 fi

 if [[ $pip_cache -eq 3 ]] ; then
 rm -rfv ${HOME}/.cache/pip
 cache=""
 break
 fi
 
elif [[ $check == "no" ]] ; then
continue

else
echo "You must enter yes or no."
continue

fi
done

#Valid ECP cookie to clone
echo "Enter your LIGO.ORG password to get a cookie to clone lalsuite."
ecp-cookie-init LIGO.ORG https://versions.ligo.org/git $directory
echo
echo
#Create a Virtual Environment
echo "--- creating virtual environment --------------------------------"
if [[ $OSG_PYTHON == "no" ]] ; then
  unset PYTHONPATH
  unset LD_LIBRARY_PATH
fi
virtualenv $NAME

#Enter Virtual Environment
source $NAME/bin/activate
mkdir -p $VIRTUAL_ENV/src

echo
echo "--- setting up optimized libraries ------------------------------"
echo
echo

#Add script that sets up the MKL environment to virtualenv activate script
if [[ ! -z "${intel_path}" ]] ; then 
  echo "source ${intel_path}" >> ${VIRTUAL_ENV}/bin/activate
  set +e
  source ${intel_path}
  set -e
fi

#Installing lalsuite into Virtual Environment
#Install unitest2, python-cjson, and numpy
echo "--- installing required packages --------------------------------"
if [[ $IS_BUNDLE_ENV == "yes" ]] ; then
  pip $cache install "numpy==1.9.3" unittest2 "python-cjson==1.1.0" "Cython==0.23.2"
else
  pip $cache install "numpy>=1.6.4" unittest2 python-cjson Cython
fi
pip $cache install git+http://github.com/ligo-cbc/pyinstaller.git@pycbc_install#egg=pyinstaller

#Install HDF5
echo "--- installing HDF5 libraries -----------------------------------"
cd $VIRTUAL_ENV/src
pip install "nose>=1.0.0"
curl https://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.12/src/hdf5-1.8.12.tar.gz > hdf5-1.8.12.tar.gz
tar -zxvf hdf5-1.8.12.tar.gz
rm hdf5-1.8.12.tar.gz
cd hdf5-1.8.12
./configure --prefix=$VIRTUAL_ENV
make -j $nproc install
if [[ $IS_BUNDLE_ENV == "yes" ]] ; then
  HDF5_DIR=${VIRTUAL_ENV} pip $cache install "h5py==2.5.0"
else
  HDF5_DIR=${VIRTUAL_ENV} pip $cache install h5py
fi

#Authenticate with LIGO Data Grid services, install M2Crypto
SWIG_FEATURES="-cpperraswarn -includeall -I/usr/include/openssl" pip $cache install M2Crypto
echo
echo
echo "--- cloning lalsuite repository -----------------------------------"
echo
#Valid ECP cookie to clone

#Tell git the location of the cookie
git config --global http.cookiefile /tmp/ecpcookie.u`id -u`

#Get Copy of LalSuite Repository
LALSUITE_BUILD_DIR=`mktemp --tmpdir -d -t lalsuite-XXXXXXXXXX`
ln -sf ${LALSUITE_BUILD_DIR} $VIRTUAL_ENV/src/lalsuite
cd $VIRTUAL_ENV/src/lalsuite
git clone https://versions.ligo.org/git/lalsuite.git
 
#Obtaining source code and checking version
#Change to lalsuite directory
cd lalsuite

#Determine which version of the code you want to install
echo
git checkout $lalbranch

#Building and Installing into your Virtual Environment
#Use the master configure script to build and install all the components
echo
echo "--- building lalsuite -------------------------------------------"
echo
./00boot
./configure --prefix=${VIRTUAL_ENV}/opt/lalsuite --enable-swig-python --disable-lalstochastic --disable-lalxml --disable-lalinference --disable-laldetchar --disable-lalapps --with-hdf5=no
make -j $nproc
make install

#Add to virtualenv activate script
echo 'source ${VIRTUAL_ENV}/opt/lalsuite/etc/lalsuiterc' >> ${VIRTUAL_ENV}/bin/activate
source ${VIRTUAL_ENV}/opt/lalsuite/etc/lalsuiterc

#Build a static lalapps_inspinj and install it
cd $VIRTUAL_ENV/src/lalsuite/lalsuite/lalapps
./configure --prefix=${VIRTUAL_ENV}/opt/lalsuite --enable-static-binaries --disable-mpi --disable-bambi --disable-cfitsio --disable-pss --disable-gds --disable-lalstochastic --disable-lalxml --disable-lalinference --disable-laldetchar --disable-lalpulsar
cd $VIRTUAL_ENV/src/lalsuite/lalsuite/lalapps/src/lalapps
make -j $nproc
cd $VIRTUAL_ENV/src/lalsuite/lalsuite/lalapps/src/inspiral
make lalapps_inspinj
cp lalapps_inspinj $VIRTUAL_ENV/bin

#Check that lalsuite is installed
echo "LAL installed into $LAL_PREFIX"

#Installing pyCBC to Virtual Environment
echo
echo "--- installing pegasus and dqsegdb ------------------------------"
echo

#Install Pegasus WMS python libraries
pip $cache install http://download.pegasus.isi.edu/pegasus/4.6.1/pegasus-python-source-4.6.1.tar.gz

#Install dqsegb from Duncan's repository
pip $cache install git+https://github.com/ligo-cbc/dqsegdb.git@pypi_release#egg=dqsegdb

#Install gracedb client tools
pip install ligo-gracedb

#Install pycbc and glue from non-cached versions to get the rpaths correct
pip $cache install pycbc-glue pycbc-pylal

#ROM Data Download

if [[ ${install_rom} == "no" ]] ; then
    pushd ${LAL_DATA_PATH}
    echo "--- Downloading ROM DATA ---------------------------------"
    svn co https://svn.ligo.caltech.edu/svn/lalsuite-extra/
    pushd lalsuite-extra
    ./00boot
    ./configure --prefix=${LAL_DATA_PATH}
    make install
    popd
    popd
    #echo "export LAL_DATA_PATH=${LAL_DATA_PATH}" >> ${VIRTUAL_ENV}/bin/activate
fi

#ROM Data Validation
rom_hash=('f82ddc5dc0b6fdc75122e767bd5e78c8' '62afa5351d6b775ac33cb4d898f0016b' 'a6829fa05437cc0aad81e3f8dae839cc' '98ea14b01e729d15ff666caa25afaed6' 'b41f0f7fbaf8be1d1848de7ee702bc67' '20ee260c870109766a6f048e20c7e10f' '96c384617edd8375ceaa03f9b7456467' '67d4f206fe19104fbc98b923b37318bb' 'd0bf97b4e17b5c9a7cfd222aaaafd742' 'c2ea5d296fee01abe16c0dd9e5f71f04' '412953726ca4bc72a810b27b810831c7' '4d5378935a7fba5e96f671581bce99fb' '31f48cb651a60837a3e99ee050aa9bc2' '727d31f6dc678aba8539817c8d0ae930' 'd0e1601c7cf4bd727d03e6cf7d2f722b' 'e6c243f76609cada55612cfe53f82e41' '08186a21682d2e73cb00a3ef35aa5c9c' '1ef7953a977a1fb551f59585c5d63d7a' 'b5923860bf021e6a2a23d743e5724bee' '2947032d0ad7ffde9704e24bf9e676f5')

if [[ ${install_rom} == "yes" ]] ; then
md5sum_hash=$( echo -n 'test' | md5sum ${LAL_DATA_PATH}/SEOBNRv2ROM*.dat | cut -d' ' -f1)
fi

if [[ ${install_rom} == "no" ]] ; then
md5sum_hash=$( echo -n 'test' | md5sum ${LAL_DATA_PATH}/share/lalsimulation/SEOBNRv2ROM*.dat | cut -d' ' -f1)
fi

for j in "${rom_hash[@]}" ; do
  if [[ ${rom_hash[*]} == ${md5sum_hash[*]} ]] ; then
    echo "All files are in ${LAL_DATA_PATH}." 
    if [[ ${install_rom} == "yes" ]] ; then
        echo "export LAL_DATA_PATH=${LAL_DATA_PATH}" >>  ${VIRTUAL_ENV}/bin/activate
    elif [[ ${install_rom} == "no" ]] ; then
        echo "export LAL_DATA_PATH=${LAL_DATA_PATH}/share/lalsimulation" >>  ${VIRTUAL_ENV}/bin/activate
    fi
    break
  fi

  if [[ ${rom_hash[*]} != ${md5sum_hash[*]} ]] ; then
    echo "The files are not the same."
    echo
    exit 1
  fi
done

#Released or Development
echo
echo "--- downloading PyCBC -------------------------------------------"
echo
echo

while true ; do

  if [[ $dev_or_rel -eq 1 ]] ; then
    #Installing a released version of pyCBC
    curl https://raw.githubusercontent.com/ligo-cbc/pycbc/${reltag}/requirements.txt > ${VIRTUAL_ENV}/requirements.txt
    pip install -r ${VIRTUAL_ENV}/requirements.txt
    #Install Version
    pip $cache install -e git+https://github.com/ligo-cbc/pycbc@${reltag}#egg=pycbc --process-dependency-links
    rm -f ${VIRTUAL_ENV}/src/pip-delete-this-directory.txt

    # continue with install
    break

  elif [[ $dev_or_rel -eq 2 ]] ; then

    #Fork pyCBC to your Github account
    echo 
  
    #Install PyCBC source code from GitHub URL
    pip $cache install -e git+git@github.com:${github}/pycbc.git#egg=pycbc --process-dependency-links
  
    #Prevent Pip from removing source directory
    rm -f ${VIRTUAL_ENV}/src/pip-delete-this-directory.txt
  
    #Connect the user to the main PyCBC repo and bring in the objects
    cd ${VIRTUAL_ENV}/src/pycbc
    git remote add upstream git@github.com:ligo-cbc/pycbc.git
    git fetch upstream
    
    #Continue with install
    break
  else
    echo "You must enter 1 or 2."
  fi

done

#Build static sbank
echo
echo "--- building static version of lalapps sbank tools --------------"
echo
pushd $VIRTUAL_ENV/src/lalsuite/lalsuite/lalapps/src/inspiral
for prog in lalapps_cbc_sbank_choose_mchirp_boundaries lalapps_cbc_sbank_merge_sims lalapps_cbc_sbank_pipe lalapps_cbc_sbank_plot_sim lalapps_cbc_sbank lalapps_cbc_sbank_sim lalapps_cbc_sbank_splitter
do
  pyinstaller ${prog}.py                       \
    --hidden-import scipy.linalg.cython_blas   \
    --hidden-import scipy.linalg.cython_lapack \
    --hidden-import scipy.special._ufuncs_cxx  \
    --hidden-import scipy.integrate            \
    --strip                                    \
    --onefile
    cp dist/${prog} $VIRTUAL_ENV/bin
done


if [[ $dev_or_rel -eq 2 ]] ; then
  #Building and Installing Documentation
  #Install Sphinx and the helper tools
  echo
  echo "--- downloading documentation tools -----------------------------"
  echo
  pip $cache install "Sphinx>=1.3.1"
  pip $cache install sphinxcontrib-programoutput
  pip $cache install numpydoc
  
  #patch the bug in numpydoc for python 2.6
  cat <<EOF > ${VIRTUAL_ENV}/plot_directive.patch
--- lib/python2.6/site-packages/matplotlib/sphinxext/plot_directive.py 2015-09-26 13:48:46.000000000 -0400
+++ lib/python2.6/site-packages/matplotlib/sphinxext/plot_directive.py 2015-09-24 20:59:35.843029957 -0400
@@ -333,8 +333,13 @@
     """
     Remove the coding comment, which six.exec_ doesn't like.
     """
-    return re.sub(
+    try:
+      return re.sub(
         "^#\s*-\*-\s*coding:\s*.*-\*-$", "", text, flags=re.MULTILINE)
+    except TypeError:
+      return re.sub(
+        "^#\s*-\*-\s*coding:\s*.*-\*-$", "", text, re.MULTILINE)
+
 
 #------------------------------------------------------------------------------
 # Template
EOF
  
  set +e
  patch -f -p0 ${VIRTUAL_ENV}/lib/python2.6/site-packages/matplotlib/sphinxext/plot_directive.py < ${VIRTUAL_ENV}/plot_directive.patch
  if [[ $? -eq 0 ]] ; then
    echo "patched sphinxext/plot_directive.py successfully"
  else
    echo "patch to sphinxext/plot_directive.py failed"
    echo "plots in documentation build may not work"
  fi
  rm ${VIRTUAL_ENV}/plot_directive.patch
  set -e
fi
  

echo

if [[ ! -z "$create_agent" ]] ; then
  echo "Terminating ssh agent $SSH_AGENT_PID"
  kill -TERM $SSH_AGENT_PID
fi

echo
echo "=================================================================="
echo
echo "PyCBC has been installed in a virtual environment in the directory"
echo
echo "  ${VIRTUAL_ENV}"
echo
echo "To use this virtual environment run the command"
echo
echo "  source ${VIRTUAL_ENV}/bin/activate"
echo

if [[ $dev_or_rel -eq 2 ]] ; then
  echo "A clone of your PyCBC repository has been placed in the directory"
  echo
  echo "   ${VIRTUAL_ENV}/src/pycbc"
  echo
  echo "and connected to the ligo-cbc/pycbc repository as 'upstream.'"
  echo
  echo "You can use this repository to edit and make changes to the code."
  echo "To install your updated code in your virtual environment, first"
  echo "make sure you have run the activate script, then from the directory"
  echo
  echo "   ${VIRTUAL_ENV}/src/pycbc"
  echo
  echo "run the command"
  echo
  echo "   python setup.py install"
  echo
fi

if [[ $IS_BUNDLE_ENV == "yes" ]] ; then
  echo
  echo "=================================================================="
  echo
  echo "Making bundled executables."
  echo 
  #Make the dag that makes the bundles
  cd ${VIRTUAL_ENV}/src/pycbc/tools/static
  ./build_dag.sh

  #Run the dag that makes the bundles
  condor_submit_dag build_static.dag

  # Wait for the dag to finish
  condor_wait -status -echo build_static.dag.dagman.log

  #Copy the existing static files into the dist directory
  cp -v ${VIRTUAL_ENV}/bin/lalapps_inspinj dist/
  for prog in lalapps_cbc_sbank_choose_mchirp_boundaries lalapps_cbc_sbank_merge_sims lalapps_cbc_sbank_pipe lalapps_cbc_sbank_plot_sim lalapps_cbc_sbank lalapps_cbc_sbank_sim lalapps_cbc_sbank_splitter
    do cp -v ${VIRTUAL_ENV}/bin/$prog dist/
  done

  #Copy the completed build
  cp -v dist/ ${VIRTUAL_ENV}/../${UNIQUE_ID}-dist
  cd ${VIRTUAL_ENV}/..
  tar -zcvf ${UNIQUE_ID}.tar.gz ${VIRTUAL_ENV}

  deactivate
  echo 
  rm -rf ${VIRTUAL_ENV}

  echo
  echo "Complete."
  echo
  echo "Bundles are in ${UNIQUE_ID}-dist"
  echo "Environment saved in ${UNIQUE_ID}.tar.gz"
  echo 
  echo "=================================================================="
else
  #Leave the virtual environment and exit
  deactivate
  echo "PyCBC setup complete"
fi

echo

exit 0
