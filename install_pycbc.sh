#!/bin/bash
# Copyright 2015 Amber Lenon, Steve Reyes, Duncan Brown.

# Exit if any command fails
set -e 

#Installing pyCBC

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
  exit 1
fi

if [[ -d $NAME ]] ; then
   echo "ERROR: the directory $NAME already exists."
   echo "If you want to use this path, remove the directory and try again."
   exit 1
fi

read -p "Enter the number of processors that you want to use for builds: " nproc

if [[ -z $nproc ]] ; then
  nproc=1
fi

if [[ $nproc -lt 1 ]] ; then
   echo "ERROR: invalid number of processors specified: $nproc"
   exit 1
fi

if [[ $nproc -gt 24 ]] ; then
   echo "ERROR: please do not use more than 24 CPUs for the build. You asked for $nproc"
   exit 1
fi

echo "Using $nproc processors for parallel build."
echo

#Create a Virtual Environment
echo "--- creating virtual environment --------------------------------"
unset PYTHONPATH
virtualenv $NAME

#Enter Virtual Environment
source $NAME/bin/activate
mkdir -p $VIRTUAL_ENV/src

#Install HDF5
echo "--- installing HDF5 libraries -----------------------------------"
cd $VIRTUAL_ENV/src
curl https://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.12/src/hdf5-1.8.12.tar.gz > hdf5-1.8.12.tar.gz
tar -zxvf hdf5-1.8.12.tar.gz
rm hdf5-1.8.12.tar.gz
cd hdf5-1.8.12
./configure --prefix=$VIRTUAL_ENV/opt/hdf5-1.8.12
make -j $nproc install
HDF5_DIR=${VIRTUAL_ENV}/opt/hdf5-1.8.12 pip install h5py

#Installing lalsuite into Virtual Environment
#Install unitest2, python-cjson, and numpy
echo "--- installing required packages --------------------------------"
pip install "numpy>=1.6.4" unittest2 python-cjson Cython

#Authenticate with LIGO Data Grid services, install M2Crypto
SWIG_FEATURES="-cpperraswarn -includeall -I/usr/include/openssl" pip install M2Crypto

echo
echo
echo "--- cloning lalsuite repository ---------------------------------"
echo
read -p "Enter your LIGO.ORG username in (e.g. albert.einstein): " directory

#Valid ECP cookie to clone
echo "Enter your LIGO.ORG password to get a cookie to clone lalsuite."
ecp-cookie-init LIGO.ORG https://versions.ligo.org/git $directory

#Tell git the location of the cookie
git config --global http.cookiefile /tmp/ecpcookie.u`id -u`

#Get Copy of LalSuite Repository
cd $VIRTUAL_ENV/src
git clone https://versions.ligo.org/git/lalsuite.git
 
#Obtaining source code and checking version
#Change to lalsuite directory
cd lalsuite

#Determine which version of the code you want to install
echo
echo
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
git checkout $lalbranch

#Building and Installing into your Virtual Environment
#Use the master configure script to build and install all the components
echo
echo "--- building lalsuite -------------------------------------------"
echo
./00boot
./configure --prefix=${VIRTUAL_ENV}/opt/lalsuite --enable-swig-python --disable-lalstochastic --disable-lalxml --disable-lalinference --disable-laldetchar
make -j $nproc
make install

#Add to virtualenv activate script
echo 'source ${VIRTUAL_ENV}/opt/lalsuite/etc/lalsuiterc' >> ${VIRTUAL_ENV}/bin/activate
source ${VIRTUAL_ENV}/opt/lalsuite/etc/lalsuiterc

#Check that lalsuite is installed
echo "LAL installed into $LAL_PREFIX"

#Installing pyCBC to Virtual Environment
echo
echo "--- installing pegasus and dqsegdb ------------------------------"
echo

#Install Pegasus WMS python libraries
pip install http://download.pegasus.isi.edu/pegasus/4.5.2/pegasus-python-source-4.5.2.tar.gz

#Install dqsegb from Duncan's repository
pip install git+https://github.com/duncan-brown/dqsegdb.git@pypi_release#egg=dqsegdb

#Released or Development
echo
echo "--- downloading PyCBC -------------------------------------------"
echo
echo

while true ; do

  echo "Would you like to install a released version of PyCBC or a development copy?"
  echo "Please choose either"
  echo "  1. Release version"
  echo "  2. Development version"
  echo
  read -rp "Enter either 1 or 2: " dev_or_rel

  if [[ $dev_or_rel -eq 1 ]] ; then
    #Installing a released version of pyCBC
    #Choose release version
    echo "Please enter the name of a valid release tag. These can be found at"
    echo "   https://github.com/ligo-cbc/pycbc/releases"
    read -rp "Enter the tag name (e.g. v1.1.5): " reltag

    #Install Version
    pip install git+https://github.com/ligo-cbc/pycbc@${reltag}#egg=pycbc --process-dependency-links

    # continue with install
    break

  elif [[ $dev_or_rel -eq 2 ]] ; then

    #Fork pyCBC to your Github account
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
  
    #Install PyCBC source code from GitHub URL
    pip install -e git+git@github.com:${github}/pycbc.git#egg=pycbc --process-dependency-links
  
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

if [[ $dev_or_rel -eq 2 ]] ; then
  #Building and Installing Documentation
  #Install Sphinx and the helper tools
  echo
  echo "--- downloading documentation tools -----------------------------"
  echo
  pip install "Sphinx>=1.3.1"
  pip install sphinxcontrib-programoutput
  pip install numpydoc
  
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
echo "--- setting up optimized libraries ------------------------------"
echo
echo
echo "To run MKL optimized code, you need to enter the path and architecture"
echo "for the Intel optimized toolkit on your cluster. For example:"
echo "on sugar, enter"
echo "     /opt/intel/bin/compilervars.sh intel64"
echo "on atlas, enter"
echo "     /opt/intel/2015/bin/compilervars.sh intel64"
echo "If you do not have these tools installed, just press return."
echo 
read -p "Enter path and architecture for Intel compilervars.sh or press return: " intel_path

#Add script that sets up the MKL environment to virtualenv activate script
if [[ ! -z "${intel_path}" ]] ; then 
  echo "source ${intel_path}" >> ${VIRTUAL_ENV}/bin/activate
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

#Leave the virtual environment and exit
deactivate
echo "PyCBC setup complete"
echo
exit 0
