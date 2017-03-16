.. _lalsuite_install:

##############################################
Installing lalsuite for PyCBC
##############################################

The following page describes how to build lalsuite from source for use with PyCBC. 

.. note::

    These instructions assume that you already have the required dependencies for lalsuite installed on your system. You must have ``git``, ``FFTW``, ``GSL``, ``FrameL``, and ``MetaIO`` installed before continuing. 
    
    These packages are installed by default on a LIGO Data Grid cluster. If you are not on a LIGO Data Grid cluster, then you can download them from the `lscsoft repository <https://wiki.ligo.org/DASWG/SoftwareDownloads>`_ or see the `lscsoft install instructions <https://www.lsc-group.phys.uwm.edu/daswg/docs/howto/lscsoft-install.html>`_ for instructions on building the dependencies from source. You can also contact your system administrator to determine if these packages are available using the ``module load`` command.

====================================================
Obtaining the source code and checking out a version
====================================================

In your virtual environment, make a directory for the lalsuite source. Enter
your virtual environment with the command

.. code-block:: bash

    source ~/src/pycbc/bin/activate

changing the path to the activate script appropriately.  Clone the lalsuite repository into the directory ``${VIRTUAL_ENV}/src/lalsuite`` by running the commands

.. code-block:: bash

    mkdir -p ${VIRTUAL_ENV}/src
    cd ${VIRTUAL_ENV}/src
    git clone https://github.com/lscsoft/lalsuite.git

Note that this checks out a read-only repository. If you want a git repository
that you can edit, you can either fork this repository to your own GitHub
account or, if you have ``LIGO.ORG`` credentials, you can follow the 
`instructions on the DASWG pages for cloning lalsuite <https://www.lsc-group.phys.uwm.edu/daswg/docs/howto/advanced-lalsuite-git.html#clone>`_.

.. note::

    Since building lalsuite is very disk intensive, you may want to store the lalsuite git repository on a local disk rather than an NSF-mounted directory. If this is the case, change the path in the ``mkdir`` and ``cd`` above to a directory on a non-NFS mounted filesystem. This ie not required, as lalsuite will build on an NFS disk, it is just slower.

Once you have the repository cloned, you should change your working directory to the top-level of the repository with 

.. code-block:: bash

    cd ${VIRTUAL_ENV}/src/lalsuite

Now determine which version of the code you want to install. To run the latest (possibly unstable) version of the code, use the ``master`` branch by running the command:

.. code-block:: bash

    git checkout master

If you want to build a specific release, replace ``master`` with a release tag, for example ``lalsuite-v6.30``. See the `list of lalsuite tags <https://ligo-vcs.phys.uwm.edu/cgit/lalsuite/refs/tags>`_ for available tags. You can also check out a branch from the `list of lalsuite branches <https://ligo-vcs.phys.uwm.edu/cgit/lalsuite/refs/heads>`_ in the same way; replace ``master`` with the branch name, e.g. ``lalsuite_o1_branch``.  Once you have checked out a tag or a branch, you can build and install lalsuite.


=====================================================
Building and installing into your virtual environment
=====================================================

.. note::

    The install instructions below install lalsuite into a directory called ``opt/lalsuite`` under your virtual environment. You can remove lalsuite by removing this directory with ``rm -rf ${VIRTUAL_ENV}/opt/lalsuite`` from inside your virtual environment. If you want to install multiple versions of lalsuite in your virtual environment, you can chose different directories by specifying a different ``--prefix`` to configure below.

First make sure you are in your virtual environment by activiating it. If you
are not already in your virtual environment run the command

.. code-block:: bash

    source ~/src/pycbc/bin/activate

changing the string ``~/src/pycbc`` appropriately to point to your
environment. 

From the top-level lalsuite directory, you can use the master configure script to build and install all of the components that you need with the commands 

.. code-block:: bash

    ./00boot 
    ./configure --prefix=${VIRTUAL_ENV}/opt/lalsuite --enable-swig-python --disable-lalstochastic --disable-lalxml --disable-lalinference --disable-laldetchar --disable-lalapps

Next make the software and install it. If you are on a multicore machine, you
can speed this up by running ``make -j N`` where ``N`` is the number of
processors you want to use for the build (e.g. 16).

.. code-block:: bash

    make
    make install

The install process creates a shell script called ``lalsuite-user-env.sh`` that sources all of the ``${VIRTUAL_ENV}/opt/lalsuite/etc/lal*-user-env.sh`` scripts that set up the environment for lalsuite. You can add this to your virtualenv ``activate`` script so that it gets set up when you enter your virtual environment. To do this, run the commands

.. code-block:: bash

    echo 'source ${VIRTUAL_ENV}/opt/lalsuite/etc/lalsuite-user-env.sh' >> ${VIRTUAL_ENV}/bin/activate
    deactivate

You can reenter your virtual environment with the usual command

.. code-block:: bash

    source ~/pycbc/src/bin/activate

changing ``~/pycbc/src`` as appropriate for your virtual environment path.

.. note::

    If you want to manage multiple versions of lalsuite, it is not reccommended to source the ``lalsuite-user-env.sh`` script from your activate script.  You should just source it when you enter your virtual environment with the command ``source ${VIRTUAL_ENV}/opt/lalsuite/etc/lalsuite-user-env.sh``

lalsuite is now installed in your virtual environment. You can check this with the command

.. code-block:: bash

    echo $LAL_PREFIX

which should return the path to the installation under your virtual environment.

If you are running a pipeline that uses the old LALApps programs ``lalapps_inspinj`` or ``lalapps_coh_PTF_inspiral`` then you can optionally build and install these by running the commands

.. note::

    On the AEI atlas cluster, the HDF5 libraries are installed in a non-stardard location so you will need to add ``-L/usr/lib/x86_64-linux-gnu/hdf5/serial`` to the ``LIBS`` variable before configuring lalapps.

.. code-block:: bash

    cd $VIRTUAL_ENV/src/lalsuite/lalapps
    LIBS="-lhdf5_hl -lhdf5 -ldl -lz" ./configure --prefix=${VIRTUAL_ENV}/opt/lalsuite --enable-static-binaries --disable-lalinference --disable-lalburst --disable-lalpulsar --disable-lalstochastic --disable-lalxml
    cd $VIRTUAL_ENV/src/lalsuite/lalapps/src/lalapps
    make
    cd $VIRTUAL_ENV/src/lalsuite/lalapps/src/inspiral
    make lalapps_inspinj
    cp lalapps_inspinj $VIRTUAL_ENV/bin
    cd $VIRTUAL_ENV/src/lalsuite/lalapps/src/ring
    make lalapps_coh_PTF_inspiral
    cp lalapps_coh_PTF_inspiral $VIRTUAL_ENV/bin

.. note::

    The LALApps build above builds static binaries, so you will need static libraries for fftw, glibc, etc. installed on your system to do this. These libraries are present by default in a LIGO Data Grid environment. If you do not wish to build static LALApps programs, the omit the ``--enable-static-binaries`` option to the configure script.


In addition to lalsuite, the generation of certain template waveforms (e.g.  the reduced order model implementations of SEOBNRv2 and SEOBNRv4) requires addition data files from the `lalsuite-extra repository <https://svn.ligo.caltech.edu/svn/lalsuite-extra/>`_. These data can either be obtained by downloading and installing lalsuite-extra into your virtual environment or using a copy of the data from the CERN virtual filesystem.

To install the data into your virtual environment, run the commands

.. code-block:: bash

    cd ${VIRTUAL_ENV}/src
    svn co https://svn.ligo.caltech.edu/svn/lalsuite-extra/
    cd lalsuite-extra
    ./00boot
    ./configure --prefix=${VIRTUAL_ENV}/opt/lalsuite-extra
    make install
    echo 'export LAL_DATA_PATH=${VIRTUAL_ENV}/opt/lalsuite-extra/share/lalsimulation' >> ${VIRTUAL_ENV}/bin/activate

Then deactivate and activate your virtual environment.

Alternatively, follow the `instructions for installing CVMFS for OSG
<https://twiki.grid.iu.edu/bin/view/Documentation/Release3/InstallCvmfs>`_ and
run the command

.. code-block:: bash

    echo 'export LAL_DATA_PATH=/cvmfs/oasis.opensciencegrid.org/ligo/sw/pycbc/lalsuite-extra/current/share/lalsimulation' >> $VIRTUAL_ENV/bin/activate

to add the appropriate path to your virtual environment's ``activate`` script.
Then deactivate and activate your virtual environment.


