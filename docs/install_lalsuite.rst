.. _lalsuite_install:

##############################################
Installing lalsuite for PyCBC
##############################################

The following page describes how to build lalsuite from source for use with PyCBC. 

.. note::

    These instructions assume that you already have the required dependencies for lalsuite installed on your system. You must have ``git``, ``FFTW``, ``GSL``, ``FrameL``, and ``MetaIO`` installed before continuing. 
    
    These packages are installed by default on a LIGO Data Grid cluster. If you are not on a LIGO Data Grid cluster, then you can download them from the `lscsoft repository <https://www.lsc-group.phys.uwm.edu/daswg/download/repositories.html>`_ or see the `lscsoft install instructions <https://www.lsc-group.phys.uwm.edu/daswg/docs/howto/lscsoft-install.html>`_ for instructions on building the dependencies from source. You can also contact your system administrator to determine if these packages are available using the ``module load`` command.

====================================================
Obtaining the source code and checking out a version
====================================================

Clone the lalsuite repository into the directory ``${VIRTUAL_ENV}/src/lalsuite`` by following the `instructions on the DASWG pages <https://www.lsc-group.phys.uwm.edu/daswg/docs/howto/advanced-lalsuite-git.html#clone>`_. 

.. note::

    Since building lalsuite is very disk intensive, you may want to store the lalsuite git repository on a local disk rather than an NSF-mounted directory. 

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
    ./configure --prefix=${VIRTUAL_ENV}/opt/lalsuite --enable-swig-python --disable-lalstochastic --disable-lalxml --disable-lalinference --disable-laldetchar --disable-lalapps --with-hdf5=no

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

.. code-block:: bash

    cd $VIRTUAL_ENV/src/lalsuite/lalapps
    LIBS=-lz ./configure --prefix=${VIRTUAL_ENV}/opt/lalsuite --enable-static-binaries --disable-lalinference --disable-lalburst --disable-lalpulsar --disable-lalstochastic
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


