.. _lalsuite_install:

##############################################
Installing lalsuite for PyCBC
##############################################

The following page describes how to build lalsuite from source for use with PyCBC. 

====================================================
Obtaining the source code and checking out a version
====================================================

Clone the lalsuite repository by following the `instructions on the DASWG pages <https://www.lsc-group.phys.uwm.edu/daswg/docs/howto/advanced-lalsuite-git.html#clone>`_. 

.. note::

    Since building lalsuite is very disk intensive, you may want to store the lalsuite git repository on a local disk rather than an NSF-mounted directory. 

Once you have the repository cloned, you should change your working directory to the top-level of the repository with 

.. code-block:: bash

    cd /path/to/your/lalsuite

Now determine which version of the code you want to install. To run the latest (possibly unstable) version of the code, use the ``master`` branch by running the command:

.. code-block:: bash

    git checkout master

If you want to build a specific release, replace ``master`` with a release tag, for example ``lalsuite-v6.30``. See the `list of lalsuite tags <https://ligo-vcs.phys.uwm.edu/cgit/lalsuite/refs/tags>`_ for available tags. You can also check out a branch from the `list of lalsuite branches <https://ligo-vcs.phys.uwm.edu/cgit/lalsuite/refs/heads>`_ in the same way; replace ``master`` with the branch name, e.g. ``lalsuite_o1_branch``.  Once you have checked out a tag or a branch, you can build and install lalsuite.


=====================================================
Building and installing into your virtual environment
=====================================================

Set the shell variable ``NAME`` to the path to your the virtual environment that you created for PyCBC and activate your environment, for example

.. code-block:: bash

    NAME=${HOME}/pycbc-dev
    source $NAME/bin/activate

From the top-level lalsuite directory, you can use the master configure script to build and install all of the components that you need with the commands 

.. code-block:: bash

    ./00boot 
    ./configure --prefix=$NAME --enable-swig-python --disable-lalstochastic --disable-lalxml --disable-lalinference --disable-laldetchar --disable-lalburst
    make
    make install

The install process creates a shell script called ``$NAME/etc/lal-user-env.sh`` that sets up the environment for lalsuite. You can add this to your virtualenv ``activate`` script so that it gets set up when you enter your virtual environment. To do this, run the commands

.. code-block:: bash

    echo 'source ${VIRTUAL_ENV}/etc/lal-user-env.sh' >> $NAME/bin/activate
    deactivate
    source $NAME/bin/activate

lalsuite is now installed in your virtual environment. You can check this with the command

.. code-block:: bash

    echo $LAL_PREFIX

which should return the path to the installation under your virtual environment.


