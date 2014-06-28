.. _lalsuite_install:

##############################################
Installing lalsuite, glue and pylal for PyCBC
##############################################

The following page describes how to build lalsuite, glue and pylal from source for use with PyCBC. 

---------------------------------
Choose a lalsuite install method
---------------------------------

A version of lalsuite and glue (but not pylal) is installed system-wide on LDG and XSEDE clusters. During science runs this version is recommended for use, but between science runs, and for development, it is better to install lalsuite, glue and pylal from source. To use this system version and just install pylal follow :ref:`systeminstall`. To install lalsuite, glue and pylal from source follow :ref:`sourceinstall`.

.. _systeminstall:

===============================================
Using the system installed version of lalsuite
===============================================

The system installation includes all executables except the pylal programs in ``/usr/bin``. You will need to build pylal from source against the system installed libraries. 

.. note::

   On the TACC XSEDE cluster, you will need to run module load git to make sure that you are running the correct version of git to clone lalsuite.

.. note::

   On the TACC XSEDE cluster it is recommended to install code into the $WORK directory not $HOME because of space limitations.

To do this, clone the lalsuite repository by following the `instructions on the DASWG pages <https://www.lsc-group.phys.uwm.edu/daswg/docs/howto/advanced-lalsuite-git.html#clone>`_

If you want to install a specific version of pylal, for example the ER5 release (v0.5.0) you can run

.. code-block:: bash

    cd /path/to/your/lalsuite
    git checkout ${TAG_NAME} # For example git checkout pylal-v0.5.0-v1

if this step is not done then the latest version of pylal (master) will be used.

You can then build pylal with the commands

.. code-block:: bash

    cd pylal
    python setup.py install --prefix=${HOME}/local/pylal-v0.5.0-v1

This installs pylal in ``${HOME}/local/pylal-v0.5.0-v1`` but you can change this to another directory, if you prefer.

Set up your environment to use this pylal with

.. code-block:: bash

    source ${HOME}/local/pylal-v0.5.0-v1/etc/pylal-user-env.sh

and you are done! Ignore the next section.

.. _sourceinstall:

===============================================
Building and installing your own lalsuite
===============================================
 
.. note::

    On the TACC XSEDE cluster, you will need to run module load git to make sure that you are running the correct version of git to clone lalsuite.

.. note::

    On the TACC XSEDE cluster, you will need to run module load condor to have access to condor_compile

.. note::

    On the TACC XSEDE cluster it is recommended to install code into the $WORK directory not $HOME because of space limitations.

Clone the lalsuite repository by following the `instructions on the DASWG pages <https://www.lsc-group.phys.uwm.edu/daswg/docs/howto/advanced-lalsuite-git.html#clone>`_. Once you have the repository cloned, you will need to checkout the master branch by running

.. code-block:: bash

    cd /path/to/your/lalsuite
    git checkout master

.. note::

    If you want to build a specific version change master for whatever version you want. For example if installing the ER5 release this would be ``git checkout lalsuite-v6.22``

Once you are on the branch, you can build and install lalsuite in the normal way.

The attached :download:`example script <resources/build_new_lalsuite.sh>` can be used to build and install the code. This script should be run from the directory containing your lalsuite git directory. Run the script with ./build_new_lalsuite.sh BRANCH SOURCEDIR INSTALLDIR NUMCORES where:

* BRANCH is the name of the branch that you want to install (e.g. master).
* SOURCEDIR is the path to the source directory. This is the directory that contains your lalsuite folder.
* INSTALLDIR is the directory where you want the code installed. *This must be under your NFS-mounted home directory so it is accessible to the cluster nodes running your jobs.*
* NUMCORES is the number of cores for a parallel build (e.g. 8).

For example, the following will build lalsuite, glue and pylal and install it in /home/$USER/local/master/

.. code-block:: bash

    sh ./build_new_lalsuite.sh master /home/$USER/ /home/$USER/local/master 8

This script will create a file INSTALLDIR/etc/lscsoftrc that can be sourced to set up your environment to used the installed code. In the example above, you would do

.. code-block:: bash

    source /home/$USER/local/master/etc/lscsoftrc

to set up your environment to use the installed code.

Congratulations, you now have lalsuite, glue and python set up and ready to use!
