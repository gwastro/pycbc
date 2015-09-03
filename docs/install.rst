################
Installing PyCBC
################

There are three typical use cases for PyCBC:

1. Installing a release of PyCBC from GitHub for an end user to run the tools.
2. Installing an editable version from GitHub for development.
3. Production LIGO astrophysical searches.

This page documents the first two use cases. For production analysis, users must obtain the pre-built binaries from the PyCBC server. 

.. note::

    PyCBC uses the `fork and pull <https://help.github.com/articles/using-pull-requests/>`_ model for development. If you wish to develop PyCBC, then you will need an account on `GitHub <https://www.github.com>`_. Once you have set up your account you should follow the instructions to `fork a repository <https://help.github.com/articles/fork-a-repo/>`_ to fork the `ligo-cbc/pycbc <https://github.com/ligo-cbc/pycbc>`_ repository into your own account. From your own fork, you can follow the `GitHub flow model <https://help.github.com/articles/github-flow-in-the-browser/>`_ to develop and maintain the code. For each new feature or bug fix, you should `create a new branch <https://help.github.com/articles/creating-and-deleting-branches-within-your-repository/>`_ to develop the feature. You can then `create a pull request <https://help.github.com/articles/creating-a-pull-request/>`_ so that the PyCBC maintainers can review and merge your changes into the official repository.

These instructions walk you through the process of

    * Initial Setup to install PyCBC
        * `Setting up a virtual environment for installing PyCBC`_.
        * `Installing lalsuite into a virtual environment`_.
    * `Installing PyCBC in a virtual environment`_.
        * `Installing a released version of PyCBC`_.
        * `Installing source from GitHub for development`_.
    * Optional additional installation steps
        * `Building and Installing Documentation`_.
        * `Modifying pycbc-glue and pycbc-pylal`_.
        * `Use of Intel MKL Optimized FFT libraries`_.
        * `Graphics Processing Unit support with CUDA`_

=====================================================
Setting up a virtual environment for installing PyCBC
=====================================================

The recommended way of installing PyCBC is to use `pip <https://pip.pypa.io/en/stable/>`_ within a `Python Virtual Environment <https://virtualenv.pypa.io/en/latest/>`_. Virtualenv isolates PyCBC and its dependencies from the system environment and installing with pip ensures that PyCBC picks up the correct dependencies. The following instructions will create a working virtual environment into which you can install PyCBC. 

Make sure that you have at least version 13.1.1 of virtualenv by running 

.. code-block:: bash

    virtualenv --version
    
If this returns ``virtualenv: command not found`` (as is the case with LIGO Data Grid Scientific Linux 6 systems) or the command returns a lower version than 13.1.1, then follow the instructions for setting virtualenv at:

.. toctree::
    :maxdepth: 1

    install_virtualenv

Once you have virtualenv installed, unset your current ``PYTHONPATH`` so that you do not inherit any packages that may conflict with your installation. To do this, run the command:

.. code-block:: bash

    unset PYTHONPATH

By default, virtualenv will modify your shell prompt so that it prepends the name of the virtual environment. This can be useful to make sure that you are developing in the virtual environment, or if you have several virtual environments. However, if you do not want this, then set

.. code-block:: bash

    VIRTUAL_ENV_DISABLE_PROMPT=True
    
Before running the command to create the new virtual environment.

Next, you need to choose a directory name where you'd like to make your virtual environment, and then make it. In this example, we use ``${HOME}/pycbc-dev`` but this can be changed to any path that you you have write access to, except for ``${HOME}/.local`` as that will cause conflicts with pip.

Set up the new virtual environment with the command

.. code-block:: bash

    NAME=${HOME}/pycbc-dev
    virtualenv $NAME
    
To enter your virtual environment run the command

.. code-block:: bash
    
    source $NAME/bin/activate
    
You will now be in your virtual environment, and so you can install packages, etc, without conflicting with either the system build, or other builds that you may have sitting around. You may install other programs and libraries, such as lalsuite (:ref:`lalsuite_install`), into this virtual environment.

To leave this virtual environment type

.. code-block:: bash

    deactivate
    
which will return you to a regular shell.

==============================================
Installing lalsuite into a virtual environment
==============================================

Enter the virtual environment that you wish to use for PyCBC development by sourcing the activate script, for example

.. code-block:: bash

    source $NAME/bin/activate

First install unittest2, python-cjson, and numpy with the command:

.. code-block:: bash

    pip install "numpy>=1.6.4" unittest2 python-cjson Cython

To authenticate with LIGO Data Grid services, you need M2Crypto which you should install with 

.. code-block:: bash

    SWIG_FEATURES="-cpperraswarn -includeall -I/usr/include/openssl" pip install M2Crypto

Once you have these packages installed, you can now install lalsuite following the instructions at:

.. toctree::
    :maxdepth: 1

    install_lalsuite

=========================================
Installing PyCBC in a virtual environment
=========================================

Enter the virtual environment that you wish to use for PyCBC development by sourcing the activate script, for example

.. code-block:: bash

    source $NAME/bin/activate

Next install the Pegasus WMS python libraries needed to build the workflows with the command:

.. code-block:: bash

    pip install http://download.pegasus.isi.edu/pegasus/4.5.2/pegasus-python-source-4.5.2.tar.gz

You now need to decide whether you want to install a release of PyCBC or an editable version of the source code from a git repository for development. 

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Installing a released version of PyCBC
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

.. note::

    Make sure you have run the commands in the section :ref:`Installing lalsuite into a virtual environment` above to install unittest2 before installing PyCBC.

To install a release of the code, determine the tag of the release that you want to install from the `list of PyCBC tags <https://github.com/ligo-cbc/pycbc/tags>`_. This example installs the v1.1.0 release. If you want to install a different release, change the command below accordingly:

.. code-block:: bash

    pip install git+https://github.com/ligo-cbc/pycbc@v1.1.0#egg=pycbc --process-dependency-links

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Installing source from GitHub for development
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

.. note::

    Make sure you have run the command in the section :ref:`Installing lalsuite into a virtual environment` above to install unittest2 before installing PyCBC.

To install and editable version of PyCBC you need to have `forked PyCBC to your own account <https://help.github.com/articles/fork-a-repo/>`_ and know the URL of your fork. This can be obtained from the clone URL on your GitHub repository page. This example uses the URL ``git@github.com:your-username-here/pycbc.git`` which you should change to the URL of your fork of PyCBC on GitHub. 

You can also read the `pip git instructions <https://pip.pypa.io/en/latest/reference/pip_install.html#git>`_ for more details on how to install a branch or a specific tag.

Install the PyCBC source code from the GitHub URL using the command:

.. code-block:: bash

    pip install -e git+git@github.com:duncan-brown/pycbc.git#egg=pycbc --process-dependency-links

This will fetch the PyCBC source and will also install all the listed dependencies. The ``-e`` option to pip creates a directory called ``$NAME/src/pycbc`` with a git checkout which is fully editable. To prevent pip from removing this source directory run the command

.. code-block:: bash

    rm -f $NAME/src/pip-delete-this-directory.txt

You can then make changes to your PyCBC source code in the directory ``$NAME/src/pycbc``

You can also use the repository created by pip as your working repository, creating branches, commits, and `pull requests <https://help.github.com/articles/creating-a-pull-request/>`_ as you need to. To keep your repository in sync with the ligo-cbc/pycbc repository, you can read the GitHub instructions that explaion how to `sync a fork of a repository to keep it up-to-date with the upstream repository. <https://help.github.com/articles/syncing-a-fork/>`_.

.. note::

    The version of PyCBC that is checked out will be on the master branch. To track a remote branch from your GitHub repository, run the command
    
    ``git branch --set-upstream branch_name origin/branch_name``

    ``git checkout branch_name``
        
    where branch_name is the name of the branch that you want to track.

To build and install any changes that you make to the source code in your virtual environment, run the command

.. code-block:: bash

    python setup.py install
    
from the PyCBC source directory in ``$NAME/src/pycbc``

=====================================
Building and Installing Documentation
=====================================

To build the documentation from your virtual environment, first make sure that you have `Sphinx <http://sphinx-doc.org/>`_ and the required helper tools installed with

.. code-block:: bash

    pip install "Sphinx>=1.3.1"
    pip install sphinxcontrib-programoutput
    pip install numpydoc
    
To generate the documentation, from the top level of the PyCBC source tree run

.. code-block:: bash

    python setup.py build_docs
    
This will build the documentation in the directory docs/_build/html which can be copied to a web-accessible directory. For example

.. code-block:: bash

    cp -a docs/_build/html/ ~/public_html/pycbc-docs
    
will copy the documentation to a directory called ``pycbc-docs`` under your public html pages.

To maintain the documentation under GitHub project pages, see

.. toctree::
    :maxdepth: 1

    build_gh_pages


====================================
Modifying pycbc-glue and pycbc-pylal
====================================

PyCBC depends on the packages pycbc-glue and pycbc-pylal which are forks of the lalsuite development of these packages. The correct versions are automatically installed by pip as part of the main PyCBC install. If you are developing code in these packages, then you can clone them from GitHib into your virtual environment's source directory and build and install them from there.

.. note::

    If you want to develop pycbc-glue and pycbc-pylal, you should follow the instructions to `fork a repository <https://help.github.com/articles/fork-a-repo/>`_ to fork the `ligo-cbc/pycbc-glue <https://github.com/ligo-cbc/pycbc-glue>`_ and `ligo-cbc/pycbc-pylal <https://github.com/ligo-cbc/pycbc-pylal>`_ repositories into your own account.

You can obtain these repositories in the standard way using git, replacing ``ligo-cbc`` with your GitHub user account name

.. code-block:: bash

    cd $NAME/src
    git clone git@github.com:ligo-cbc/pycbc-glue.git
    git clone git@github.com:ligo-cbc/pycbc-pylal.git

Once you have the source code cloned, you can run 

.. code-block:: bash

    python setup.py install

to install each of them into your virtual environment.

========================================
Use of Intel MKL Optimized FFT libraries
========================================

PyCBC has the ability to use optimized FFT libraries such as FFTW and MKL. If MKL is the correct library for your platform, you can add the script that sets up the MKL environment to you virtualenv ``activate`` script with the command

.. code-block:: bash

    echo 'source /opt/intel/bin/compilervars.sh intel64' >> $NAME/bin/activate

changing the path to the ``compilervars.sh`` script approriately for your cluster. 

==========================================
Graphics Processing Unit support with CUDA
==========================================

PyCBC has the ability to accelerate its processing using CUDA. To take advantage of this, follow the instructions linked below to install the CUDA dependencies.

.. toctree::
    :maxdepth: 1

    install_cuda
