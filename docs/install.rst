.. _installing_pycbc:

################
Installing PyCBC
################

This document explains how to set up a virtual environment to install PyCBC
either for development or use in a production analysis with a release. The
code build will be a standard Python install which requires that the
installation directory containing the Python libraries is accessible at
runtime. Some executables also use weave for just-in-time compilation of code
at runtime. These executables require a gcc and Python build environment on
the execution machine.

If you wish to run PyCBC executables on a machine that does not have the
required environment, then you must use PyInstaller to build bundled versions
of the executables. Documentation on doing this is available on the page
:ref:`building_bundled_executables`.

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
        * `Modifying pycbc-glue`_.
        * `Use of Intel MKL Optimized FFT libraries`_.
        * `Graphics Processing Unit support with CUDA`_

=====================================================
Setting up a virtual environment for installing PyCBC
=====================================================

The recommended way of installing PyCBC is to use `pip <https://pip.pypa.io/en/stable/>`_ within a `Python Virtual Environment <https://virtualenv.pypa.io/en/latest/>`_. Virtualenv isolates PyCBC and its dependencies from the system environment and installing with pip ensures that PyCBC picks up the correct dependencies. The following instructions will create a working virtual environment into which you can install PyCBC. 

Make sure that you have at least version 13.1.1 of virtualenv. To do this, you can run the command below and check the return result. Run this command and note the result:

.. code-block:: bash

    virtualenv --version
    
If this returns ``virtualenv: command not found`` or the command returns a lower version than ``13.1.1`` (as is the case with Scientific Linux 6 systems) then follow the instructions for setting virtualenv at:

.. toctree::
    :maxdepth: 1

    install_virtualenv

Once you have a version virtualenv installed that is at least as new as ``13.1.1`` you can continue with these instructions.

Unset your current ``PYTHONPATH`` so that you do not inherit any packages that may conflict with your installation. To do this, run the command:

.. code-block:: bash

    unset PYTHONPATH

By default, virtualenv will modify your shell prompt so that it prepends the name of the virtual environment. This can be very useful to make sure that you are developing in the virtual environment, or if you have several virtual environments, so it is not recommended to disable this. However, if you do not want your prompt changed, then you can set

.. code-block:: bash

    VIRTUAL_ENV_DISABLE_PROMPT=True
    
Before running the command to create the new virtual environment.

Next, you need to choose a directory name where you'd like to make your virtual environment, and then make it. In this example, we use ``${HOME}/pycbc-dev`` but this can be changed to any path that you you have write access to, except for ``${HOME}/.local`` as that will cause conflicts with pip.

.. note::

    It is very important that your ``virtualenv`` version is at least 13.1.1 before continuing. Read the preceding instructions if you are not sure how to check this, or need to upgrade virtualenv.

You first set up a new virtual environment. A virtual environment is defined by a directory path that will hold the contents of the virtual environment. In this example, we chose ``~/src/pycbc`` although you can change this to whatever you like by changing the command below. Initialize the virtual environment by running the command:

.. code-block:: bash
    
    virtualenv ~/src/pycbc

You can create as many different virtual environments as you like, as long as
they all have different paths. 

.. note::

    Do not run ``virtualenv`` twice with the same path as it will overwrite the existing virtual environment with a new one, destroying the environment.

To enter your virtual environment run the command (replacing the string ``~/src/pycbc/`` if you chose a different path)

.. code-block:: bash
    
    source ~/src/pycbc/bin/activate

After running the ``activate`` script, you will now be in your virtual environment, and so you can install packages without conflicting with either the system build, or other builds that you may have sitting around. You may install other programs and libraries, such as lalsuite (:ref:`lalsuite_install`), into this virtual environment. The ``activate`` script also sets the environment variable ``${VIRTUAL_ENV}`` to the full path to your virtual environment.

.. note::

    Python implements a `per user site packages directory <https://www.python.org/dev/peps/pep-0370/>`_ to install packages in the user's home directory. By default this is ``~/.local`` but the location is controlled by the ``PYTHONUSERBASE`` environment variable. If you make use of the ``~/.local`` directory, you should add the line ``export PYTHONUSERBASE=${VIRTUAL_ENV}/.local`` to the end of your virtual environment's ``activate`` script to prevent conflicts. Similarly, pip caches data in the directory ``~/.caches/pip``. To prevent conflicts with this directory, you can add the line ``export XDG_CACHE_HOME=${VIRTUAL_ENV}/.cache`` to the virtual environment ``activate`` script so that pip uses a cache that is specific to the virtual environment.

To leave this virtual environment type

.. code-block:: bash

    deactivate
    
which will return you to a regular shell.


.. _installinglalsuite:

==============================================
Installing lalsuite into a virtual environment
==============================================

Enter the virtual environment that you wish to use for PyCBC development by sourcing the activate script, for example

.. code-block:: bash

    source ~/src/pycbc/bin/activate

.. note::

   CentOS 6 provides a buggy version of the HDF5 library, so you will need to install a newer version into your virtual environment. If you are using a CentOS 6 cluster, you must install HDF5. If you are using another cluster, then this step is optional.  

If you are running on a Scientific Linux 6 cluster, you need to install the HDF5 library. To do this, run the commands:

.. code-block:: bash

    mkdir -p $VIRTUAL_ENV/src
    cd $VIRTUAL_ENV/src
    pip install "nose>=1.0.0"
    curl https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.12/src/hdf5-1.8.12.tar.gz > hdf5-1.8.12.tar.gz
    tar -zxvf hdf5-1.8.12.tar.gz
    cd hdf5-1.8.12
    ./configure --prefix=$VIRTUAL_ENV/opt/hdf5-1.8.12
    make install
    HDF5_DIR=$VIRTUAL_ENV/opt/hdf5-1.8.12 pip install h5py


.. note::

    On a Scientific Linux 6 or Scientific Linux 7 system, if you do not upgrade the version of setuptools that is in your virtual environment, then the PyCBC install will fail to process the required dependencies.

Make sure your versions of ``pip`` and ``setuptools`` are up to date by running the command:

.. code-block:: bash

    pip install --upgrade pip
    pip install --upgrade setuptools

Install unittest2, python-cjson, and numpy with the command:

.. code-block:: bash

    pip install "numpy>=1.6.4" unittest2 python-cjson Cython decorator

To authenticate with LIGO Data Grid services, you need M2Crypto which you should install with 

.. code-block:: bash

    SWIG_FEATURES="-cpperraswarn -includeall -I/usr/include/openssl" pip install M2Crypto

On MacOS using homebrew to install openssl you may need to set the following extra environment variables to install M2Crypto:

.. code-block:: bash

   CFLAGS="-I/usr/local/opt/openssl/include -L/usr/local/opt/openssl/lib" SWIG_FEATURES="-cpperraswarn -includeall -I/usr/local/opt/openssl/include/" pip install M2Crypto

Once you have these packages installed, you can now install lalsuite following the instructions at:

.. toctree::
    :maxdepth: 1

    install_lalsuite

=========================================
Installing PyCBC in a virtual environment
=========================================

Enter the virtual environment that you wish to use for PyCBC development by sourcing the activate script, for example

.. code-block:: bash

    source ~/src/pycbc/bin/activate

Next install the Pegasus WMS python libraries needed to build the workflows with the command:

.. code-block:: bash

    pip install http://download.pegasus.isi.edu/pegasus/4.7.3/pegasus-python-source-4.7.3.tar.gz

To query the new Advanced LIGO and Advanced Virgo Segment Database, you will need to install the ``dqsegdb`` tools. Install the 1.4.1 pre-release of these tools, run the command:

.. code-block:: bash

    pip install dqsegdb

For uploading triggers to GraceDB at the end of the workflow you will need to have the gracedb client tools installed. The latest release is in pip

.. code-block:: bash

    pip install ligo-gracedb

You now need to decide whether you want to install a release of PyCBC or an editable version of the source code from a git repository for development. 

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Installing a released version of PyCBC
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

.. note::

    Make sure you have run the commands in the section :ref:`installinglalsuite` above to install unittest2 before installing PyCBC.

Releases of PyCBC are available from the `PyPi PyCBC page <https://pypi.python.org/pypi/PyCBC/1.6.0>`_. To install the latest release run the command

.. code-block:: bash

    pip install PyCBC

To install an older version, use the `pip version specifier <https://packaging.python.org/glossary/#term-version-specifier>`_.

To install a release of the code from GitHub, determine the tag of the release that you want to install from the `list of PyCBC tags <https://github.com/ligo-cbc/pycbc/tags>`_. This example installs the v1.1.0 release. If you want to install a different release, change the command below accordingly:

.. code-block:: bash

    pip install git+https://github.com/ligo-cbc/pycbc@v1.1.0#egg=pycbc

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Installing source from GitHub for development
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

.. note::

    Make sure you have run the command in the section :ref:`installinglalsuite` above to install unittest2 before installing PyCBC.

To install and editable version of PyCBC you need to have `forked PyCBC to your own account <https://help.github.com/articles/fork-a-repo/>`_ and know the URL of your fork. This can be obtained from the clone URL on your GitHub repository page. This example uses the URL ``git@github.com:your-username-here/pycbc.git`` which you should change to the URL of your fork of PyCBC on GitHub. 

.. note:: 
    There are two main authentication schemes for GitHub: SSH and HTTPS. The examples below use URLs containing ``git@github.com``, which assumed that you are usingg SSH authentication. If you have not already enabled ssh keys on your GitHub account, you can follow the `GitHub instructions for setting up ssh keys <https://help.github.com/articles/generating-ssh-keys/>`_ to set this up. Alternatively, you can use the HTTPS connection method, where ``git@github.com:`` is replaced with ``https://github.com/``. See the `GitHub documentation on URLs <https://help.github.com/articles/which-remote-url-should-i-use/>`_ for more information. You can also read the `pip git instructions <https://pip.pypa.io/en/latest/reference/pip_install.html#git>`_ for more details on how to install a branch or a specific tag.

Install the PyCBC source code from the GitHub URL using the command:

.. code-block:: bash

    pip install -e git+git@github.com:your-username-here/pycbc.git#egg=pycbc

This will fetch the PyCBC source and will also install all the listed dependencies. The ``-e`` option to pip creates a directory called ``${VIRTUAL_ENV}/src/pycbc`` with a git checkout which is fully editable. To prevent pip from removing this source directory run the command

.. code-block:: bash

    rm -f ${VIRTUAL_ENV}/src/pip-delete-this-directory.txt

You can then make changes to your PyCBC source code in the directory ``${VIRTUAL_ENV}/src/pycbc``

You can also use the repository created by pip as your working repository, creating branches, commits, and `pull requests <https://help.github.com/articles/creating-a-pull-request/>`_ as you need to. To keep your repository in sync with the ligo-cbc/pycbc repository, you can read the GitHub instructions that explain how to `sync a fork of a repository to keep it up-to-date with the upstream repository. <https://help.github.com/articles/syncing-a-fork/>`_.

.. note::

    The version of PyCBC that is checked out will be on the master branch. To track a remote branch from your GitHub repository, run the command
    
    ``git branch --set-upstream branch_name origin/branch_name``

    ``git checkout branch_name``
        
    where branch_name is the name of the branch that you want to track.

To build and install any changes that you make to the source code in your virtual environment, run the command

.. code-block:: bash

    python setup.py install
    
from the PyCBC source directory in ``${VIRTUAL_ENV}/src/pycbc``

.. note::

    PyCBC no longer requires an installation of the ``pycbc-pylal`` package.
    If you are running a workflow that needs code from this package, you will
    need to manually install it with ``pip install pycbc-pylal``.


=====================================
Building and Installing Documentation
=====================================

To build the documentation from your virtual environment, first make sure that you have `Sphinx <http://sphinx-doc.org/>`_ and the required helper tools installed with

.. code-block:: bash

    pip install "Sphinx>=1.5.0"
    pip install sphinx-rtd-theme
    pip install git+https://github.com/ligo-cbc/sphinxcontrib-programoutput.git#egg=sphinxcontrib-programoutput
    
To generate the documentation and push it to your personal GitHub pages, first create a branch names ``gh-pages``, if you do not already have one. Follow the `GitHub branch <https://help.github.com/articles/creating-and-deleting-branches-within-your-repository/>`_ instructions to do this.

To build and publish the documentation, run the following commands from the
top-level of your PyCBC source tree, replacing ``github-username`` with your
GitHub user name:

.. code-block:: bash

    git clone git@github.com:github-username/pycbc.git _gh-pages
    cd _gh-pages
    git checkout gh-pages
    git rm -rf *
    git commit -a -m "flush documentation"
    cd ..
    python setup.py build_gh_pages
    cd _gh-pages
    git add --all
    git commit -a -m "documentation update"
    git push origin gh-pages

The documentation will then be visible at http://github-username.github.io/pycbc/latest/html where ``github-username`` should be replaced with your GitHub username.

.. note::

    Be careful with the ``git rm -rf *`` command as if you run it in the wrong
    directory you can delete the contents of your git repository. If you do
    this by accident, you can use ``git reset`` to undo the commit.

For more details on building and maintaining the documentation under GitHub project pages, see

.. toctree::
    :maxdepth: 1

    build_gh_pages


====================
Modifying pycbc-glue
====================

PyCBC depends on the package pycbc-glue which is a fork of the lalsuite packages. The correct version is automatically installed by pip as part of the main PyCBC install. If you are developing code in pycbc-glue, then you can clone them from GitHib into your virtual environment's source directory and build and install them from there.

.. note::

    If you want to develop pycbc-glue, you should follow the instructions to `fork a repository <https://help.github.com/articles/fork-a-repo/>`_ to fork the `ligo-cbc/pycbc-glue <https://github.com/ligo-cbc/pycbc-glue>`_ repository into your own account.

You can obtain these repositories in the standard way using git, replacing ``ligo-cbc`` with your GitHub user account name

.. code-block:: bash

    cd ${VIRTUAL_ENV}/src
    git clone git@github.com:ligo-cbc/pycbc-glue.git

Once you have the source code cloned, you can run 

.. code-block:: bash

    python setup.py install

to install pycbc-glue into your virtual environment.

========================================
Use of Intel MKL Optimized FFT libraries
========================================

PyCBC has the ability to use optimized FFT libraries such as FFTW and MKL. If MKL is the correct library for your platform, you can add the script that sets up the MKL environment to you virtualenv ``activate`` script with the command

.. code-block:: bash

    echo 'source /opt/intel/bin/compilervars.sh intel64' >> ${VIRTUAL_ENV}/bin/activate

changing the path to the ``compilervars.sh`` script appropriately for your cluster. 

==========================================
Graphics Processing Unit support with CUDA
==========================================

PyCBC has the ability to accelerate its processing using CUDA. To take advantage of this, follow the instructions linked below to install the CUDA dependencies.

.. toctree::
    :maxdepth: 1

    install_cuda
