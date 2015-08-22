################
Installing PyCBC
################

The reccomended way of installing PyCBC is to use `pip <https://pip.pypa.io/en/stable/>`_ within a `Python Virtual Envionment <https://virtualenv.pypa.io/en/latest/>`_. Virtualenv isolates PyCBC and its dependencies from the system environment and installing with pip ensures that PyCBC picks up the correct dependencies. 

There are three typical use cases for PyCBC

1. Installing a release of PyCBC from GitHub for an end user to run the tools.
2. Installing an editable version from GitHub for development.
3. Production LIGO analyses.

This page documents the first two use cases. For production analysis, users must obtain the pre-built binaries from the PyCBC server. 

If you wish to develop PyCBC, then you will need an account on `GitHub <https://www.github.com>`_ Once you have set up your account you should follow the instructions to `fork a repository <https://help.github.com/articles/fork-a-repo/>`_ to fork the `ligo-cbc/pycbc <https://github.com/ligo-cbc/pycbc>`_ repository into your own account.

=============
Setting up pip and virtualenv
=============

LIGO Data Grid users on Scientic Linux 6 systems need to perform a one-time setup of pip and virtualenv, as the versions installed on Scientic Linux 6 are too old for use with PyCBC. Note that once this part of the installation has been performed, you should not need to do it again.

First install pip in your home directory. In this example, we install in the directory ${HOME}/local/

.. code-block:: bash

    mkdir -p ${HOME}/local/pip-7.1.0/lib/python2.6/site-packages
    export PYTHONPATH=${HOME}/local/pip-7.1.0/lib/python2.6/site-packages
    export PATH=${HOME}/local/pip-7.1.0/bin:${PATH}
    
Note that when setting PYTHONPATH we have excluded all other directories from the path. This is to prevent the install process from picking up any existing libraries that may be in your PYTHONPATH. When setting PATH, the pip install directory is placed at the start of the PATH so that it is found first by the shell.

Now use the system version of easy_install to install pip 7.1.0 in this directory with the command

.. code-block:: bash

    easy_install --prefix=${HOME}/local/pip-7.1.0 https://pypi.python.org/packages/source/p/pip/pip-7.1.0.tar.gz#md5=d935ee9146074b1d3f26c5f0acfd120e

Next check that you are using the correct version of pip by running the command

.. code-block:: bash

    pip --version
    
This should report pip 7.1.0. If it returns an older version, check that you set your PATH and PYTHONPATH correctly. 

Finally install virtual env using the version of pip that you just inst

.. code-block:: bash

    pip install virtualenv --upgrade --user
    
This installs virtualenv into your ${HOME}/.local directory where pip installs user packages. (Note the period at the start of .local).

Add virtualenv to your path by running the command

.. code-block:: bash

    export PATH=$HOME/.local/bin:$PATH
    
You may want to add this command to your .bash_profile so that virtualenv is available when you log in.

===========================
Creating a virtualenv
===========================

Installing PyCBC into a virtual environment provides isolation between different sets of python packages. The following instructions will create a working PyCBC environment on an LDG cluster. 

Make sure that you have at least version 13.1.1 of virtualenv by running 

.. code-block:: bash

    virtualenv --version
    
If this returns virtualenv: command not found or a lower version, see the instructions above for installing virtualenv.

The first task is to clear out your current PYTHONPATH in case you had been using that before. To do this, run the command:

.. code-block:: bash

    unset PYTHONPATH

By default, virtualenv will modify your shell prompt so that it prepends the name of the virtual environment. This can be useful to make sure that you are developing in the virtual environment, or if you have several virtual environments. However, if you do not want this, then set

.. code-block:: bash

    VIRTUAL_ENV_DISABLE_PROMPT=True
    
Before running the command to create the new virtual environment.

Next, you need to choose a directory name where you'd like to make your virtual environment, and then make it. In this example, we use ${HOME}/pycbc-dev but this can be changed to any path other than a directory under ${HOME}/.local: 

.. code-block:: bash

    NAME=${HOME}/pycbc-dev
    virtualenv $NAME
    
To enter your virtual environment run the command

.. code-block:: bash
    
    source $NAME/bin/activate
    
You will now be 'inside' your virtual environment, and so you can install packages, etc, without conflicting with either the system build, or other builds that you may have sitting around. You may install c-dependencies such as lalsuite (:ref:`lalsuite_install`), or rely on the system versions.

To leave this virtual environment type

.. code-block:: bash

    deactivate
    
which will return you to a regular shell.

===========================
Installing PyCBC in a virtualenv
===========================

Enter the virtual enviornment that you wish to use for PyCBC development by sourcing the activate script, as shown in the previous section.

Install pycbc from source as follows. First install unittest2 and numpy with the command:

.. code-block:: bash

    pip install "numpy>=1.6.4" unittest2
    
You now need to decide whether you want to install a release of PyCBC for end-use, or an editable git repository for development. 

To install a release of the code, determine the tag of the relase that you want to install from the `list of PyCBC tags <https://github.com/ligo-cbc/pycbc/tags>`_. This example installs the v1.1.0 release. If you want to install a different release, change the command below accordingly:

.. code-block:: bash

    pip install git+https://github.com/ligo-cbc/pycbc@v1.1.0#egg=pycbc --process-dependency-links

To install and editable version of PyCBC you need to have `forked PyCBC to your own account <https://help.github.com/articles/fork-a-repo/>`_ and know the URL of your fork. This can be obtained from the clone URL on your GitHub repository page. This example uses the URL git@github.com:duncan-brown/pycbc.git which you should change as appropriate. You can read the `pip git instructions <https://pip.pypa.io/en/latest/reference/pip_install.html#git>`_ for more details on how to install a branch or a specific tag.

.. code-block:: bash

    pip install -e git+git@github.com:duncan-brown/pycbc.git#egg=pycbc --process-dependency-links

This will create a directory called $NAME/src/pycbc git checkout which is fully edittable, and will also install all the listed dependencies.

You may now enter your environement by sourcing 'activate', and leave by running
the command 'deactivate'.


===============================
Optional GPU acceleration
===============================
PyCBC has the ability to accelerate its processing using CUDA. 

.. toctree::
    :maxdepth: 1

    cuda_install
