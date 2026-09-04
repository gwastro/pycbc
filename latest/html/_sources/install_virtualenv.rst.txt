#####################
Setting up virtualenv
#####################

.. note::

    LIGO Data Grid users on Scientic Linux 6 systems will need to perform a one-time setup of pip and virtualenv, as the versions installed on Scientic Linux 6 are too old for use with PyCBC. Note that once this part of the installation has been performed, you should not need to do it again.

You need at least version 7.1.0 of pip to manage the installation process for virtualenv. Run the command

.. code-block:: bash

    pip --version

to determine the version of pip that you have installed. If you have version 7.1.0 or greater, you can continue to the section `Install virtualenv`_ without installing pip. If this command returns ``pip: command not found`` or a version of pip less than 7.1.0, then follow the instruction in the section to `Install pip`_

================
Install pip
================

First install pip in your home directory. In this example, we install in the directory ``${HOME}/local/``. Run the commands

.. code-block:: bash

    VERSION=`python -c 'import sys; print "%d.%d" % (sys.version_info[0], sys.version_info[1])'`
    mkdir -p ${HOME}/local/pip-7.1.0/lib/python${VERSION}/site-packages
    export PYTHONPATH=${HOME}/local/pip-7.1.0/lib/python${VERSION}/site-packages
    export PATH=${HOME}/local/pip-7.1.0/bin:${PATH}
    
to set up your environment to install pip. Note that when setting ``PYTHONPATH`` we have excluded all other directories from the path. This is to prevent the install process from picking up any existing libraries that may be in your ``PYTHONPATH``. When setting ``PATH``, the pip install directory is placed at the start of the ``PATH`` so that it is found first by the shell.

Now use the system version of easy_install to install pip 7.1.0 in this directory with the command

.. code-block:: bash

    easy_install --prefix=${HOME}/local/pip-7.1.0 https://pypi.python.org/packages/source/p/pip/pip-7.1.0.tar.gz#md5=d935ee9146074b1d3f26c5f0acfd120e

Next check that you are using the correct version of pip by running the command

.. code-block:: bash

    pip --version
    
This should report pip 7.1.0. If it returns an older version, check that you set your ``PATH`` and ``PYTHONPATH`` correctly. 

==================
Install virtualenv
==================

Install virtual env using the version of pip that you just inst

.. code-block:: bash

    pip install virtualenv --upgrade --user
    
This installs virtualenv into your ``${HOME}/.local`` directory where pip installs user packages. (Note the period at the start of ``.local``).

.. note:: 

    If this command returns ``Can not perform a '--user' install. User site-packages are not visible in this virtualenv.`` then you are already in a Python virtual environment. You will need to deactivate this environment or contact your system administrator to install virtualenv for you in this environment.

Add virtualenv to your path by running the command

.. code-block:: bash

    export PATH=${HOME}/.local/bin:${PATH}
    
You may want to add this command to your ``.bash_profile`` so that virtualenv is available when you log in by running

.. code-block:: bash

    echo 'export PATH=${HOME}/.local/bin:${PATH}' >> ${HOME}/.bash_profile

You can now continue installing PyCBC.
