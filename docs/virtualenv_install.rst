################
Setting up pip and virtualenv
################

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

