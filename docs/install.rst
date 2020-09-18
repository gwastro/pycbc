.. _installing_pycbc:

################
Installing PyCBC
################

%%%%%%%%%%%%%%%%%%%%%%%%
Simple Installation
%%%%%%%%%%%%%%%%%%%%%%%%

PyCBC is available through the PyPI. For straightforward use of the PyCBC library
and executables, we recommend installing with the following command. If you
are not running in a specialized computing environment, this is probably the
appropriate thing to do.

.. code-block:: bash

    pip install lalsuite pycbc

PyCBC depends on `lalsuite` for a lot of functionality, however, if you are
getting lalsuite through another means, you may ommit this part of the command.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Full Virtualenv for Development and Production
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This document explains how to set up a virtual environment to install PyCBC
either for development or use in a production analysis with a release. The
code build will be a standard Python install which requires that the
installation directory containing the Python libraries is accessible at
runtime.

PyCBC uses the `fork and pull <https://help.github.com/articles/using-pull-requests/>`_ model for development. If you wish to develop PyCBC, then you will need an account on `GitHub <https://www.github.com>`_. Once you have set up your account you should follow the instructions to `fork a repository <https://help.github.com/articles/fork-a-repo/>`_ to fork the `gwastro/pycbc <https://github.com/gwastro/pycbc>`_ repository into your own account. From your own fork, you can follow the `GitHub flow model <https://help.github.com/articles/github-flow-in-the-browser/>`_ to develop and maintain the code. For each new feature or bug fix, you should `create a new branch <https://help.github.com/articles/creating-and-deleting-branches-within-your-repository/>`_ to develop the feature. You can then `create a pull request <https://help.github.com/articles/creating-a-pull-request/>`_ so that the PyCBC maintainers can review and merge your changes into the official repository.


Create a virtual environment for your development or production environment and
do a clean install. The following will create a python3 environment.

.. code-block:: bash

    virtualenv -p python3 env
    source env/bin/activate
    pip install --upgrade pip setuptools

We can then make a fresh clone of the repository.

.. code-block:: bash

    git clone git@github.com/gwastro/pycbc.git

Finally, install the most common pycbc develeopment environment packages
as follows.

.. code-block:: bash

    cd pycbc
    pip install -r requirements.txt
    pip install -r companion.txt
    python setup.py install

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Other scenarios
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

=====================================
Building the Documentation
=====================================

To build the documentation from your virtual environment, first make sure that you have `Sphinx <http://sphinx-doc.org/>`_ and the required helper tools installed with

.. code-block:: bash

    pip install "Sphinx>=1.5.0"
    pip install sphinx-rtd-theme
    pip install sphinxcontrib-programoutput

You can then build the documentation locally as

.. code-block:: bash

    python setup.py build_docs

The documentation will show up locally in 'docs/_build/html'.

========================================
Use of Intel MKL Optimized FFT libraries
========================================

PyCBC has the ability to use optimized FFT libraries such as FFTW and MKL. If MKL is the correct library for your platform,
you can add the script that sets up the MKL environment to you virtualenv ``activate`` script with the command. Typically,
the MKl source command appears as follows but may vary depending on your cluster / environment.

.. code-block:: bash

    echo 'source /opt/intel/bin/compilervars.sh intel64'

==========================================
Graphics Processing Unit support with CUDA
==========================================

PyCBC has the ability to accelerate its processing using CUDA. To take advantage of this, follow the instructions linked below to install the CUDA dependencies.

.. toctree::
    :maxdepth: 1

    install_cuda
