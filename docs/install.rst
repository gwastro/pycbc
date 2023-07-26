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

    pip install pycbc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Full Virtualenv for Development and Production
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This document explains how to set up a virtual environment to install PyCBC
either for development or use in a production analysis with a release. The
code build will be a standard Python install, which requires that the
installation directory containing the Python libraries is accessible at
runtime.

PyCBC uses the `fork and pull <https://help.github.com/articles/using-pull-requests/>`_ model for development. If you wish to develop PyCBC, then you will need an account on `GitHub <https://www.github.com>`_. Once you have set up your account you should follow the instructions to `fork a repository <https://help.github.com/articles/fork-a-repo/>`_ to fork the `gwastro/pycbc <https://github.com/gwastro/pycbc>`_ repository into your own account. From your own fork, you can follow the `GitHub flow model <https://help.github.com/articles/github-flow-in-the-browser/>`_ to develop and maintain the code. For each new feature or bug fix, you should `create a new branch <https://help.github.com/articles/creating-and-deleting-branches-within-your-repository/>`_ to develop the feature. You can then `create a pull request <https://help.github.com/articles/creating-a-pull-request/>`_ so that the PyCBC maintainers can review and merge your changes into the official repository.


Create a virtual environment for your development or production environment and
do a clean install. The following will create a python3 environment: currently PyCBC requires python3.7 or higher.

.. code-block:: bash

    virtualenv -p python3 env
    source env/bin/activate
    pip install --upgrade pip setuptools

We can then make a fresh clone of the repository.

.. code-block:: bash

    git clone git@github.com:gwastro/pycbc.git

Finally, install the most common pycbc develeopment environment packages
as follows.

.. code-block:: bash

    cd pycbc
    pip install -r requirements.txt
    pip install -r companion.txt
    pip install .

========================================
Development build on LDG / IGWN clusters
========================================

The above instructions require some adjustment when working on a LIGO or other GW collaboration
compute cluster (eg CIT).  The main issue is that the default environment may not include a
sufficiently recent python version (>=3.7).  The standard workaround is to use a python executable
available in a 'IGWN Conda' environment.  To see what environments are available, you can run

.. code-block:: bash

    conda info --envs

This should yield ``igwn-py37`` as one choice.  The output of this command will also
tell you the location of the environment in the file system.  Then, the location of the
python3.7 executable is for instance ``/cvmfs/oasis.opensciencegrid.org/ligo/sw/conda/envs/igwn-py37/bin/python``
and you will create the virtualenv via the command

.. code-block:: bash

    virtualenv -p /cvmfs/oasis.opensciencegrid.org/ligo/sw/conda/envs/igwn-py37/bin/python env

Once the virtualenv has been created you can install PyCBC from PyPI or a local
copy with the `[igwn]` extra specifier to install the optional extra requirements
recommended for IGWN users:

.. code-block:: bash

    pip install -r requirements-igwn.txt
    pip install .[igwn]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Other scenarios
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

===========
Docker
===========
For more details, including instructions on starting a docker container  see:

    .. toctree::
       :maxdepth: 1

       docker

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

==============
Related Help
==============

.. toctree::
   :maxdepth: 1

   install_virtualenv
   install_lalsuite
