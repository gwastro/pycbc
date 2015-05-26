################
Installing PyCBC
################

=============
Prerequisites
=============

In order to install PyCBC, you need to have installed the following prerequisite packages:

* Python 2.6 or 2.7
* `NumPy <http://www.numpy.org>`_ >= 1.4.1 and `SciPy <http://www.scipy.org>`_ >= 0.7.2
* `decorator <https://pypi.python.org/pypi/decorator>`_
* `argparse <https://pypi.python.org/pypi/argparse>`_ >= 1.2.0
* `LALSuite <https://www.lsc-group.phys.uwm.edu/daswg/projects/lalsuite.html>`_ (with swig bindings enabled)
* `GLUE <https://www.lsc-group.phys.uwm.edu/daswg/projects/glue.html>`_
* `pylal <https://www.lsc-group.phys.uwm.edu/daswg/projects/pylal.html>`_

.. note::
    
    Python, numpy, decorator and argparse should already be installed on LDG clusters.

.. note::

    A version of lalsuite and glue are installed on LDG clusters, but you may want to build your own version. Please see :ref:`lalsuite_install` for instructions and details about building your own version of lalsuite, glue and pylal.

===========================================================
Additional Dependencies for HDF post processing and Plots
===========================================================
In order to run the HDF post processing and plotting codes
that are in active development and have not yet been reviewed, the following
additional dependencies are needed. Eventually these will become
mandatory dependencies. 

* numpy>=1.6.4
* matplotlib>=1.3.1
* `mpld3>=0.3.0 <https://github.com/jakevdp/mpld3/tarball/master>`_
* `PIL <https://pypi.python.org/pypi/PIL>`_

===================
Installing from git
===================

The source for PyCBC is under ``git`` version control, hosted on `github <https://github.com/ligo-cbc/pycbc>`_

You can install the package by first cloning the repository, either read-only:

.. code-block:: bash

    git clone https://github.com/ligo-cbc/pycbc.git


You can specify the install path directory, using the ``--prefix`` option as follows.

.. code-block:: bash

    python setup.py install --prefix=/location/to/install/pycbc
    
Alternatively, you can then run ``setup.py`` with the ``--user`` option to install the package in the default user location:

.. code-block:: bash

    cd pycbc
    python setup.py install --user

The ``--user`` option tells the installer to copy codes into the standard user library paths, on linux machines this is

.. code-block:: bash

    ~/.local/lib

while on Mac OS this is

.. code-block:: bash

    ~/Library/Python/X.Y/lib

where ``X.Y`` is the python major and minor version numbers, e.g. ``2.7``. In either case, python will autmatically know about these directories, so you don't have to fiddle with any environment variables.


===============================
Setting up the user environment
===============================

Add the following to your ``.bash_profile``

.. code-block:: bash

   source /path/to/pycbc/install/directory/etc/pycbc-user-env.sh
   
===============================
Optional GPU acceleration
===============================
PyCBC has the ability to accelerate its processing using either CUDA or
OpenCL. 

.. toctree::
    :maxdepth: 1

    cuda_install
    opencl_install
