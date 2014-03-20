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

.. note::
    
    These dependencies should already be installed on LDG clusters. 

===================
Installing from git
===================

The source for PyCBC is under ``git`` version control, hosted by the University of Wisconsin-Milwaukee.

You can install the package by first cloning the repository, either read-only:

.. code-block:: bash

    git clone git://ligo-vcs.phys.uwm.edu/pycbc.git

or with a LIGO.ORG albert.einstein-style credential:

.. code-block:: bash

    git clone albert.einstein@ligo-vcs.phys.uwm.edu:/usr/local/git/pycbc.git

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
