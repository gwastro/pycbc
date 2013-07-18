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
* `LAL <https://www.lsc-group.phys.uwm.edu/daswg/projects/lalsuite.html>`_
* `GLUE <https://www.lsc-group.phys.uwm.edu/daswg/projects/glue.html>`_

If you would like to use GPU acceleration of PyCBC through CUDA you will require these additional packages:

* `NumPy <http://www.numpy.org>`_ >= 1.5.1
* `Nvidia CUDA <http://www.nvidia.com/object/cuda_home_new.html>`_ >= 4.0 (driver and libraries)
* `PyCUDA <http://mathema.tician.de/software/pycuda>`_ >= 2012.1
* `SciKits.cuda <http://scikits.appspot.com/cuda>`_ >= 0.041
* `Mako <http://www.makotemplates.org/>`_ >= 0.7.2

==============================
Installing CUDA Python modules
==============================

These packages may not be available via the distribution packaging system, at least in the required versions. As described below, most of these packages are available via the python package installer `pip <http://www.pip-installer.org>`_, however custom installation instructions are given where required.

.. note::
    When using ``pip`` as described below, to install system-wide, simply do not give the ``--user`` option (although you might need to add ``sudo`` as a command prefix).

------
PyCUDA
------

PyCUDA should be installed from source, so that the latest updates are applied:

.. code-block:: bash

    git clone http://git.tiker.net/trees/pycuda.git
    cd pycuda
    git submodule init
    git submodule update
    ./configure.py
    python setup.py build
    python setup.py install --user

If your CUDA installation is in a non-standard location X, pass ``–cuda-root=X`` to ``configure.py``.

------------
SciKits.cuda
------------

.. code-block:: bash

   pip install scikits.cuda --user

----
Mako
----

.. code-block:: bash

   pip install Mako --user

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

You can then run ``setup.py`` to install the package:

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
