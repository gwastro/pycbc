=============================================
Instructions to add CUDA support (optional)
=============================================

If you would like to use GPU acceleration of PyCBC through CUDA you will require these additional packages:

* `NumPy <http://www.numpy.org>`_ >= 1.5.1
* `Nvidia CUDA <http://www.nvidia.com/object/cuda_home_new.html>`_ >= 4.0 (driver and libraries)
* `PyCUDA <http://mathema.tician.de/software/pycuda>`_ >= 2013.1.1
* `SciKits.cuda <http://scikits.appspot.com/cuda>`_ >= 0.041
* `Mako <http://www.makotemplates.org/>`_ >= 0.7.2

These packages may not be available via the distribution packaging system, at least in the required versions. As described below, most of these packages are available via the python package installer `pip <http://www.pip-installer.org>`_, however custom installation instructions are given where required.

.. note::
    When using ``pip`` as described below, to install system-wide, simply do not give the ``--prefix`` option (although you might need to add ``sudo`` as a command prefix).

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
    python setup.py install --prefix=/path/to/install/location

If your CUDA installation is in a non-standard location X, pass ``-–cuda-root=X`` to ``configure.py``.

------------
SciKits.cuda
------------

.. code-block:: bash

   pip install scikits.cuda --prefix=/path/to/install/location

----
Mako
----

.. code-block:: bash

   pip install Mako --prefix=/path/to/install/location
   
-------------------------------
Setting up the user environment
-------------------------------

Make sure that the chosen install locations are added to your PYTHONPATH variable.
