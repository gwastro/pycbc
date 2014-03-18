=============================================
Instructions to add OpenCL support (optional)
=============================================

If you would like to use GPU acceleration of PyCBC through OpenCL you will require these additional packages:

* `NumPy <http://www.numpy.org>`_ >= 1.5.1
* `OpenCL drivers (various sources)`_ 
* PyFFT <https://pypi.python.org/pypi/pyfft>`_ >= 0.3.9
* `PyOpenCL <http://mathema.tician.de/software/pyopencl>`_ >= 2013.1.1
* `Mako <http://www.makotemplates.org/>`_ >= 0.7.2

These packages may not be available via the distribution packaging system, at least in the required versions. As described below, most of these packages are available via the python package installer `pip <http://www.pip-installer.org>`_, however custom installation instructions are given where required.

.. note::
    When using ``pip`` as described below, to install system-wide, simply do not give the ``--prefix`` option (although you might need to add ``sudo`` as a command prefix).

------
PyOpenCL
------

PyOpenCL should be installed from source, so that the latest updates are applied:

.. code-block:: bash

    git clone http://git.tiker.net/trees/pyopencl.git
    cd pyopencl
    git submodule init
    git submodule update
    ./configure.py
    python setup.py build
    python setup.py install --prefix=/path/to/install/location
    
------------
pyfft
------------

.. code-block:: bash

   pip install pyfft --prefix=/path/to/install/location

----
Mako
----

.. code-block:: bash

   pip install Mako --prefix=/path/to/install/location
   
-------------------------------
Setting up the user environment
-------------------------------

Make sure that the chosen install locations are added to your PYTHONPATH variable.
