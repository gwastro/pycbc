===========================================
Instructions to add CUDA support (optional)
===========================================

If you would like to use GPU acceleration of PyCBC through CUDA you will require these additional packages:

* `NumPy <http://www.numpy.org>`_ >= 1.5.1
* `Nvidia CUDA <http://www.nvidia.com/object/cuda_home_new.html>`_ >= 6.5 (driver and libraries)
* `PyCUDA <http://mathema.tician.de/software/pycuda>`_ >= 2015.1.3
* `SciKits.cuda <http://scikits.appspot.com/cuda>`_ >= 0.041
* `Mako <http://www.makotemplates.org/>`_ >= 0.7.2

These packages may not be available via the distribution packaging system, at least in the required versions. As described below, most of these packages are available via the python package installer `pip <http://www.pip-installer.org>`_, however custom installation instructions are given where required.

If you are currently in your virtual environment, leave it by running ``deactivate`` as you need to add some additional environment variables before continuing.  Set the shell variable ``NAME`` to the location of your virtual environment. Here we assume that your virtual environment is installed in ``${HOME}/pycbc-dev``. If it is in different location, you will need to change this as appropriate.

.. code-block:: bash

    NAME=${HOME}/src/pycbc

The install requires that you set the environment variable ``CUDA_ROOT``, make sure that the CUDA ``bin`` directory is in your path, and add the CUDA library path to your ``LD_LIBRARY_PATH``. You can do this by adding these commands to your ``activate`` script by running the commands:

.. code-block:: bash

    echo 'export CUDA_ROOT=/usr/local/cuda' >> $NAME/bin/activate
    echo 'export PATH=${CUDA_ROOT}/bin:${PATH}' >> $NAME/bin/activate
    echo 'export LD_LIBRARY_PATH=${CUDA_ROOT}/lib64:${LD_LIBRARY_PATH}' >> $NAME/bin/activate

Now activate your virtual environment. 

.. code-block:: bash

    source ${NAME}/bin/activate

--------------------------------
Installing the CUDA dependencies
--------------------------------

Install the dependencies PyCUDA, SciKits.cuda and Mako with by running the commands

.. code-block:: bash

    pip install pycuda
    pip install scikit-cuda
    pip install Mako
   
You should now be able to use the CUDA features in PyCBC.
