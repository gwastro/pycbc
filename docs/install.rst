################
Installing PyCBC
################

=============
Prerequisites
=============

In order to install PyCBC, you need to have installed the following prerequisite packages:

* Python 2.6 or 2.7
* `LALSuite <https://www.lsc-group.phys.uwm.edu/daswg/projects/lalsuite.html>`_ (with swig bindings enabled)
* `NumPy <http://www.numpy.org>`_ >= 1.6.4 and `SciPy <http://www.scipy.org>`_ >= 0.13.0
* `decorator <https://pypi.python.org/pypi/decorator>`_ >=3.4.2
* `argparse <https://pypi.python.org/pypi/argparse>`_ >= 1.3.0
* `pycbc-glue <https://github.com/ligo-cbc/pycbc-glue>`_ >=0.9.1
* `pycbc-pylal <https://github.com/ligo-cbc/pycbc-glue>`_ >=0.9.2
* `jinja2 <https://pypi.python.org/pypi/jinja2>`_
* `mako <https://pypi.python.org/pypi/mako>`_
* `h5py <https://pypi.python.org/pypi/h5py>`_ >=2.5

.. note::

    A version of lalsuite is installed on LDG clusters, but you may want to build your own version. Please see :ref:`lalsuite_install` for instructions and details about building your own version of lalsuite.

===========================================================
Additional Dependencies for HDF post processing and Plots
===========================================================
In order to run the HDF post processing and plotting codes
that are in active development and have not yet been reviewed, the following
additional dependencies are needed. Eventually these will become
mandatory dependencies. 

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


Setting up the user environment
********************************

Add the following to your ``.bash_profile``

.. code-block:: bash

   source /path/to/pycbc/install/directory/etc/pycbc-user-env.sh
   
===========================
Installing in a virtualenv
===========================

Installing PyCBC into a virtual environment provides isolation between different sets of 
python packages. The following instructions will create a working PyCBC environment on an LDG
cluster. 

The first task is to clear out your current PYTHONPATH in case you had been using that
before, and to make sure you have the latest virtualenv code installed.:

.. code-block:: bash

    export PYTHONPATH=””
    pip install virtualenv --upgrade --user
    export PATH=$HOME/.local/bin:$PATH
    
The previous step may be skipped if you do not have conflicting python packages in 
your PYTHONPATH.

Next, you need to choose a directory name where you'd like to make your virtual environment, and
then make it.: 

.. code-block:: bash

    virtualenv $NAME
    source $NAME/bin/activate
    
You will now be 'inside' your virtual environment, and so you can install packages, etc, without
conflicting with either the system build, or other builds that you may have sitting around. You may
install c-dependencies such as lalsuite (:ref:`lalsuite_install`), or rely on the system versions.

Install pycbc from source as follows. This will create a $NAME/src/pycbc git checkout
which is fully edittable, and will also install all the listed dependencies.:

.. code-block:: bash

    pip install “numpy>=1.6.4” unittest2
    pip install -e git+https://github.com/ligo-cbc/pycbc#egg=pycbc --process-dependency-links

Finally, if you need to use system python libraries (such as Pegasus.DAX3), etc you can
add that to your own activation script.:

.. code-block:: bash

    echo 'source $NAME/bin/activate' > activate
    echo 'source $NAME/etc/pycbc-user-env.sh' >> activate
    echo 'source $NAME/etc/glue-user-env.sh' >> activate
    echo 'export PYTHONPATH=/usr/lib64/python2.6/site-packages/:$PYTHONPATH' >> activate
    chmod 755 activate
    
You may now enter your environement by sourcing 'activate', and leave by running
the command 'deactivate'.


===============================
Optional GPU acceleration
===============================
PyCBC has the ability to accelerate its processing using CUDA. 

.. toctree::
    :maxdepth: 1

    cuda_install
