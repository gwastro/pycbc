##########################
Creating Releases of PyCBC
##########################

.. note::

    Only authorized maintainers of PyCBC should use these instructions.


==============================
Creating the release on GitHub
==============================

To create a new PyCBC release:

#. Make sure that the setup.py file contains the correct version number, which should be in the format ``x.y.x`` (where x, y, and z are the major, minor, and patch levels of the release) in PyCBC's setup.py file.
#. Set ``release = True`` in the PyCBC's setup.py file.
#. Commit the changed setup.py file and push to commit to the repository.
#. Go to the `PyCBC release page <https://github.com/ligo-cbc/pycbc/releases>`_ and click on ``Draft a new release``. 
#. Enter a tag version in the format ``vx.y.z``. Note the ``v`` in front of the major, minor, and patch numbers. 
#. Enter a descriptive release title and write a desciption of the release in the text box provided.
#. Click on ``Publish release`` to create the release.
#. Update the setup.py file with an incremented major or minor version number and make sure that the string ``dev`` appears in that version. For example, if you just released ``1.2.1`` then change the string to ``1.3.dev0`` or ``2.0.dev0`` as appropriate. This is need to ensure that if someone is building from source, it always takes precedence over an older release version.
#. Set ``release = False`` in PyCBC's setup.py file.
#. Commit these changes and push them to the repository.

=====================================
Uploading to the Python Package Index
=====================================

.. note::

    Keep in mind that once you register a version number with the package index and upload a package, you can never change it. You can only deactivate the package from the index, and increment the version number.

If this is the first time that you have pushed to PyPI, you will need to create a configuration file. Create a file in your home directory called ``.pypirc`` that contains the lines::
    [distutils] # this tells distutils what package indexes you can push to
    index-servers =
    pypi
    pypitest

    [pypi]
    repository: https://pypi.python.org/pypi
    username: your_username
    password: your_password

    [pypitest]
    repository: https://testpypi.python.org/pypi
    username: your_username
    password: your_password

Replace ``your_username`` and ``your_password`` with your PyPI and PyPI testing usernames and passwords.

Once you have tagged the release of PyCBC a tarball should be available from the `GitHib PyCBC releases page <https://github.com/ligo-cbc/pycbc/releases>`_. Download this tarball, untar it, and change into the source directory of the tarball. 

First check the release with the PyPI test repository. Register the package with 

.. code:: bash

    python setup.py register -r pypitest

If you get no errors, you can then upload the package with

.. code:: bash

    python setup.py sdist upload -r pypitest

You should now see the package uploaded in the `PyPI Test Repository <https://testpypi.python.org/pypi>`_. If this is successful, you can puhlish it to the main repository with the commands

.. code:: bash

    python setup.py register -r pypi
    python setup.py sdist upload -r pypi

The package should then be available in PyPI for download.

============================================
Creating static binaries for production runs
============================================

Static binaries are entirely self-contained programs that

#. Can be run by anyone, whether or not they have set up a development environment
#. Ensure that each program runs with the right versions of all dependancies
#. can be distributed from a central location

A new set of such binaries should be generated for every release.

-------------
Preliminaries
-------------

To set up an environment for the release first follow the standard installation
instructions here:

.. toctree::
    :maxdepth: 1

    install

If the default Python is 2.6 it will be useful to name your virtual environment
something like ``pycbc_dev2.6``, since a subsequent step will require a 2.7
environment.

In order to distribute the binaries it will first be necessary to have a copy of the git
repository

.. code::

    git clone git@code.pycbc.phy.syr.edu:ligo-cbc/pycbc-software.git

Note that this will take a very long time.  Once the repository has been cloned
create a directory for the new release.

.. code:: bash

   mkdir -p v1.2.5/x86_64/composer_xe_2015.0.090

changing version number and architecture as appropriate.

A number of binaries change very rarely between releases and so can be copied
from a previous one:

.. code:: bash

    cp v1.2.4/x86_64/composer_xe_2015.0.090/* v1.2.5/x86_64/composer_xe_2015.0.090

Edit the ``README.md`` file for the new release.

Update the ``executables.ini`` file

.. code-block:: bash

    curl https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-config/download/master/O1/pipeline/executables.ini  -o executables.ini
    sed -i "s+\${which:\\(.*\\)}+http://code.pycbc.phy.syr.edu/pycbc-software/v1.2.5/x86_64/composer_xe_2015.0.090/\1+" executables.ini
    for exe in `grep http://code.pycbc executables.ini | awk '{print $1}'`
    do
      echo -e "[pegasus_profile-${exe}]\npycbc|installed = False\nhints|execution.site = local\n\n"
    done > profiles
    cat executables.ini profiles > tmp
    mv tmp executables.ini
    rm profiles

Again, change version number and architecture as needed.


---------------
lalapps_inspinj
---------------

This is the one C program from lalsuite that is still needed in the current
workflow.  This almost never changes and so copying from the previous release
will generally be fine.  If it does need to be rebuilt follow the instructions
for installing lalsuite but configure with

.. code::

    --enable-static-binaries --enable-static --disable-swig --disable-lalstochastic --disable-lalxml --disable-lalinference --disable-laldetchar --disable-lalpulsar --disable-framec

Then ``make`` and ``make install`` as usual.  The static exexecutable will be placed in the 
target ``bin`` directory.

----------------------
Other lalapps programs
----------------------

There are a few stochastic bank programs in lalsuite needed by pycbc.  These change infrequently and
can usually be copied from a previous release.  If they do need to be rebuilt the process is:

.. code-block:: bash

    cd /path/to/your/lalsuite
    cd lalapps/src/inspiral/

    for prog in \*sbank\*.py
    do
      pyinstaller ${prog}                          \
        --hidden-import scipy.linalg.cython_blas   \
        --hidden-import scipy.linalg.cython_lapack \
        --hidden-import scipy.special._ufuncs_cxx  \
        --hidden-import scipy.integrate            \
        --strip                                    \
        --onefile
    done

The resulting bundles will be placed in the ``dist`` directory.


----------------------
Segment database tools
----------------------

Client tools for the segment database change infrequently and can usually be
copied from the previous release.  If they do need to be rebuilt the process is:

.. code-block:: bash

    pyinstaller ligolw_segment_query_dqsegdb --strip --onefile
    pyinstaller ligolw_segments_from_cats_dqsegdb --strip --onefile


--------------
pyCBC binaries
--------------

The program used to create a static binary from a Python program is 
`PyInstaller <http://www.pyinstaller.org/>`_.  To set up PyInstaller:

.. code-block:: bash

    source /your/virtual/environment/bin/activate
    cd /your/virtual/environment/src
    git clone https://github.com/pyinstaller/pyinstaller.git
    cd pyinstaller
    git checkout 9d0e0ad4c1c02964bbff86edbf7400cd40958b1a

.. note::
   The fix to ``bootloader/common/pyi_utils.c`` only needs to be done when installing the python 2.7 virtual environment. 

By default programs built with PyInstaller will ignore the ``LD_LIBRARY_PATH``
environment variable, which causes problems in some environments.  To fix this
edit ``bootloader/common/pyi_utils.c`` and replace the function ``set_dynamic_library_path`` with the following

.. code-block:: c

    static int set_dynamic_library_path(const char* path)
    {
        int rc = 0;
    
    #ifdef AIX
        /* LIBPATH is used to look up dynamic libraries on AIX. */
        setenv("LIBPATH", path, 1);
        VS("%s\n", path);
    #else
        /* LD_LIBRARY_PATH is used on other *nix platforms (except Darwin). */
        char * curpath = getenv("LD_LIBRARY_PATH");
        if ( ! curpath ) { /* Use required path only */
            rc = setenv("LD_LIBRARY_PATH", path, 1);
            VS("%s\n", path);
        } else { /* Append current path onto required path */
            char apath[ strlen(path) + strlen(curpath) + 2 ];
            strcpy(apath, path);
            strcat(apath, ":");
            strcat(apath, curpath);
            rc = setenv("LD_LIBRARY_PATH", apath, 1);
            VS("%s\n", apath);
        }
    #endif /* AIX */
    
        return rc;
    }

.. Closing slash-star to keep vim happy /*

Then configure the bootloader and install as usual:

.. code-block:: bash

    cd bootloader
    ./waf configure build install --no-lsb
    cd ..
    python setup.py install


To ensure that pyCBC is set up properly prior to running pyinstaller, first
clean out the pip cache

.. code:: bash

    rm -rfv ${HOME}/.cache/pip

Then checkout the official release

.. code-block:: bash

    cd /your/virtual/environment/src
    git clone git@github.com:ligo-cbc/pycbc.git
    cd pycbc
    git checkout v1.2.5

replacing ``v1.2.5`` with the version to be built.  Then install pyCBC, using
the requirements file to ensure the correct version of all dependancies is
installed

.. code:: bash

  pip install -r requirements.txt


Then ensure that the installation went as expected:

.. code::

    pycbc_inspiral --version

    Branch: None Tag: v1.2.5 Id: 51dcf08cc6016a7574c3baf2efff2bb60ed6ce4f Builder:
    Larne Pekowsky <larne.pekowsky@ligo.org> Build date: 2015-10-31 14:48:20 +0000
    Repository status is CLEAN: All modifications committed


The tag should match the one that was checked out and the repository status should
report as ``CLEAN``.

To build static executables:

.. code-block:: bash

   cd tools/static
   bash build_dag.sh

This will construct a condor dag with a pyinstaller job for each binary.
Submit as usual:

.. code:: bash

   condor_submit_dag build_static.dag

Assuming everything goes well the resulting binaries will be placed in the
``dist`` directory.  As a final test, check the version again

.. code::

    dist/pycbc_inspiral --version

    Branch: None Tag: v1.2.5 Id: 51dcf08cc6016a7574c3baf2efff2bb60ed6ce4f Builder:
    Larne Pekowsky <larne.pekowsky@ligo.org> Build date: 2015-10-31 14:48:20 +0000
    Repository status is CLEAN: All modifications committed

In principle jobs could fail if pyinstaller fails to build the executable,
although this has never been seen in practice.  A job can also fail if
pyinstaller succeeds but the resulting program throws an error when invoked
with ``--help``.  Most of the time this happens it is because a new program has
been added and pyinstaller needs to be told that it needs scipy.  This is done
by adding the name of the new program to the ``needs_full_build`` file in the
``tools/static`` directory.  


--------------
pycbc_inspiral
--------------

This program needs to be built in a special environment so that it can be run on Open Science Grid sites.
All the necessary elements are available through CVMFS:

.. code-block:: bash

  source /cvmfs/oasis.opensciencegrid.org/osg/modules/lmod/current/init/bash
  module load python/2.7

Unfortunately it is now necessary to build an entire parallel develoment environment.  Move existing directories out of the way

.. code-block:: bash

    mv ${HOME}/local ${HOME}/local_2.6
    mv ${HOME}/.local ${HOME}/.local_2.6

and repeat the installation instructions starting from installing pip.  You
will also need to rebuild lalsuite.  For clarity you may want to call the
virtual environment something like ``pycbc_devel.2.7``. 

Once this second environment is set up in order to build:

.. code-block:: bash

    cd tools/static
    ./build_one.sh ../../bin/pycbc_inspiral 

and ensure that the resulting executable has the correct version

.. code::

    dist/pycbc_inspiral --version

    Branch: None Tag: v1.2.5 Id: 51dcf08cc6016a7574c3baf2efff2bb60ed6ce4f Builder:
    Larne Pekowsky <larne.pekowsky@ligo.org> Build date: 2015-10-31 14:48:20 +0000
    Repository status is CLEAN: All modifications committed


------------
Finishing up
------------

Once everything has been built and moved to the git repositry, commit and push as usual

.. code:: bash

    git commit -a -m "Version 1.2.5"
    git push

This will take some time.

Finally, the binaries will need to be put into place on the staging server, which must be done by someone with
root access.


