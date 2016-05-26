##########################
Creating Releases of PyCBC
##########################

.. note::

    Only authorized maintainers of PyCBC should use these instructions.


==============================
Creating the release on GitHub
==============================

To create a new PyCBC release:

#. Make sure that the setup.py file contains the correct version number, which should be in the format ``x.y.z`` (where x, y, and z are the major, minor, and patch levels of the release) in PyCBC's setup.py file.
#. Set ``release = True`` in the PyCBC's setup.py file.
#. Commit the changed setup.py file and push to commit to the repository.
#. Go to the `PyCBC release page <https://github.com/ligo-cbc/pycbc/releases>`_ and click on ``Draft a new release``. 
#. Enter a tag version in the format ``vx.y.z``. Note the ``v`` in front of the major, minor, and patch numbers. 
#. Enter a descriptive release title and write a desciption of the release in the text box provided.
#. Click on ``Publish release`` to create the release.
#. Update the setup.py file with an incremented major or minor version number and make sure that the string ``dev`` appears in that version. For example, if you just released ``1.2.1`` then change the string to ``1.3.dev0`` or ``2.0.dev0`` as appropriate. This is needed to ensure that if someone is building from source, it always takes precedence over an older release version.
#. Set ``release = False`` in PyCBC's setup.py file.
#. Commit these changes and push them to the repository.

.. note::

    Releases of PyCBC should be made from the master. Do not make a branch
    unless you are back-porting a bug fix from a new release series to an
    old production release series.

------------------------------------------------
Backporting Bug Fixes to Previous Release Series
------------------------------------------------

Branches should only be created when bug fixes from master need to be back
ported to an old release series (e.g. adding a bug fix from the 1.4 series to
the 1.3 series. 

To create a branch for the bug fix, make a new branch from the last release
tag and cherry pick the changes to that branch. For example, to create a
``1.3.8`` branch to add bug fixes to the ``1.3.7`` release, you should

.. code:: bash

  git checkout -b b1.3.8 v1.3.7
  git cherry-pick [hash of first change]
  git cherry-pick [hash of second change]
  [continue using git cherry-pick for all new changes]
  git commit
  git push

Then go to the `PyCBC release page <https://github.com/ligo-cbc/pycbc/releases>`_ and click on ``Draft a new release`` following the instructions above.

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

You should now see the package uploaded in the `PyPI Test Repository <https://testpypi.python.org/pypi>`_. If this is successful, you can publish it to the main repository with the commands

.. code:: bash

    python setup.py register -r pypi
    python setup.py sdist upload -r pypi

The package should then be available in PyPI for download.

=============
Preliminaries
=============

If ``bundle_env2.6.tgz`` and ``bundle_env2.7.tgz`` exist then you do not need to follow these steps as they help set up the virtual environments.

.. code:: bash
            
  cd ~
  tar -xf bundle_env2.6.tgz
  tar -xf bundle_env2.7.tgz

.. note::

   Every time you log into ``sugar-dev2`` you will need to run the commands ``source /cvmfs/oasis.opensciencegrid.org/osg/modules/lmod/current/init/bash`` and ``module load python/2.7``.

You will need to have a copy of lalsuite to distribute the binaries.

.. code:: bash

   mkdir -p /path/to/home/src
   cd /path/to/home/src
   git clone git://versions.ligo.org/lalsuite.git

In order to distribute the binaries it will first be necessary to have a copy of the git repository

.. code::
    
    cd ${HOME}
    git clone git@code.pycbc.phy.syr.edu:ligo-cbc/pycbc-software.git

To set up an environment for the release first follow the standard installation
instructions here or use the install script to create ``bundle_env2.6``: 

.. toctree::
       :maxdepth: 1

    install

If the default Python is 2.6 it will be useful to name your virtual environment
something like ``bundle_env2.6``, since a subsequent step will require a 2.7
environment.

.. note::
    
  The bundles can still pull files from the virtual environment that you create. To prevent this it is necessary to ensure that the build environment does not persist after you create the bundles.

A ``bundle_env2.6``  environment needs to be created on ``sugar-dev3``. Once you have created an environment on ``sugar-dev3`` another environment ``bundle_env2.7`` needs to be created on ``sugar-dev2``. The ``pycbc_inspiral`` bundle will need to be generated in the 2.7 environment.

.. code:: bash

    cd ~/
    ./build_bundle_env.sh 

.. note::

      The build_bundle_env.sh can be found at https://github.com/ligo-cbc/pycbc.

Before creating your 2.7 virtual environment:

.. code:: bash

    mv ${HOME}/local ${HOME}/local_2.6
    mv ${HOME}/.local ${HOME}/.local_2.6

Once your 2.7 virtual environment is built, run the same move commands and change 2.6 to 2.7. These steps will configure your virtual environments to be able to build bundles.

Install pyCBC using the requirements file to ensure the correct version of all dependancies is installed. The steps in the following box has to be done on both ``sugar-dev3`` and ``sugar-dev2`` with the correct bundle environment name in both cases. 

.. code:: bash

    source ${HOME}/bundle_env2.6/bin/activate
    cd ${HOME}/bundle_env2.6/src
    rm -rf pycbc
    git clone git@github.com:[github username]/pycbc.git
    cd pycbc
    git remote add upstream git@github.com:ligo-cbc/pycbc.git
    git remote -vvv

Running ``git remote -vvv`` will allow you to check that you origin points to your version of pycbc.

The program used to create a static binary from a Python program is
`PyInstaller <http://www.pyinstaller.org/>`_.  To set up PyInstaller do the following in both ``sugar-dev3`` and ``sugar-dev2`` :

.. code:: bash

    cd $VIRTUAL_ENV/src
    git clone https://github.com/pyinstaller/pyinstaller.git
    cd pyinstaller
    git checkout 9d0e0ad4c1c02964bbff86edbf7400cd40958b1a
    vim bootloader/common/pyi_utils.c

By default programs built with PyInstaller will ignore the ``LD_LIBRARY_PATH``
environment variable, which causes problems in some environments.  To fix this
edit ``bootloader/common/pyi_utils.c`` and replace the function ``set_dynamic_library_path`` with the following

.. note::

    The fix to ``bootloader/common/pyi_utils.c`` only needs to be done when installing the python 2.7 virtual environment on ``sugar-dev2``. 

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

On ``sugar-dev3`` since you do not need to make the change to bootloader, you just need to follow this step. However this will also need to be done after you make the change to bootloader on ``sugar-dev2``.
    
.. code-block:: bash   

    python setup.py install
    which pyinstaller

Running ``which pyinstaller`` prints something like: ``~/bundle_env2.6/bin/pyinstaller``. If you see this then it means that pyinstaller is successfully installed.


====================
Updates to Lalsuite
====================

This will need to be done every time there is a change to lalsuite.

.. code:: bash

   rm -rf ${HOME}/lalsuite
   mkdir -p ${HOME}/src
   cd ${HOME}/src
   git clone git://versions.ligo.org/lalsuite.git
   cd lalsuite
   git checkout [lalsuite_tag]
   git cherry-pick [commit]

.. note::
      Always verify the correct lalsuite tag or branch before updating lalsuite.

============================================
Creating static binaries for production runs
============================================

Static binaries are entirely self-contained programs that

#. Can be run by anyone, whether or not they have set up a development environment
#. Ensure that each program runs with the right versions of all dependancies
#. can be distributed from a central location

A new set of such binaries should be generated for every release.

--------------------
Bundle Preliminaries
--------------------

To update your virtual environments with the one on github:

Log in to ``sugar-dev3`` and continue with the following steps :

.. code::

   cd $VIRTUAL_ENV/src/pycbc
   git fetch upstream
   git rebase upstream/master
   git push
   git checkout [new version]
   python setup.py install

Create a directory for the new release.

.. code:: bash

   mkdir -p ${HOME}/pycbc-software/v1.2.5/x86_64/composer_xe_2015.0.090

changing version number and architecture as appropriate.

A number of binaries change very rarely between releases and so can be copied
from a previous one:

.. code:: bash

    cp ${HOME}/pycbc-software/v1.2.4/x86_64/composer_xe_2015.0.090/* ${HOME}/pycbc-software/v1.2.5/x86_64/composer_xe_2015.0.090

Edit the ``README.md`` file and update the ``executables.ini`` file for the new release.

.. note::
      The ``README.md`` file and the ``executables.ini`` file live in ``${HOME}/pycbc-software/v1.2.5/x86_64/composer_xe_2015.0.090/`` with the version being the one you are building.

To find the commit to add to the ``README.md``:
 
.. code:: bash

   git log

The commit you should add will look like:

.. code::

  commit 3a98bf98cb95e363b218a23734e7d0f32a58807f
  Author: Larne Pekowsky larne.pekowsky@ligo.org
  Date:   Fri Mar 25 20:29:05 2016 -0400

       v1.3.7

To update the ``executables.ini`` file do the following :

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

An executables.ini will also need to be created for running on OSG.

.. code-block:: bash

   cp executables.ini osg_executables.ini
   vim osg_executables.ini

Once you are in the ``osg_executables.ini`` file. You need to make some changes to the locations of the executables.
In vim after pressing ``shift+;``, run the command 

.. code::

   %s#http://code.pycbc.phy.syr.edu#/opt#gc

After this change you will have to manually change inpsiral back to: ``inspiral = http://code.pycbc.phy.syr.edu/pycbc-software/v1.2.5/x86_64/composer_xe_2015.0.090/pycbc_inspiral`` with the correct version number.

.. note::

   lalapps_inspinj, Other lalapps programs, Segment database tools, and pyCBC binaries only needs to be built on ``sugar-dev3``.

---------------
lalapps_inspinj
---------------

This is the one C program from lalsuite that is still needed in the current
workflow.  This almost never changes and so copying from the previous release
will generally be fine.  If it does need to be rebuilt follow the instructions
for installing lalsuite but configure with

.. code::

    --enable-static-binaries --enable-static --disable-swig --disable-lalstochastic --disable-lalxml --disable-lalinference --disable-laldetchar --disable-lalpulsar --disable-framec

Then ``make`` and ``make install`` as usual.  The static executable will be placed in the 
target ``bin`` directory.

----------------------
Other lalapps programs
----------------------

There are a few stochastic bank programs in lalsuite needed by pycbc.  

.. note::

    These change infrequently and can usually be copied from a previous release.  If they do need to be rebuilt the process is:

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

Then ensure that the installation went as expected:

.. code::

    pycbc_inspiral --version

    Branch: None Tag: v1.2.5 Id: 51dcf08cc6016a7574c3baf2efff2bb60ed6ce4f Builder:
    Larne Pekowsky <larne.pekowsky@ligo.org> Build date: 2015-10-31 14:48:20 +0000
    Repository status is CLEAN: All modifications committed


The tag should match the one that was checked out and the repository status should report as ``CLEAN``.

To build static executables:

.. code-block:: bash

   cd tools/static
   bash build_dag.sh

This will construct a condor dag with a pyinstaller job for each binary.
Submit as usual:

.. code:: bash

   condor_submit_dag build_static.dag

If you want to keep track of the build process run:

.. code::

   tail -f build_static.dag.dagman.out

Assuming everything goes well the resulting binaries will be placed in the
``dist`` directory.  As a final test, check the version again

.. code::

    dist/pycbc_inspiral --version

    Branch: None Tag: v1.2.5 Id: 51dcf08cc6016a7574c3baf2efff2bb60ed6ce4f Builder:
    Larne Pekowsky <larne.pekowsky@ligo.org> Build date: 2015-10-31 14:48:20 +0000
    Repository status is CLEAN: All modifications committed


In principle jobs could fail if pyinstaller fails to build the executable, although this has never been seen in practice. A job can also fail if pyinstaller succeeds but the resulting program throws an error when invoked with --help. There are generally two ways in which a build can fail; it may need the full set of scipy hooks, or for some reason pyinstaller is just unable to bundle it. The first case might happen if a new program has been added and pyinstaller needs to be told that it needs scipy, which can be resolved by adding the name of the executable to the ``needs_full_build`` file in the ``tools/static`` directory and running the condor dag again. The second case can be handled by adding the name of the executable to the ``cant_be_built`` file in the ``tools/static`` directory, provided the program is not needed by the workflow and running the condor dag again. The executables.ini lists all the files needed by the workflow. Therefore, anything not there can be skipped.

When the build finishes you will have to move the bundles to the directory under ``pycbc-software``.

.. code:: bash

   mv dist/* ${HOME}/pycbc-software/v1.2.5/x86_64/composer_xe_2015.0.090

Making the change to the correct version.

To clean up the build directory on ``sugar-dev3`` after you move the bundles you will need to be in ``tools/static`` directory:

.. code:: bash

   rm -rf *err *out *spec dist/* build/* build_static.dag*

--------------
pycbc_inspiral
--------------

.. note::

   ``pycbc_inspiral`` can only be built on ``sugar-dev2``.

This program needs to be built in a special environment so that it can be run on Open Science Grid sites.
All the necessary elements are available through CVMFS:

.. code-block:: bash

  source /cvmfs/oasis.opensciencegrid.org/osg/modules/lmod/current/init/bash
  module load python/2.7

Activate your virtual environment

.. code:: bash

   source ${HOME}/bundle_env2.7/bin/activate
   cd ${HOME}/bundle_env2.7/src/pycbc

To update your virtual environment with the one on github.

.. code:: bash

   git fetch upstream
   git rebase upstream/master
   git push
   git checkout [new version]
   python setup.py install

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

When the build finishes you will have to move ``pycbc_inspiral`` to the directory under ``pycbc-software``.

.. code:: bash

      mv dist/pycbc_inspiral ${HOME}/pycbc-software/v1.2.5/x86_64/composer_xe_2015.0.090

To clean up the build directory on ``sugar-dev2`` after you move the bundles you will need to be in ``tools/static`` directory:

.. code:: bash 

   rm -rf *spec dist/* build/*

------------
Finishing up
------------

Once everything has been built and moved to ``${HOME}/pycbc-software``, commit and push as usual to the git binaries repository doing the following steps, changing the version to the one you have built. This can be done either from ``sugar-dev2`` or from ``sugar-dev3``. 

.. code:: bash
 
    cd ${HOME}/pycbc-software
    git add v1.2.5
    git commit -a -m "Version 1.2.5"
    git push

This will take some time.

Finally, the binaries will need to be put into place on the staging server, which must be done by someone with root access.

--------------------------
Saving Virtual Environment
--------------------------

To save your virtual environment as a .tar file

.. code:: bash
  
  cd ~
  tar -zcf bundle_env2.6.tgz bundle_env2.6
  tar -zcf bundle_env2.7.tgz bundle_env2.7
  rm -r bundle_env2.6 bundle_env2.7


Before creating the next release

.. code:: bash
   
  cd ~
  tar -xf bundle_env2.6.tgz
  tar -xf bundle_env2.7.tgz


