################
Installing PyCBC
################

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

============
Dependencies
============

PyCBC is heavily dependent on the `LIGO Algorithm Library (LAL) <https://www.lsc-group.phys.uwm.edu/daswg/projects/lalsuite.html>`_.
