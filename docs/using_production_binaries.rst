.. _using_production_binaries:

#####################################
Using Production PyCBC Binary Bundles
#####################################

Production binaries for various releases of PyCBC are available built and
stored in a git repository from:

 * `<https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-software/tree/master>`_.

These binaries are built with `PyInstaller <http://www.pyinstaller.org/>`_ and
so provide completely self-contained execution environments. No local
installation of PyCBC is necessary to run the software.

===================
Workflow generators
===================

To run a workflow generator, download a copy of the binary bundle and make it
executable. For example, to download a copy of the program
``pycbc_make_coinc_search_workflow`` from the 1.2.0 release of PyCBC that has
been linked against the Intel Math Kernel Library, run the command::

    curl https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-software/raw/master/v1.2.0/x86_64/composer_xe_2015.0.090/pycbc_make_coinc_search_workflow > pycbc_make_coinc_search_workflow
    chmod +x pycbc_make_coinc_search_workflow

The executable can then be used in the normal way. For example, run::

    ./pycbc_make_coinc_search_workflow --help

to print the help message.

========================
Executables in workflows
========================

All of the PyCBC workflow generators are able to configure Pegasus to fetch
the programs from the remote server. To do this, simply pass the workflow
generator an on of the ``ini`` files from the pycbc-software repository. For
example using the file::

    https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-software/downlaod/master/v1.2.0/x86_64/composer_xe_2015.0.090/executables.ini

will configure the workflow to get and install pre-build binaries for the
1.2.0 release of PyCBC.
