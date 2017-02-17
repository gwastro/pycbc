.. _using_pyinstaller_binaries::

######################################
Building and Using Bundled Executables
######################################

For most applications, PyCBC executables can be built and run from a standard
Python install following the instructions in :ref:`installing_pycbc`. The
standard installation requires that the install directory is available at
runtime so that Python can find the libraries. It also requires that the
build and runtime environments are compatible (e.g. compatible versions of
glibc, gcc, and Python). Some PyCBC executables (e.g. ``pycbc_inspiral``)
require runtime compilation of code using weave, and so the execution
environment must have a full installation of gcc and the Python development
libraries.

When running on, e.g. the Open Science Grid, these requirements may not be
satisfied. Although CVMFS can provide access to the PyCBC libraries at
runtime on the OSG, many OSG execute machines do not have the required
environment to weave-compile code at runtime. Running PyCBC using
Einstein@Home places even stricter limitations on the runtime environment of
the code, as execute machines do not have access to CVMFS and may be running a
variety of operating systems.

For both OSG and Einstein@Home, a **bundled** executable must be built using
PyInstaller. This bundle contains all of the Python and user-space C libraries
required, as well as a Python interperter to run the code. This bundle must
also contain pre-compiled objects for all of the weave code that is needed at
runtime. The script ``pycbc_build_eah.sh`` can be used to build a
self-contained PyInstaller bundle for ``pycbc_inspiral`` so that it can be run
on OSG and Einstein@Home.

Note that the PyInstaller bundles are not completely static and still require
a dynamically linked version of glibc at runtime. Since Linux systems are
backwards, but not forwards compatible the bundle must be built on the
lowest-common denominator operating system for the execution platform. For OSG
this is RHEL6, as the OSG executale machines are typically RHEL6 or RHEL7.
RHEL6 executables also work on Debian Wheezy and Jessie. Einstein@Home
requires an even older build platforms, typically Debian Etch.  Since these
older build platforms may not have the required software installed (e.g. FFTW,
Python, etc.) the ``pycbc_build_eah.sh`` downloads and builds a complete
installation environment of **all** of the required software. It can therefore
be used to build PyCBC on a machine without the standard LIGO Data Grid
software installed.

======================
Using the build script
======================

The command-line arguments for the ``pycbc_build_eah.sh`` build script are:

.. command-output:: pycbc_build_eah.sh --force-debian4 --help

.. note::

    The build script checks that it is being run on one of the lowest-common
    denominator platforms that it knows about. To run the script on another
    platform, pass the command-line argument ``--force-debian4`` to the script
    as the **first** argument in the list of arguments.

The minimal set of command line options require to ``pycbc_inspiral`` is
typically the hash of the versions of LALSuite and PyCBC required::

    pycbc_build_eah.sh --lalsuite-commit=a2a5a476d33f169b8749e2840c306a48df63c936 --pycbc-commit=b68832784969a47fe2658abffb3888ee06cd1be4

To include extra run-time libraries in the bundle, e.g. to add the Intel MKL
libraries, specify them with the command-line argument::

    --with-extra-libs=file://opt/intel/composer_xe_2015.0.090.tar.gz

The script executes ``pycbc_inspiral`` as part of the build process. This may
require LAL data at build time. The LAL data can be given with the command
line argument::
    
    --with-lal-data=/cvmfs/oasis.opensciencegrid.org/ligo/sw/pycbc/lalsuite-extra/11/share/lalsimulation

The default command line arguments clone PyCBC from the standard GitHub
repository.  If you would like to build a bundle using code from your own
GitHub repository ior branch you can use the arguments::

    --pycbc-remote=soumide1102 --pycbc-branch=comp_wave_in_search

You may also tell the script to run ``pycbc_inspiral`` with additional
template banks or approximants to ensure that all of the necessary weave code
is compiled into the executable with the arguments::

    --with-extra-approx='SPAtmplt:mtotal<4'
    --with-extra-approx='SEOBNRv4_ROM:else'
    --with-extra-approx=--use-compressed-waveforms
    --with-extra-bank=/home/soumi.de/projects/cbc/SEOBNRROM-proj/testbank_TF2v4ROM.hdf


===========================
Building Releases for CVMFS
===========================

To build a release of ``pycbc_inspiral`` for installation in CVMFS, run the
script with the arguments::

    ./pycbc_build_eah.sh --lalsuite-commit=a2a5a476d33f169b8749e2840c306a48df63c936 --pycbc-commit=b68832784969a47fe2658abffb3888ee06cd1be4 --with-extra-libs=file://`pwd`/composer_xe_2015.0.090.tar.gz

changing the ``--lalsuite-commit`` and ``--pycbc-commit`` to the appropriate
hashes for the release.
