.. PyCBC documentation master file, created by
   sphinx-quickstart on Tue Jun 11 17:02:52 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

###################
PyCBC documentation
###################

PyCBC is a python toolkit for analysis of data from gravitational-wave laser interferometer detectors with the goal of detecting and studying signals from compact binary coalescences (CBCs).

The goals of the PyCBC project are to:

- Provide tools for building gravitational-wave search workflows for CBCs
- Create a flexible, extensible production code for CBC analysis that can be released for the public
- Enable simple, easy and transparent access for various many-core architectures like GPUs

===========================
PyCBC Contributors
===========================

Stanislav Babak (1), Christopher M. Biwer (2), Duncan Brown (2), Collin Capano (3), Tito Dal Canton (4), Gergely Debreczeni (5), Thomas Dent (4), Steve Fairhurst (10), Henning Fehrmann (4), Ian Harry (1,2), Marcel Kehl (11), Drew Keppel (4), Badri Krishnan (4), Prayush Kumar (2,11), Amber Lenon (2), Andrew Lundgren (4), Duncan Macleod (6), Thomas Massinger (2), Adam Mercer (7), Andrew Miller (8), Saeed Mirshekari (9), Alex Nitz (1,2), Laura Nuttall (2), Francesco Pannarale (10), Harald Pfeiffer (11), Samantha Usman (2), Karsten Wiesner (4), Andrew Williamson (10), Josh Willis (8).

#. AEI, Golm, Germany
#. Syracuse University, NY, USA
#. University of Maryland, MD, USA
#. AEI, Hannover, Germany
#. Wigner RCP, Budapest, Hungary
#. Louisiana State University, LA, USA
#. University of Wisconsin-Milwaukee, WI, USA
#. Abilene Christian University, TX, USA
#. Instituto de Fisico Teorica, Sao Paulo, Brazil
#. Cardiff University, Cardiff, UK
#. Canadian Institute for Theoretical Astrophysics, Toronto, Canada

=========================
Documentation for Users
=========================

To install PyCBC and get started, follow the links at:

.. toctree::
   :maxdepth: 1

   install

Users who want to create and run scientific workflows to search for compact
binaries should read the documentation in the links at:

.. toctree::
   :maxdepth: 1

   workflow/pycbc_make_psd_estimation_workflow
   workflow/pycbc_make_coinc_search_workflow
   workflow/pycbc_make_sngl_workflow

Users who are interested in tools that PyCBC provides for various other
analysis tasks (e.g. template bank generation, hardware injections, and testing
template banks) should read the documentation at:

.. toctree::
   :maxdepth: 1

   tmpltbank
   hwinj
   uberbank_verify
   banksim
   faithsim

Users who are intersted in using PyCBC for investigation of CBC waveforms
should read the documentation at:

.. toctree::
   :maxdepth: 2

   waveform

=============================
Documentation for Developers
=============================

PyCBC developers should read the pages below which explain how to write
documentation, develop the code, and create releases:

.. toctree::
   :maxdepth: 1
    
   documentation
   release

Developers who are interested in file I/O, data storage, and access should
read the documentation at:

.. toctree::
   :maxdepth: 1

   frame
   formats/hdf_format

Developers who are interested in creating new scientific workflow generation
scripts should read the documentation at:

.. toctree::
   :maxdepth: 1

   workflow

Full Module Documentation is avaialable at:

.. toctree::
   :maxdepth: 1

   modules

===================
Indexes and Tables
===================

* :ref:`modindex`
* :ref:`genindex`
