.. PyCBC documentation master file, created by
   sphinx-quickstart on Tue Jun 11 17:02:52 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

###################
PyCBC documentation
###################

PyCBC is a python toolkit for analysis of data from gravitational-wave laser interferometer detectors with the goal of detecting and studying signals from compact binary coalescences (CBCs).

The goals of the PyCBC project are to

- Create a flexible, extensible sw production code for CBC analysis that can be released for the public
- Enable simple, easy and transparent access for various many-core architectures like GPUs
- Ultimately become the data analysis tool of the 'advanced era'

===========================
PyCBC Contributors
===========================

Stanislav Babak (1), Christopher M. Biwer (2), Duncan Brown (2), Collin Capano (3), Tito Dal Canton (4), Gergely Debreczeni (5), Thomas Dent (4), Steve Fairhurst (10), Henning Fehrmann (4), Ian Harry (2), Drew Keppel (4), Badri Krishnan (4), Prayush Kumar (2), Andrew Lundgren (4), Duncan Macleod (6), Adam Mercer (7), Andrew Miller (8), Saeed Mirshekari (9), Alex Nitz (2), Karsten Wiesner (4), Josh Willis (8)

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

=========================
Links
=========================

Getting Started:

.. toctree::
   :maxdepth: 1

   install

Executable/Package documentation:

.. toctree::
   :maxdepth: 1

   workflow
   tmpltbank
   frame
   banksim
   faithsim

Walkthroughs and Tutorials:

.. toctree::
   :glob:
   :maxdepth: 1

   workflow/*walkthrough*

Code Examples:

.. toctree::
   :maxdepth: 2

   waveform

For Developers:

.. toctree::
    :maxdepth: 1
    
    documentation
    release
    
Format Specifications:

.. toctree::
   :glob:
   :maxdepth: 1

   formats/*

===================
Indexes and Tables
===================

* :ref:`modindex`
* :ref:`genindex`
