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

=======================================
Use of PyCBC in Scientific Publications
=======================================

If you use any code from PyCBC in a scientific publication, then we ask that
it is cited in the following way:

::

    These results were generating using the PyCBC software package
    \cite{Canton:2014ena,Usman:2015kfa,pycbc-software}

For the citation ``pycbc-software``,  please use a bibtex entry and DOI for the
appropriate release of the PyCBC software (or the latest available release).
A bibtex key and DOI for each release is avaliable from `Zendo <http://zendo.org/>`_.
A key for the latest release is available at:

.. image:: https://zenodo.org/badge/31596861.svg
   :target: https://zenodo.org/badge/latestdoi/31596861

Bibtex keys for the citations ``Canton:2014ena`` and ``Usman:2015kfa`` are::

    @article{Canton:2014ena,
      author         = "Dal Canton, Tito and others",
      title          = "{Implementing a search for aligned-spin neutron
                        star-black hole systems with advanced ground based
                        gravitational wave detectors}",
      journal        = "Phys. Rev.",
      volume         = "D90",
      year           = "2014",
      number         = "8",
      pages          = "082004",
      doi            = "10.1103/PhysRevD.90.082004",
      eprint         = "1405.6731",
      archivePrefix  = "arXiv",
      primaryClass   = "gr-qc",
      reportNumber   = "LIGO-P1400053",
      SLACcitation   = "%%CITATION = ARXIV:1405.6731;%%"
    }

    @article{Usman:2015kfa,
      author         = "Usman, Samantha A. and others",
      title          = "{The PyCBC search for gravitational waves from compact
                        binary coalescence}",
      journal        = "Class. Quant. Grav.",
      volume         = "33",
      year           = "2016",
      number         = "21",
      pages          = "215004",
      doi            = "10.1088/0264-9381/33/21/215004",
      eprint         = "1508.02357",
      archivePrefix  = "arXiv",
      primaryClass   = "gr-qc",
      reportNumber   = "LIGO-P1500086",
      SLACcitation   = "%%CITATION = ARXIV:1508.02357;%%"
    }


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
   workflow/pygrb.rst

Users who want to create and run parameter estimation workflows should read the
documentation at:

.. toctree::
      :maxdepth: 1

   workflow/pycbc_make_inference_workflow

Users who are interested in tools that PyCBC provides for various other
analysis tasks (e.g. template bank generation, hardware injections, and testing
template banks) should read the documentation at:

.. toctree::
   :maxdepth: 1

   tmpltbank
   inference
   hwinj
   banksim
   faithsim
   upload_to_gracedb

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
