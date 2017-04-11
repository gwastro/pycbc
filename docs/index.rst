.. PyCBC documentation master file, created by
   sphinx-quickstart on Tue Jun 11 17:02:52 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

###################
PyCBC documentation
###################

PyCBC is a software package used to explore astrophysical sources of gravitational waves. It contains algorithms that can detect signals using LIGO and measure the astrophysical parameters of detected sources. PyCBC was used in the first direct detection of gravitational waves (GW150914) and is used in the ongoing analysis of LIGO and Virgo data.

The easiest way to start using PyCBC is to install one of our `Docker containers <https://hub.docker.com/u/pycbc/>`_. First, install the `Docker Community Edition <https://www.docker.com/community-edition>`_ for your `Mac <https://store.docker.com/editions/community/docker-ce-desktop-mac?tab=description>`_ or `Windows <https://store.docker.com/editions/community/docker-ce-desktop-windows?tab=description>`_ desktop. Docker CE installations for `Linux platforms <https://www.docker.com/community-edition#/download>`_ are also available.

After installing and starting Docker, type the commands below to download and run a PyCBC container. This example downloads the ``v1.7.0`` release in a container. To get the current version of the code from the `GitHub master branch <https://github.com/ligo-cbc/pycbc>`_ replace the string ``v1.7.0`` with ``latest``. The container includes all of the required software and dependencies to run PyCBC, including a compatible version of LALSuite.

.. raw:: html

    <script type="text/javascript">
        document.addEventListener("DOMContentLoaded", function(){
            Typed.new(".element", {
                strings: ["^500docker pull pycbc/pycbc-el7:v1.7.0<br>$ ^500docker run -it pycbc/pycbc-el7:v1.7.0 /bin/bash -l<br>&#40;pycbc-software&#41;&#91;pycbc@37184573e664 &#126;&#93;$ ^500python<br>Python 2.7.5 &#40;default, Nov  6 2016, 00:28:07&#41;<br>&#91;GCC 4.8.5 20150623 &#40;Red Hat 4.8.5-11&#41;&#93; on linux2<br>&gt;&gt;&gt; ^500import pycbc.version<br>&gt;&gt;&gt; ^500print pycbc.version.git_tag<br>v1.7.0<br>&gt;&gt;&gt; ^500import lal.git_version<br>&gt;&gt;&gt; ^500print lal.git_version.id<br>539c8700af92eb6dd00e0e91b9dbaf5bae51f004<br>&gt;&gt;&gt; "],
                typeSpeed: 0
            });
        });
    </script>
    <div class="text-editor-wrap">
        <div class="title-bar"><span class="title">pycbc &mdash; bash &mdash; 80x<span class="terminal-height">25</span></span></div>
        <div class="text-body">
            $ <span class="element"></span>
        </div>
    </div>
    <br>
    <br>

Full installation instructions for users who want to install and develop PyCBC are available at:

.. toctree::
   :maxdepth: 1

   install

===================
PyCBC Project Goals
===================

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
A bibtex key and DOI for each release is avaliable from `Zenodo <http://zenodo.org/>`_.
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

Users who want to create and run scientific workflows to search for compact
binaries should read the documentation in the links at:

.. toctree::
   :maxdepth: 1

   workflow/pycbc_make_psd_estimation_workflow
   workflow/pycbc_make_coinc_search_workflow
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

Users who are intersted in using PyCBC utilities and functions should take a look at the
short code snippets below.

.. toctree::
   :maxdepth: 2

   gw150914
   frame
   psd
   noise
   waveform
   filter
   distributions

=============================
Documentation for Developers
=============================

Documentation on building stand-alone bundled executables with PyInstaller is available at:

.. toctree::
   :maxdepth: 1

   building_bundled_executables

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
