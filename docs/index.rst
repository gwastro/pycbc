=======
PyCBC
=======

PyCBC is a software package used to explore astrophysical sources of gravitational waves.
It contains algorithms that can detect coalescing compact binaries and measure
the astrophysical parameters of detected sources. PyCBC was used
in the `first direct detection of gravitational waves (GW150914) by
LIGO <https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.116.061102>`_ and
is used in the ongoing analysis of LIGO and Virgo data.
If you use PyCBC in your scientific publications or projects,
we ask that you acknowlege our work by citing the papers described on the page:

.. toctree::
   :maxdepth: 1

   credit

===============
Getting Started
===============

 - Use the PyCBC Library within your Browser

   `Try out our tutorials <https://github.com/gwastro/PyCBC-Tutorials>`_.

=====================
Installation
=====================

Note, if you are a LIGO / Virgo member with access to LDG resources, PyCBC is *already*
installed on your cluster through CVMFS! Instructions to source any release of PyCBC
is available from the `releases page <https://github.com/gwastro/pycbc/releases>`_.

You may also install PyCBC directly with pip. You may ommit `lalsuite` if you have
your own build.

.. code-block:: bash

   pip install lalsuite pycbc

Full detailed installation instructions which covers other installation cases:

.. toctree::
   :maxdepth: 1

   install
   install_virtualenv
   install_lalsuite

==================================================
Parameter Estimation of Gravitational-wave Sources
==================================================

Users who want to create and run parameter estimation workflows should read the
documentation at:

.. toctree::
   :maxdepth: 2

   inference

========================================
Searching for Gravitational-wave Signals
========================================

Users who want to create and run scientific workflows to search for compact
binaries should read the documentation in the links at:

.. toctree::
   :maxdepth: 1

   workflow/pycbc_make_psd_estimation_workflow
   workflow/pycbc_make_coinc_search_workflow
   workflow/pygrb.rst

==================
Indexes and Tables
==================

* :ref:`modindex`
* :ref:`genindex`
