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
   :hidden:
   :maxdepth: 1

   inference
   apps
   tutorials
   extend
   devs
   credit
   modules
   genindex
   modindex

.. card:: Core Library Examples
    :link:  tutorials
    :link-type: ref

    The core library provides the tools to do GW data analysis.
    See examples of how to use PyCBC core tools to read gravitational-wave
    data, query detector status, filter data, generate waveforms, estimate
    PSDs, calculate SNRs, and much more.


.. card:: PyCBC Inference
    :link:  inference
    :link-type: ref

    Easy-to-use, configurable, and robust Parameter Estimation for Gravitational-wave Astronomy and Multi-messenger analysis.


.. card:: Detect Signals in Archival Data
    :link:  search_workflow
    :link-type: ref

    The deep analysis used to detect gravitational-wave sources in
    archival data.

.. card:: Detect Signals in association with GRBs
    :link:  pygrb
    :link-type: ref

    Targetted analyssi to detect gravitational-wave sources in association
    with gamma-raty bursts and other transient sources.

================
Getting Started
================

 - Use the PyCBC Library within your Browser

   `Try out our tutorials <https://github.com/gwastro/PyCBC-Tutorials>`_.

=====================
Installation
=====================

Note, if you are a LIGO / Virgo member with access to LDG resources, PyCBC is *already*
installed on your cluster through CVMFS! Instructions to source any release of PyCBC
is available from the `releases page <https://github.com/gwastro/pycbc/releases>`_.

You may also install PyCBC directly with pip.

.. code-block:: bash

   pip install pycbc

Full detailed installation instructions which covers other installation cases:

.. toctree::
   :maxdepth: 1

   install

==================
Indexes and Tables
==================

* :ref:`modindex`
* :ref:`genindex`
