=============================================
PyCBC: Powering Gravitational-wave Astronomy
=============================================

PyCBC is the result of a community effort to build
a set of core libraries and application suites
used to study gravitational-wave data and astrophysics.
It contains algorithms that can detect coalescing compact binaries and measure
the astrophysical parameters of detected sources. PyCBC was used
in the `first direct detection of gravitational waves (GW150914) by
LIGO <https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.116.061102>`_ and
is used in the ongoing analysis of LIGO and Virgo data.

If you are interesting in building community tools for
gravitational-wave astronomy, please consider `contributing <https://github.com/gwastro/pycbc>`_, whether it is
providing feedback, examples, documentation or helping to improve the core
library and application suite.

.. carousel::
    :show_controls:
    :show_indicators:
    :show_dark:
    :show_shadows:
    :show_captions_below:

    .. figure:: https://pycbc.org/pycbc/latest/html/_images/data.png
        :height: 400px
        :target: catalog.html

        Working with gravitational wave data

    .. figure:: https://pycbc.org/pycbc/latest/html/_images/plot_detwaveform.png
        :height: 400px
        :target: waveform.html

        Your interface to generating gravitational wave signals

    .. figure:: https://collincapano.com/wp-content/uploads/2020/02/posterior3d-gw150914_masses-e1581860222179.png
        :height: 400px
        :target: inference.html

        Flexible, easy-to-use, parameter estimation for GW Astronomy

    .. figure:: https://pycbc.org/pycbc/latest/html/_images/demarg_150914.png
        :height: 400px
        :target: inference/examples/margtime.html
        
        Fast Parameter Estimation with Advanced Marginalizations
        

.. toctree::
   :hidden:
   :maxdepth: 1

   install
   credit

   Modules <modules>
   genindex

.. toctree::
   :hidden:
   :caption: User Guides
   :maxdepth: 1

   tutorials
   search
   inference
   apps

.. toctree::
   :caption: Dev Guides
   :hidden:
   :maxdepth: 1

   extend
   devs

.. toctree::
   :hidden:
   :maxdepth: 1

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
    +++
    Documentation may be out of data. Please help us improve it!

.. card:: Detect Signals in association with GRBs
    :link:  pygrb
    :link-type: ref

    Targeted analysis to detect gravitational-wave sources in association
    with gamma-ray bursts and other transient sources.
    +++
    Expected for full implementation in PyCBC in the next year.

.. card:: Other applications
    :link:  apps
    :link-type: ref

    Documentation for a select sample of the pycbc software suite. These
    include for generating template banks, hardware injections, etc.


================
Getting Started
================

 - Use the PyCBC Library within your Browser

   `Try out our tutorials <https://github.com/gwastro/PyCBC-Tutorials>`_.

=====================
Installation
=====================

You may also install PyCBC directly with pip or conda.

.. code-block:: bash

   pip install pycbc

Detailed instructions are found :ref:`here <installing_pycbc>`.

Note, if you are a LIGO / Virgo member with access to IGWN resources, PyCBC is *already*
installed on your cluster through CVMFS! Instructions to source any release of PyCBC
is available from the `releases page <https://github.com/gwastro/pycbc/releases>`_.
