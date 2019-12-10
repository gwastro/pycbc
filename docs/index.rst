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

 -  Get the full PyCBC software suite with one of our Docker Images

    If you have `Docker <https://www.docker.com/community-edition>`_ installed, you can get started using PyCBC with just two commands:

    .. raw:: html

        <script type="text/javascript">
            document.addEventListener("DOMContentLoaded", function(){
                Typed.new(".element", {
                    strings: ["^500<strong>docker pull pycbc/pycbc-el7:latest</strong><br>$ ^500<strong>docker run -it pycbc/pycbc-el7:latest</strong><br>&#40;pycbc-software&#41;&#91;pycbc@37184573e664 &#126;&#93;$ ^500<strong>python</strong><br>Python 2.7.5 &#40;default, Nov  6 2016, 00:28:07&#41;<br>&#91;GCC 4.8.5 20150623 &#40;Red Hat 4.8.5-11&#41;&#93; on linux2<br>&gt;&gt;&gt; ^500<strong>execfile&#40;&quot;/opt/pycbc/src/pycbc/examples/waveform/match_waveform.py&quot;&#41;</strong><br>^1000The match is: 0.953<br>&gt;&gt;&gt; ^500<strong>from pycbc.waveform import td_approximants</strong><br>&gt;&gt;&gt; ^500<strong>print td_approximants&#40;&#41;&#91;20:24&#93;</strong><br>['SEOBNRv3', 'SEOBNRv2', 'SpinTaylorT1', 'SEOBNRv4']<br>&gt;&gt;&gt; "],
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

    For more details, including instructions on starting a container that can display graphics, see:

    .. toctree::
       :maxdepth: 1

       docker


 - Use the PyCBC Library within your Browser

   We have interactive jupyter notebooks on azure where you can try some of our
   examples, and use the pycbc library. All from within the browser.
   `Give it a try <https://notebooks.azure.com/nitz/libraries/pycbc>`_

=====================
Installation
=====================

Note, if you are a LIGO / Virgo member with access to LDG resources, PyCBC is *already*
installed on your cluster through CVMFS! Instructions to source any release of PyCBC
is available from the `releases page <https://github.com/ligo-cbc/pycbc/releases>`_.

You may also install PyCBC directly with pip. You may ommit `lalsuite` if you have
your own build.

.. code-block:: bash

   pip install lalsuite pycbc

Full detailed installation instructions for users who want to use and develop PyCBC are available at:

.. toctree::
   :maxdepth: 1

   install

=======================
Documentation for Users
=======================

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
   :maxdepth: 2

   inference

Users who are interested in tools that PyCBC provides for various other
analysis tasks (e.g. template bank generation, hardware injections, and testing
template banks) should read the documentation at:

.. toctree::
   :maxdepth: 1

   tmpltbank
   hwinj
   banksim
   faithsim
   upload_to_gracedb

==========================================
Library Examples and Interactive Tutorials
==========================================

We have interactive tutorials and examples of using the pycbc library hosted
on azure jupyter notebooks. `Please give them a try. <https://notebooks.azure.com/nitz/libraries/pycbc>`_

In addition we have some examples below.

.. toctree::
   :maxdepth: 2

   catalog
   dataquality
   frame
   gw150914
   detector
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

Full Module Documentation is available at:

.. toctree::
   :maxdepth: 1

   modules

===================
Indexes and Tables
===================

* :ref:`modindex`
* :ref:`genindex`
