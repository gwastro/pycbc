.. PyCBC documentation master file, created by
   sphinx-quickstart on Tue Jun 11 17:02:52 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

===============
Getting Started
===============

PyCBC is a software package used to explore astrophysical sources of gravitational waves. It contains algorithms that can detect coalescing compact binaries and measure the astrophysical parameters of detected sources. PyCBC was used in the `first direct detection of gravitational waves (GW150914) by LIGO <https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.116.061102>`_ and is used in the ongoing analysis of LIGO and Virgo data.  If you have `Docker <https://www.docker.com/community-edition>`_ installed, you can get started using PyCBC with just two commands:

.. raw:: html

    <script type="text/javascript">
        document.addEventListener("DOMContentLoaded", function(){
            Typed.new(".element", {
                 strings: ["^500<strong>docker pull pycbc/pycbc-el7:latest</strong><br>$ ^500<strong>docker run -it pycbc/pycbc-el7:latest /bin/bash -l</strong><br>&#40;pycbc-software&#41;&#91;pycbc@37184573e664 &#126;&#93;$ ^500<strong>python</strong><br>Python 2.7.5 &#40;default, Nov  6 2016, 00:28:07&#41;<br>&#91;GCC 4.8.5 20150623 &#40;Red Hat 4.8.5-11&#41;&#93; on linux2<br>&gt;&gt;&gt; ^500<strong>import pycbc.version</strong><br>&gt;&gt;&gt; ^500<strong>print pycbc.version.git_tag</strong><br>v1.7.0<br>&gt;&gt;&gt; ^500<strong>import lal.git_version</strong><br>&gt;&gt;&gt; ^500<strong>print lal.git_version.id</strong><br>539c8700af92eb6dd00e0e91b9dbaf5bae51f004<br>&gt;&gt;&gt; "],
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

To use If you use PyCBC in your scientific publications or projects, we ask that you acknowlege our work by citing the papers described on the page:

.. toctree::
   :maxdepth: 1

   credit

===========
About PyCBC
===========

The goals of the PyCBC project are:

- Provide reliable and robust tools for building gravitational-wave search and parameter estimation workflows for CBCs.
- Create a flexible, extensible production code for CBC analysis that can be released for the public.
- Enable simple, easy and transparent access for various many-core architectures like GPUs.

==========================
Running PyCBC under Docker
==========================

The easiest way to start using PyCBC is to install one of our `Docker containers <https://hub.docker.com/u/pycbc/>`_. First, install the `Docker Community Edition <https://www.docker.com/community-edition>`_ for your `Mac <https://store.docker.com/editions/community/docker-ce-desktop-mac?tab=description>`_ or `Windows <https://store.docker.com/editions/community/docker-ce-desktop-windows?tab=description>`_ desktop. Docker CE installations for `Linux platforms <https://www.docker.com/community-edition#/download>`_ are also available.

After installing and starting Docker, type the commands::

    docker pull pycbc/pycbc-el7:v1.7.0
    docker run -it pycbc/pycbc-el7:v1.7.0 /bin/bash -l

 in the window above to download and run a PyCBC container. This example downloads the ``v1.7.0`` release in a container. To get the current version of the code from the `GitHub master branch <https://github.com/ligo-cbc/pycbc>`_ replace the string ``v1.7.0`` with ``latest``. The container includes all of the required software and dependencies to run PyCBC, including a compatible version of LALSuite.

If you want to run a docker container with graphics for making plots, you can
start a container that runs an SSH daemon and connect to that.  If do not
already have a personal ssh public/private key, first create one with the
command::

    ssh-keygen -t rsa
    cat ~/.ssh/id_rsa.pub >> ~/.ssh/authorized_keys
    chmod 600 ~/.ssh/authorized_keys

Now start the docker container with::

    docker run --name pycbc -d -P -v ${HOME}/.ssh:/home/pycbc/.ssh -t pycbc/pycbc-el7 tail -f /dev/null i   
    docker exec -u root -it pycbc /usr/bin/pycbc-sshd

And then you can connect to it with the command::

    ssh -Y pycbc@127.0.0.1 -p `docker port pycbc 22 | awk -F: '{print $NF}'`

Full installation instructions for users who want to install and develop PyCBC are available at:

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
