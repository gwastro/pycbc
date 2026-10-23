.. _workflowhomepage:

########################################################################
Workflow: the inspiral analysis workflow generator (``pycbc.workflow``)
########################################################################

============
Introduction
============

Pycbc's workflow module is a tool used to create the workflows needed to perform coincident analyses of gravitational-wave data searching for compact-binary-coalescences using matched-filtering.

The pycbc workflow runs through a number of stages, which are designed to be as independent of each other as possible, while also allowing an integrated end-to-end workflow to be constructed.

Documentation of the workflow module and how to run it can be found below.

Please see the `following poster <https://dcc.ligo.org/LIGO-G1400223>`_, presentied at the March LVC meeting, 2014, for an introduction to the workflow module (referred to as "ahope"). Especially see the following top-level workflow generation model.

.. image:: images/workflow_planning.png
   :width: 100 %

==============================
Workflow module documentation
==============================

The following contains a list of the sub-modules in pycbc's workflow module, a bried description of what each does, and a link to documentation on how to use each module. Collectively these pages should provide complete details on how to set up a workflow from scratch.

---------------------
Basics and overview
---------------------

The following page gives a description of how the workflow module expects
configuration files to be layed out and how the basic command-line interface
functions. If writing a workflow for the first time we recommend you read
through this page!

.. toctree::
   :maxdepth: 1

   workflow/initialization

--------------------
Generating segments
--------------------

Obtain the science segments and data-quality segments from making queries to a segment database.

.. toctree::
   :maxdepth: 1

   workflow/segments

--------------------
Obtaining data
--------------------

Run queries to the datafind server to find the needed frames and test these for
consistency if desired

.. toctree::
   :maxdepth: 1

   workflow/datafind

----------------------------
Injection generation
----------------------------

Generate injection files for use later in the analysis

.. toctree::
   :maxdepth: 1

   workflow/injections

----------------------------------
Template bank
----------------------------------

Construct a template bank, or banks, of CBC waveforms that will be used to matched-filter the data with.

.. toctree::
   :maxdepth: 1 

   workflow/template_bank

----------------------------------
Split table
----------------------------------

Split an output file into numerous parts to allow parallel analysis. Normally used to split the template bank up to allow matched-filtering in parallel

.. toctree::
   :maxdepth: 1 

   workflow/splittable

----------------------------------
Matched-filtering
----------------------------------

Perform the matched-filters and calculate any signal-based consistency tests that should be calculated.

.. toctree::   :maxdepth: 1 

   workflow/matched_filter

----------------------------------
Coincidence
----------------------------------

Determine if "triggers" seen in one detector are also seen in other detectors. Also check for coincidence between time-slid triggers for background evaluation

.. toctree::   :maxdepth: 1 

   workflow/hdf_coincidence

====================
Method documentation
====================

The documentation for all functions/modules within pycbc.workflow follows, unless you are looking for a specific function, it might be easier to navigate through the section links above.

.. toctree::
    pycbc.workflow
