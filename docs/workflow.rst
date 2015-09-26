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

Each of the sections is described in detail below.

=======================
Workflow to do list
=======================

This list is now maintained `as a wiki <https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/AhopeToDoList>`_; please see that wiki to view items on the list or update their status. **PLEASE** take on an item from this list. If you do, please update the wiki and ensure that any changes are reflected in both the in-line documentation and the html web pages.

====================
Examples
====================

Some examples of how to run the workflow module
in different configurations can be found
in the source tree under::

    pycbc/examples/workflow

These examples are described in each section below

-----------------------------
``pycbc_make_coinc_workflow``
-----------------------------

This is an example of how to run a coincidence workflow, mimicking standard ihope coincidence analysis as closely as possible. It calls into pipedown to do the post-processing and using `lalapps_write_ihope_page` to make an output html page. In total this will:

* Get science and data-quality segments from the server.
* Query the datafind server for frames.
* Create the template bank jobs needed to cover the times
* Split the template banks (into whatever number is specified in the ini file). This step could be easily removed, just delete this module in the python file and send the matched-filtering code the template bank input directly.
* Run the matched-filtering code on the split template banks
* Use ligolw_add and ligolw_sstinca to generate coincidence triggers
* Do some compatibility conversions and then call pipedown to create a dag to do the post-processing
* Add a node to run write_ihope_page at the end of analysis
* Native post-processing is also present in a testing mode, this will not be seen in the default output page, but does include a full set of plots.
* Write a dax to file that can be submitted to run the workflow

This will therefore set up an *almost complete* mimic of a weekly ihope analysis and automatically generate the output webpage at the end of the analysis.

More details of how to run this is found in the following link:

.. toctree::
   :maxdepth: 1

   workflow/pycbc_make_coinc_workflow

----------------------------
``pycbc_make_sngl_workflow``
----------------------------

This is an example of how to run a single-detector workflow done in the past. This will:

* Get science and data-quality segments from the server.
* Query the datafind server for frames.
* Setup the template bank jobs needed to cover the times.
* Split the template banks. This step could be easily removed, just delete this module in the python file and send the matched-filtering code the template bank input directly.
* Setup the matched-filtering code on the split template banks.
* Setup the plotting code.
* Write a dax to file that can be submitted to run the workflow

This will therefore set up a *complete* single-detector analysis workflow and automatically generate the webpage at the end of the analysis.

More details of how to run ``pycbc_make_sngl_workflow`` is found in the following link:

.. toctree::
   :maxdepth: 1

   workflow/pycbc_make_sngl_workflow

-----------------------------
``pycbc_make_sngl_workflow``
-----------------------------

Another single-detector workflow generator that mimics the daily ihope analysis is ``pycbc_make_sngl_workflow``. More details of how to run ``pycbc_make_sngl_workflow`` is found in the following link:

.. toctree::
   :maxdepth: 1

   workflow/pycbc_make_sngl_workflow

---------------------
`data_checker`
---------------------

This is an example of using only the segment query and datafind query modules
of the workflow module, with all checking options enabled. This will do the following

* Use the supplied information in the veto-definer file query for science times
* Using the veto-definer file query for CAT1,2,3,4,5 times
* Connect to the datafind server and query for frames at the times given in science - CAT1 (frame types provided in .ini file)
* Check that frames have been returned at all times that are given as science (-CAT1) by the segment server. Fail if not (this can be changed to update the science times, but here we want to fail if data is not present)
* Run os.path.is_file() on every frame file that was returned by datafind. This can be slow if run on long periods of time and if the fileservers are slow. Not recommended for a full coincidence workflow at present, but useful for daily [a,i]hope runs.

This example can be edited for any time you are interested in and can be used to identify missing data before submitting a workflow.

---------------------
``pygrb``
---------------------

This is an example of how to generate a coherent network segment for analysis in a targeted, coherent workflow in response to an external trigger, such as a gamma-ray burst. If submitted, single detector coincidence analysis jobs will be run, however coherent matched filtering is in development. In summary, this will currently:

* Get science and data-quality segments from the server.
* Query the datafind server for frames.
* Create single detector template bank jobs needed to cover the search parameter space, or use a pregenerated bank (recommended).
* Split the template banks according to options given in the ini file.
* Run single detector coincident matched-filtering code on the split template banks.
* Write a dax to file that can be submitted to run the workflow.

More details of how to run this is found in the following link:

.. toctree::
   :maxdepth: 1

   workflow/pygrb

==============================
Workflow module documentation
==============================

The following contains a list of the sub-modules in pycbc's workflow module, a bried description of what each does, and a link to documentation on how to use each module. Collectively these pages should provide complete details on how to set up a workflow from scratch.

--------------------
Initialization
--------------------

Take in command line input and the configuration file.

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
Time-slide generation
----------------------------------

Generate the time-slide files used for describing time slides to be performed later in the analysis

.. toctree::
   :maxdepth: 1

   workflow/time_slides

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

   workflow/coincidence

----------------------------------
Post processing preparation
----------------------------------

Prepare output files for the post-processing stage

.. toctree:: :maxdepth: 1

   workflow/postprocprep

----------------------------------
Post-processing
----------------------------------

Assess the significance of all triggers and calculate any rate statements

.. toctree:: :maxdepth: 1

   workflow/postproc

====================
Method documentation
====================

The documentation for all functions/modules within pycbc.workflow follows, unless you are looking for a specific function, it might be easier to navigate through the section links above.

.. toctree::
    pycbc.workflow
