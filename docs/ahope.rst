##################################################################
Ahope: the inspiral analysis workflow generator (``pycbc.ahope``)
##################################################################

============
Introduction
============

Ahope is a tool used to create the workflows needed to perform coincident analyses of gravitational-wave data searching for compact-binary-coalescences using matched-filtering.

The ahope workflow runs through a number of stages, which are designed to be as independent of each other as possible, while also allowing an integrated end-to-end workflow to be constructed.

Documentation of the ahope executable and how to run it can be found here.

Each of the sections is described in detail below. Also refer to the `page here <https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/InspiralPipelineDevelopment/130627093845InspiralPipelineDocumentationAhope_development_plan>`_ for details. Also see the `link here <https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/InspiralPipelineDevelopment/130614025403InspiralPipelineDevelopmentThe%20evolution%20of%20ihope%20and%20inspiral>`_ for descriptions of what ahope should be and motivation for making ahope in the first place. **NOTE: REMOVE THESE LINKS ONCE AHOPE DOCUMENTATION IS SUFFICIENT THAT THIS IS NOT LONGER NEEDED**

=======================
Ahope to do list
=======================

This list is now maintained `as a wiki <https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/AhopeToDoList>`_; please see that wiki to view items on the list or update their status. **PLEASE** take on an item from this list. If you do, please update the wiki and ensure that any changes are reflected in both the in-line documentation and the html web pages. Not taking on any items will result in owing Ian beers at the Nice meeting.

-----------------------------------------------------------
Questions concerning the coincidence stage ligolw_sstinca
-----------------------------------------------------------

* IHOPE: Using .from_url() on the ihope filenames and then using the resulting segment as the time that file is valid from is incorrect as ihope does not analyse the first/last 64s of each file (sometimes even less than that if trig-start-time or trig-end-time are specified). Is this a problem?
* IHOPE: When setting start/end for thinca it basically sets them to cover *all* triggers in the input files. Sometimes, when segments are > 10000s this should be split over multiple jobs, so you might only want to read in *some* triggers from files at the end of the job. Would this cause duplicate triggers? Is this a problem at all?
* AHOPE: Sstinca doesn't work unless it is given a veto file. I think this is because it is using it to get the ifos. This is potentially dangerous with ahope, which is *not* guaranteed to have entries for all ifos!
* BOTH: Handling of veto information in sstinca : the job takes the veto segments as input, then removes them from the output xml.  This is in order to avoid massive duplication of the veto definer and segment tables in thinca output files, which causes pipedown to choke when creating and simplifying databases.  Pipedown puts the veto/segment tables back in again in a single step. 


====================
Examples
====================

Some examples of how to run ahope in different configurations can be found
in the source tree under::

    pycbc/examples/ahope

These examples are described in each section below

---------------------
data_checker
---------------------

This is an example of using only the segment query and datafind query modules
of ahope, with all checking options enabled. This will do the following

* Use the supplied information in the veto-definer file query for science times
* Using the veto-definer file query for CAT1,2,3,4,5 times
* Connect to the datafind server and query for frames at the times given in science - CAT1 (frame types provided in .ini file)
* Check that frames have been returned at all times that are given as science (-CAT1) by the segment server. Fail if not (this can be changed to update the science times, but here we want to fail if data is not present)
* Run os.path.is_file() on every frame file that was returned by datafind. This can be slow if run on long periods of time and if the fileservers are slow. Not recommended for a full ahope run at present, but useful for daily [a,i]hope runs.

This example can be edited for any time you are interested in and can be used to identify missing data before submitting a workflow. 

-------------------
er_daily_ahope
-------------------

This is an example of how to run an ahope workflow mimicing the daily_ihope analysis done in the past. This will:

* Get science and data-quality segments from the server.
* Query the datafind server for frames.
* Create the template bank jobs needed to cover the times
* Split the template banks (into 2). This step could be easily removed, just delete this module in the python file and send the matched-filtering code the template bank input directly.
* Run the matched-filtering code on the split template banks
* Call into the code that generates the dag for running daily_ihope pages and plots and adds this dag to the workflow.
* Write a dag/dax to file that can be submitted to run the workflow

This will therefore set up a *complete* daily_ahope workflow and automatically generate the webpage at the end of the analysis.

This is currently being used in ER5 as a supplement and backup to the daily_ihope runs. 

More details of how to run this is found in the following link:

.. toctree::
   :maxdepth: 1

   ahope/daily_ahope

====================
Initialization
====================

Take in command line input and the configuration file.

.. toctree::
   :maxdepth: 1

   ahope/initialization

====================
Generating segments
====================

Obtain the science segments and data-quality segments from making queries to a segment database.

.. toctree::
   :maxdepth: 1

   ahope/segments

=====================
Obtaining data
=====================

Run queries to the datafind server to find the needed frames and test these for
consistency if desired

.. toctree::
   :maxdepth: 1

   ahope/datafind

======================
Injection generation
======================

Generate injection files for use later in the analysis

.. toctree::
   :maxdepth: 1

   ahope/injections

========================
Time-slide generation
========================

Generate the time-slide files used for describing time slides to be performed later in the analysis

.. toctree::
   :maxdepth: 1

   ahope/time_slides

====================
Template bank
====================

Construct a template bank, or banks, of CBC waveforms that will be used to matched-filter the data with.

.. toctree::
   :maxdepth: 1 

   ahope/template_bank

====================
Matched-filtering
====================

Perform the matched-filters and calculate any signal-based consistency tests that should be calculated.

.. toctree::   :maxdepth: 1 

   ahope/matched_filter

====================
Coincidence
====================

Determine if "triggers" seen in one detector are also seen in other detectors. Also check for coincidence between time-slid triggers for background evaluation

.. toctree::   :maxdepth: 1 

   ahope/coincidence

====================
Post-processing
====================

Interpret the results, calculating false-alarm probabilities and rate limits.

.. toctree::   :maxdepth: 1 

   ahope/post_processing

====================
Plotting
====================

Create plots of the various outputs through the pipeline to easily enable the user to view results

.. toctree::   :maxdepth: 1

    ahope/plotting


====================
Generate html pages
====================


Tie the plots together in one html page, where *all* desired information should be available

.. toctree::   :maxdepth: 1

    ahope/web_pages


====================
Method documentation
====================

The documentation for all functions/modules within pycbc.ahope follows, unless you are looking for a specific function, it might be easier to navigate through the section links above.

.. toctree::
    pycbc.ahope
