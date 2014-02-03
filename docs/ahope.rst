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

* Want to have the ability to run ligolw_segments_from_cats in separate-categories mode (and have this be the *only* mode that ahope will use). This requires some additions to ligolw_segments_compat so that vetoes called VETO_CAT_1, VETO_CAT_2 are merged into one list called CUMULATIVE_VETO_CAT_3. (Names are probably wrong, but you get the idea). - ALEX DOING?
* Adding ligolw_segments_compat breaks pegasus support and therefore the ability to generate vetoes within the workflow (instead of at runtime) has been temporarily disabled. This is because its input and output files are the *same name*, which pegasus cannot support. This needs fixing. This to-do item probably needs to done in parallel with the above item, and may need a new code, ligolw_make_cumulative_coincident_segments, or similar. - ALEX DOING?
* Dependent on above, once fixed we need to add some code for running ligolw_add and ligolw_segments_compat in the workflow if desired. - ALEX DOING?
* Enable the dax functionality (at the moment a dax is written, but we don't know what to do with it!) and ensure that codes are pegasus-compliant. (DUNCAN B + ALEX/IAN).
* The documentation, both in code and on this page, is woefully out of date. This needs fixing ... anyone/everyone should feel free to edit and improve the existing documentation!
* Merge functionality of ligolw_add into ligolw_sstinca to produce a single coincidence code that takes multiple inputs, ideally either from the command line or from a cache file.
* Ligolw_tisi takes repeated (and somewhat opaque) command line input e.g. ligolw_tisi -i H1:0,0,0 -i L1:0,100,5. This is not supported in pipeline.py. Edit ligolw_tisi so that it can be run within pipeline.py, or edit pipeline.py so that this is possible.
* Add options to CondorDAGJob so that condor commands can be provided as macro arguments (this would allow a number of features in the nodes).
* Add a function to CondorDAGNode to create_node from CondorDAGJob ... this would then allow us to use pipeline.py in the same way as ahope creates nodes from jobs. *BUT* does not break backwards compatibility. - IAN DOING
* Add in needed pipedown plots to ahope ... *There are no plans to tie plot_hipe into ahope*. This may require rewriting parts of a number of these codes, (especially to use the SQL where possible). The codes should also be moved away from pylal. The output format of the plots should be considered, do we still want to use the ihope plotting conventions or do we need something better.
* Edit ahope modules so that output directories are automatically created if they don't exist. - IAN DOING
* Add write_ihope_page to weekly ahope example. - IAN DOING
* Add ability to create a node as usual, but instead of adding it to the workflow have it run at ahope runtime. This should fail and refuse to run if it requires input files that are made within the workflow. This should probably be added at the pipeline.py level. - ALEX DOING?
* When above is completed, move segment calls in segment_utils to this and remove duplicate code paths. Also do same with ligolw_tisi call when tisi item above is also completed. - ALEX DOING?
* Begin moving pipedown codes into ahope.
* Add fixes to pipedown in cases where pipedown looks for information based on very specific command line calls made to earlier code, which are not made either in ahope or when using pycbc code in ihope.
* There is some not-https element in the pycbc docs homepage, which won't load by default on most browsers. Not sure what it is, but it should be fixed.

----------------------------
Longer term/lower priority
----------------------------

* If we want to move to having a variable number of segments for analysis (ie. if we only have 1000s use say 4 analysis segments, but if we have 5000s use 20), then we need to fix the issues of bias in the inverse PSD causing a bias (from the expected chi-squared distribution) in the output SNR time series, *which will vary based on the number of segments used to analyse the PSD*.

------------------
Questions concerning the coincidence stage ligolw_sstinca
------------------

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
daily_ahope
-------------------

This is an example of how to run an ahope workflow as far as it is coded. This
will:

* Get segments from the server, calculating CAT3 and above in the workflow for speed.
* Query the datafind server for frames.
* Create the template bank jobs needed to cover the times
* Split the template banks (into 2). This step could be easily removed, just delete this module in the python file and send the matched-filtering code the template bank input directly.
* Run the matched-filtering code on the split template banks
* Write a dax to file that can be submitted to run the workflow

This currently matches Chris' S6 test and works. Note that running this example (ie. time for the python code to do all of these steps and exit) takes approximately 1 minute to setup a workflow to analyse all of S6D. This is with the ping all frames step turned off, with this on we are seeing run times considerably longer.

-----------------------
er_daily_ahope
-----------------------

This is a more detailed example of daily_ahope. This will do the same things
as the daily_ahope example but then call into the remaining jobs, including
the daily_page dag generation, to set up a *complete* daily_ahope workflow.
This will also automatically generate the webpage at the end of the analysis.
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
Obtaining data
====================

Determine what data is available and what veto flags are active and use this to construct analysis and veto segment lists. Also check the data exists on the machine that will be performing the analysis.

.. toctree::
   :maxdepth: 1

   ahope/segments

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
