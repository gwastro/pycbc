###################################################
Ahope: the inspiral analysis workflow generator (``pycbc.ahope``)
###################################################

============
Introduction
============

Ahope is a tool used to create the workflows needed to perform coincident analyses of gravitational-wave data searching for compact-binary-coalescences using matched-filtering.

The ahope workflow runs through a number of stages, which are designed to be as independent of each other as possible, while also allowing an integrated end-to-end workflow to be constructed.

Documentation of the ahope executable and how to run it can be found here.

Each of the sections is described in detail below. Also refer to the `Link page here <https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/InspiralPipelineDevelopment/130627093845InspiralPipelineDocumentationAhope_development_plan>`_ for details. Also see the `Link page here <https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/InspiralPipelineDevelopment/130614025403InspiralPipelineDevelopmentThe%20evolution%20of%20ihope%20and%20inspiral>`_ for descriptions of what ahope should be and motivation for making ahope in the first place. **NOTE: REMOVE THESE LINKS ONCE AHOPE DOCUMENTATION IS SUFFICIENT THAT THIS IS NOT LONGER NEEDED**

The sections are 

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

.. autofunction:: pycbc.ahope
