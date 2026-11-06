.. _workflowintervalcoincmod:

##################################
HDF5 Based Coincidence Code
##################################

=============
Introduction
=============

On this page we go over how to use the interval slides based post processing
functionality in PyCBC.

.. note:: 
    This functionality is experimental and under active development. 
    No gaurantee of stability is given at this time. Contact Alex Nitz (alex.nitz@ligo.org) 
    if you have questions, issues, comments, etc. 
    
The post processing discussed here implements coincidence, background estimation,
injection finding, and some basic plotting infrastructure based on hdf files. 

This document assumes that one is familiar with the standard pycbc workflow modules.
Only differences from the standard configuration are discussed here.


==============
Limitations
==============
There are several limitations that may not be familiar to people used to existing
post-processing codes. 

The following restrictions are not slated to be changed in the near future.

* Coincidence is exact match only

The following restrictions are slated to be be removed when time permits

* Only two detector searches are implemented
* The significance estimation does not include any parameter space binning
* Only a single fixed template bank for the entire anlysis is possible

===================
Configuration File
===================

An example configuration file is included as example_hdf_post.ini. This acts as a supplement to the standard ini files. It does not stand alone.

----------------------
Workflow Configuration
----------------------

The following are the sections used to configure the workflow.

* [workflow-coincidence]

There is one option related to the performance of the background estimation and coincidence.

* --parallelization-factor

This is an integer that determines the total number of jobs to calculate the foreground and background
coincidences for the entire run. The number is inversely proportional to the memory usage of each job, and the expected run time.
A factor of 10 is sufficient for a one month run with an SNR threshold of 5.0, 8000 templates, and a timeslide every 1.1 seconds. 
To ensure short running jobs, the factor should be increased in proportion to the number of background triggers. 
The memory usage of each jobs is only porportional to the fraction of single detector triggers that are read in. 

------------------------
Executable Configuration
------------------------

There following executables are used in the coincidence, injection finding, and 
plotting. They following needs to be in your [executables] section.::

    bank2hdf = ${which:pycbc_coinc_bank2hdf}
    trig2hdf = ${which:pycbc_coinc_trig2hdf}
    hdfinjfind = ${which:pycbc_coinc_hdfinjfind}
    coinc = ${which:pycbc_coinc_findtrigs}
    statmap = ${which:pycbc_coinc_statmap}
    plot_sensitivity = ${which:pycbc_page_sensitivity}
    plot_foundmissed = ${which:pycbc_page_foundmissed}
    plot_snrifar = ${which:pycbc_page_snrifar}
    plot_snrchi = ${which:pycbc_page_snrchi}
    page_foreground = ${which:pycbc_page_foreground}
    hdf_trigger_merge = ${which:pycbc_coinc_mergetrigs}

Executables that do not require any additional configuration.

* bank2hdf : This converts the xml based template bank file to an hdf format.
* trig2hdf : This converts the xml based inspiral trigger file to an hdf file. If you have used a split template bank they are recombined into a single hdf file.
* plot_foundmissed : Make an interactive found missed plot. 
* plot_snrifar : Makes ifar vs SNR cumulative plots

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Options to the coincidence/background executable 
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
::

    [coinc]
    ; Coincidence is based on the time-of-flight. This options sets how many
    ; seconds to add to this fixed window.
    coinc-threshold = .005
    ; Background coincident triggers are decimated. This option controls how
    ;many of the loudest triggers to keep from each group of templates.
    loudest-keep = 200
    ; For now, don't touch these options, but they control how to decimate low
    ; significane triggers.
    decimation-factor = 1000

    [coinc-inj]
    [coinc-full]
    ; The time interval in seconds between each timeslide. 
    timeslide-interval=1.9

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Options to the coincidence clustering + statistic asignment executable
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
::

    [statmap]
    ; The only option is to set the time window in seconds to cluster 
    ; coincident triggers over the entire template bank
    cluster-window = 5.0

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Options to the injection finding executable
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
::

    [hdfinjfind]
    ; Injections are associated to triggers by time alone. The one option
    ; determines the size in seconds of the window. Make sure this is smaller
    ; then the coincident clustering window, otherwise note that injections
    ; found in association with multiple triggers are simply marked as missed.
    injection-window = 1.0

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Options for sensitivity plotting
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
::

    [plot_sensitivity]
    ; The type of bins to use. Options are mchirp, total_mass, or spin.
    bin-type =  mchirp 
    ; The bin boundaries
    bins = 0.89 1.31 1.74 2.17 2.60 
    min-dist = 40 
    max-dist = 120 
    ; The number of distance bins. If there is no injection in a bin, it fails, as it should!
    ; Remove this option to use binless sensitivity estimation. The distance bins options
    ; only exists for comparison to other pipelines where the codes do distance bins.
    dist-bins = 15 

Additional plots can be made by adding a tag. This works similary to the system for injection files. For 
example you can add  [plot_sensitivity-mchirp], [plot_sensitivity-mtotal], and [plot_sensitivity-spin]
sections to make three versions of the plot.

