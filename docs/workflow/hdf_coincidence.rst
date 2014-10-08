.. _workflowintervalcoincmod:

##################################

##################################

=============
Introduction
=============

On this page we go over how to use the interval slides based post processing
functionality in PyCBC.

.. note:: This functionality is experimental and under active development. 
    No gaurantee of stability is given at this time. Contact Alex Nitz (alex.nitz@ligo.org) 
    if you have questions, issues, comments, etc. 
    
The post processing discussed here implements coincidence, background estimation,
injection finding, and some basic plotting infrastructure based on hdf files. 

This document assumes that one is familiar with the standard pycbc workflow modules.
Only differences from the standard configuration are discussed here.

Familiarize yourself with the standard modules.

.. toctree::
   :maxdepth: 1

   ../weekly_ahope

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

----------------------
Workflow Configuration
----------------------

The following are the sections used to configure the workflow.

* [workflow-coincidence]

There are two options, both of which are permorfance options.

* --number-of-groups

This is integer that sets how many groups of templates to analyze at once
for coincidence. A reasonable choice seems to be 1 per 100 templates in your bank.
This will directly correlate the to memory usage of each coincidence job.

* --groups-per-coinc

This sets how many template groups per coincidence job. This knob helps control
how long each coincidence job should take. A number of ~ 50 seems to give
fairly short running jobs. 

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
    decimation-keep = 100
    ; For now, don't touch these options, but they control how to decimate low
    ; significane triggers.
    decimation-factor = 1000
    decimation-bins = 1

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
    cluster-window = 0.2

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Options to the injection finding executable
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
::

    [hdfinjfind]
    ; Injections are associated to triggers by time alone. The one option
    ; determines the size in seconds of the window. Make sure this is smaller
    ; then the coincident clustering window, otherwise note that injections
    ; found in association with multiple triggers are simply marked as missed.
    injection-window = .05

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

=====================================
Workflow Generation an Planning
=====================================

Configure a two-detector search using a single fixed template bank for all ifos
and the entire analysis. Follow the standard module configurations up through
the inspiral stage of the workflow.

For module configuation documentation.
.. toctree::
   :maxdepth: 1

   ../weekly_ahope

To generate the workflow::

    GPS_START_TIME=966384015
    GPS_END_TIME=967384015
    export LOGPATH=/usr1/${USER}/log
    mkdir -p $LOGPATH

    pycbc_make_nitz_coinc_workflow \
    --output-dir gwanalysis \
    --local-config-files ncoinc.ini \
    --config-overrides \
    workflow:start-time:${GPS_START_TIME} \
    workflow:end-time:${GPS_END_TIME} 

To plan the workflow::

    cd gwanalysis
    pycbc_basic_pegasus_plan weekly_ahope.dax $LOGPATH

==============================================================
Reusing data from workflow that use some other post-processing
==============================================================

Assuming the list of limitations is satisfied by the previous run, then one
can simply select the 'weekly_ahope.map' file in the '--cache' option to the
pegasus planner.::

    cd gwanalysis
    pycbc_basic_pegasus_plan weekly_ahope.dax $LOGPATH --cache /path/to/prior/worklfow/weekly_ahope.map

If the prior workflow did not use the post-processing described on this page, then there is no need
to edit the map file and it can be used as is.

If you are rerunning a workflow using this post-processing, then select from your weekly_ahope.map file only
the files you want to reuse, and then point the '--cache' option to it instead.

