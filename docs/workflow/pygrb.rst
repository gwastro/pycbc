###############################################################################
PyGRB: A GRB triggered CBC analysis workflow generator
###############################################################################

===============
Introduction
===============

PyGRB is a tool used to generate a data analysis workflow for a targeted,
coherent gravitational wave search triggered by short duration gamma-ray
bursts.

When submitted, this workflow will run a pipeline to analyse data from multiple
gravitational wave detectors coherently. It will then perform various
signal-based veto cuts and data quality cuts to determine whether or not a
compact binary coalescence signal is present in the given data coming from the
same point in the sky and at the same time as an observed short duration
gamma-ray burst.

The output will be a webpage containing plots and other data files that can be
used to understand the results of the analysis.

.. _howtorunpygrb:

=======================
How to run
=======================

Here we document the stages needed to run the triggered coherent GRB search.

Once PyCBC is installed, you should be able to run the following help command for the workflow generation
script::

    pycbc_make_offline_grb_workflow --help

This should produce a help message like the following

.. command-output:: pycbc_make_offline_grb_workflow --help

This outlines the command line arguments that may be passed to the executable.
The majority of options passed to the workflow will come from configuration
files, and these are known to the executable via the option
--config-files.

----------------------
Set up a run directory
----------------------

Navigate to the directory you wish to run in::

    RUN_DIR=/path/to/run/directory
    mkdir -p $RUN_DIR
    cd $RUN_DIR

Next gather together configuration files for your run.

-----------------------------------------------------------------------
Configuration files - Do you already have configuration (.ini) file(s)?
-----------------------------------------------------------------------
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Yes, I already have configuration files
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

Copy over configuration files into your run directory::

    cp /path/to/config_file1.ini /path/to/config_file2.ini $RUN_DIR

Now add these configuration files to your environment. If you have more than
one configuration file they must be space separated::

    LOCAL_CONFIG_FILES="config_file1.ini config_file2.ini"

Now go down to :ref:`pygrbgenerate`.

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
No, I need to edit a configuration file
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

The default configuration files are found in::

    /src/pycbc/examples/workflow/pygrb/

namely::

    main.ini injections.ini postprocessing.ini

These files contain almost all options needed to construct an example PyGRB
workflow

.. note::

    If you are unfamiliar with pycbc workflows, look through these files.
    
* main.ini contains options that are used when setting up the data retreival,
  template bank generation and main analysis matched filtering parts of the
  workflow

* injections.ini contains options for generating injection sets

* postprocessing.ini contains options that are used to run the post processing
  stages of the workflow codes

Currently, the example is set up to run on S6 data, analysing only H1 and L1.

If you want to run in this default configuration please jump down to
:ref:`pygrbgenerate`.

If you want to run on non-S6 data, analyze a different set of ifos, or change
any data-type or segment length options, you will have to edit some additional
options.

========
main.ini
========

Firstly, the sections specifying some workflow-wide options include::

    [workflow]
    h1-channel-name = H1:LDAS-STRAIN
    l1-channel-name = L1:LDAS-STRAIN

    [workflow-ifos]
    ; This is the list of ifos to analyse
    h1 =
    l1 =


    [workflow-datafind]
    datafind-h1-frame-type = H1_LDAS_C02_L2
    datafind-l1-frame-type = L1_LDAS_C02_L2

    [workflow-segments]
    segments-h1-science-name = H1:DMT-SCIENCE:4
    segments-l1-science-name = L1:DMT-SCIENCE:4
    segments-database-url = https://segdb.ligo.caltech.edu
    segments-veto-definer-url = https://www.lsc-group.phys.uwm.edu/ligovirgo/cbc/public/segments/S6/H1L1V1-S6_CBC_LOWMASS_D_OFFLINE-961545543-0.xml
    segments-veto-categories = 3
    segments-minimum-segment-length = 256

    [workflow-exttrig_segments]
    ; options for the coherent search (development)
    on-before = 5
    on-after = 1
    min-before = 60
    min-after = 60
    min-duration = 256
    max-duration = 5264
    quanta = 128
    num-buffer-before = 8
    num-buffer-after = 8

To run through this

* The [workflow-ifos] section supplies which ifos will be analysed if data is
  found and available
* The [workflow-exttrig_segments] section supplies the GRB search-specific
  options for the data segment to be analyzed
* The xx-channel-name options are the h(t) channel name in the frames
* The datafind-xx-frame-type is the type of the frames for use when calling
  gw_data_find
* The segments-xx-science-name is the flag used to store science times in the
  segment database
* segments-database-url points to the segment database
* segments-veto-definer-url points to the url where the veto-definer file can
  be found.

We also set the executables to be used for these parts of the analysis::
    
    [executables]
    ; setup of condor universe and location of executables
    tmpltbank               = ${which:lalapps_tmpltbank_ahope}
    inspiral                = ${which:lalapps_coh_PTF_inspiral}
    splitbank               = ${which:pycbc_splitbank}
    segment_query           = ${which:ligolw_segment_query}
    segments_from_cats      = ${which:ligolw_segments_from_cats}
    llwadd                  = ${which:ligolw_add}
    ligolw_combine_segments = ${which:ligolw_combine_segments}

The options to be given to every job run by an executable are then given
within a secion with the relevant name, for example our inspiral jobs (in this
case, lalapps_coh_PTF_inspiral) use the options in the following section::

    [inspiral]
    ; coh_PTF_inspiral analysis parameters -- added to all inspiral jobs
    ; Note that some values are dynamically recalculated during workflow generation
    .
    .
    .

These should not be edited unless you know what you are doing. To find out
more details about the possible options for any stage of the workflow, follow
the links at :ref:`workflowhomepage`.

==============
injections.ini
==============

Multiple configuration files may be used, and in fact the same sections may be
populated from within multiple files. As an example, we might wish to have a
seperate file for the injection options. This file may contain the following::

    [executables]
    injections       = ${which:lalapps_inspinj}
    jitter_skyloc    = ${which:ligolw_cbc_jitter_skyloc}
    align_total_spin = ${which:ligolw_cbc_align_total_spin}
    split_inspinj    = ${which:pycbc_split_inspinj}

This will add to the values given in the [executables] section of the other 
file. Options for the trig_combiner code may then be given in a section 
[trig_combiner], and so on. Options in sections such as [injections-<tag>]
will be passed to the injection executable and create jobs tagged with <tag>.
This can be used, as in the example, to generate multiple injection sets.
Options in sections such as [workflow-injections-<tag>] will apply to the same
tagged injection set, but are not passed to the executable. They can instead be
used to control the behaviour of the workflow generation as applied to that
specific tagged injection set. In our example these options are the number of
injections to be included in each set.

==================
postprocessing.ini
==================

As before, the options for each of the post processing codes are given in this
configuration file. One notable section is::

    [pegasus_profile-trig_combiner]
    condor|request_memory=2000M

This is how to supply condor options that only apply to the trig_combiner jobs.
This can be generalised to any executable or tagged jobs.

===============
Copy over files
===============

Now you have a configuration file (or files) and 
can follow the same instructions as above. That is: 

Copy the configuration file into your run directory::

    cp /path/to/<file(s)>.ini .

and set the name of the configuration file in your path. If you have more than
one configuration file they must be space separated::

    LOCAL_CONFIG_FILES="main.ini injections.ini postprocessing.ini"

.. _pygrbgenerate:

-----------------------
Generate the workflow
-----------------------

When you are ready, you can generate the workflow. This may be done by setting
a number of variables in your environment before launching the generation
script.

First we need to choose a trigger time, ie. the GPS Earth-crossing time
of the GRB signal. You should also set the GRB name. For example::

    GRB_TIME=969675608
    GRB_NAME=100928A

We should next set the sky coordinates of the GRB in RA and Dec, in this
example::

    RA=223.0
    DEC=-28.5
    SKY_ERROR=0

If you are using a pregenerated template bank and do not have a path to the
bank set in your config file, set it here::

    BANK_FILE=path/to/templatebank

You also need to specify the git directory of your lalsuite install::

    export LAL_SRC=/path/to/folder/containing/lalsuite.git

If you want the results page to be moved to a location outside of your run,
provide this too::

    export HTML_DIR=/path/to/html/folder

If you are using locally editted or custom configuration files then you can
create the workflow from within the run directory using::

    pycbc_make_offline_grb_workflow \
             --config-files ${LOCAL_CONFIG_FILES} \
             --config-overrides workflow:ra:${RA} \
                                workflow:dec:${DEC} \
                                workflow:sky-error:${SKY_ERROR} \
                                workflow:trigger-name:${GRB_NAME} \
                                workflow:trigger-time:${GRB_TIME} \
                                workflow:start-time:$(( GRB_TIME - 4096 )) \
                                workflow:end-time:$(( GRB_TIME + 4096 )) \
                                workflow:html-dir:${HTML_DIR} \
                                workflow-tmpltbank:tmpltbank-pregenerated-bank:${BANK_FILE}

This may all be conveniently placed within a shell script, an example of which is given in::

    /src/pycbc/examples/workflow/pygrb/run_pygrb.sh

.. _pygrbplan:

-----------------------------------------
Planning and Submitting the Workflow
-----------------------------------------
CD into the directory where the dax was generated::

    cd GRB${GRB_NAME}

From the directory where the dax was created, run the submission script::

    pycbc_submit_dax --dax pygrb.dax --accounting-group <your.accounting.group.tag>

.. note::

    If running on the ARCCA cluster, please provide a suitable directory via
    the option --local-dir, ie. /var/tmp

-------------------------------------------------------------------------------------------------------------------------------------------
Monitor and Debug the Workflow (`Detailed Pegasus Documentation <https://pegasus.isi.edu/wms/docs/latest/tutorial.php#idm78622034400>`_)
-------------------------------------------------------------------------------------------------------------------------------------------

To monitor the above workflow, one would run::

    pegasus-status /path/to/analysis/run
    
To get debugging information in the case of failures.::

    pegasus-analyzer /path/to/analysis/run

=============================
Workflow visualization
=============================

-----------------------------
Pegasus Dashboard
-----------------------------

The `pegeasus dashboard <http://pegasus.isi.edu/wms/docs/latest/ch02s11.php>`_
is a visual and interactive way to get information about the progress, status,
etc of your workflows.

The software can be obtained from a seprate pegasus package here
<https://github.com/pegasus-isi/pegasus-service>.

-----------------------------
Pegasus Plots
-----------------------------


Pegasus has a tool called pegasus-plan to visualize workflows. To generate
these charts and create an summary html page with this information, one would
run::

    export PPLOTSDIR=${HTMLDIR}/pegasus_plots
    pegasus-plots --plotting-level all --output ${PPLOTSDIR} /path/to/analysis/run

The Invocation Breakdown Chart section gives a snapshot of the workflow. You
can click on the slices of the pie chart and it will report the number of
failures, average runtime, and max/min runtime for that type of jobs in the
workflow. The radio button labeled runtime will organize the pie chart by total
runtime rather than the total number of jobs for each job type.

The Workflow Execution Gantt Chart section breaks down the workflow how long it
took to run each job. You can click on a job in the gantt chart and it will
report the job name and runtime.

The Host Over Time Chart section displays a gantt chart where you can see what
jobs in the workflow ran on a given machine.

.. _pygrbreuse:

================================
Reuse of workflow file products
================================

One of the features of  Pegasus is to reuse the data products of prior runs.
This can be used to expand an analysis or recover a run with mistaken settings
without duplicating work.

-----------------------------------------
Generate the full workflow you want to do
-----------------------------------------

First generate the full workflow for the run you would like to do as normal,
following the instructions of this page from :ref:`howtorunpygrb`, but stop
before planning the workflow in :ref:`pygrbplan`.

-----------------------------------------------------
Select the files you want to reuse from the prior run
-----------------------------------------------------

Locate the directory of the run that you would like to reuse. There is a file
called GRB${GRB_NAME}/output.map, that contains a listing of all of the data
products of the prior workflow.

Select the entries for files that you would like to skip generating again and
place that into a new file. The example below selects all the inspiral and 
tmpltbank jobs and places their entries into a new listing called
prior_data.map.::

    # Lets get the tmpltbank entries
    cat /path/to/old/run/GRB${GRB_NAME}/output.map | grep 'TMPLTBANK' > prior_data.map
    
    # Add in the inspiral  files
    cat /path/to/old/run/GRB${GRB_NAME}/output.map | grep 'INSPIRAL' >> prior_data.map

.. note::

    You can include files in the prior data listing that wouldn't be generated
    anyway by your new run. These are simply ignored.

---------------------------
Plan the workflow
---------------------------

From the directory where the dax was created, run the planning script::

    pycbc_submit_dax --dax pygrb.dax --accounting-group <your.accounting.group.tag> --cache-file /path/to/prior_data.map

Follow the remaining :ref:`pygrbplan` instructions to submit your reduced
workflow.

