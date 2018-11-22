####################################################################################
``pycbc_make_offline_grb_workflow``: A GRB triggered CBC analysis workflow generator
####################################################################################

============
Introduction
============

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

==================
Configuration File
==================

The workflow is controlled by a configuration file, which is comprised of three
types of section: workflow, the pegasus profile and executable options. The
workflow sections control the general form of the workflow and how it is
generated. The pegasus profile sections are equivalent to lines you would have
in a condor_submit file (e.g. requirements, storage size, etc.). The executable
option sections contain those options to be fed directly to the executables
that will be used for the data analysis.

The following examples would all feature in an offline search on combined
Hanford-Livingston data from the 8th Advanced LIGO engineering run (ER8). To
find out more details about the possible options for any stage of the workflow,
follow the links at :ref:`workflowhomepage`.

-----------------
Workflow Sections
-----------------

The main section usually contains overall properties for your workflow, many of
which will be used elsewhere::

    [workflow]
    file-retention-level = all_files
    h1-channel-name = H1:GDS-CALIB_STRAIN
    l1-channel-name = L1:GDS-CALIB_STRAIN

We may define all the detectors (IFOs) that we are considering in our analysis::

    [workflow-ifos]
    ; This is the list of ifos to analyse
    h1 =
    l1 =

The data frame types and other related options are given in::

    [workflow-datafind]
    datafind-method = AT_RUNTIME_SINGLE_CACHES
    datafind-check-segment-gaps = raise_error
    datafind-check-frames-exist = raise_error
    datafind-check-segment-summary = no_test
    datafind-h1-frame-type = H1_HOFT_C00
    datafind-l1-frame-type = L1_HOFT_C00

Data segment information is given in the section::

    [workflow-segments]
    segments-method = AT_RUNTIME
    segments-h1-science-name = H1:DMT-ANALYSIS_READY:1
    segments-l1-science-name = L1:DMT-ANALYSIS_READY:1
    segments-database-url = https://segments.ligo.org
    segments-veto-categories = 3
    segments-minimum-segment-length = 256
    segments-veto-definer-url = https://code.pycbc.phy.syr.edu/detchar/veto-definitions/download/master/cbc/ER8/H1L1-HOFT_C00_ER8B_CBC.xml
    
The GRB search requires an additional set of segment-related options,
which we give in the following section::

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

-------------------
Executable Sections
-------------------

We set the executables to be used for the analysis in the following way::
    
    [executables]
    inspiral                = ${which:lalapps_coh_PTF_inspiral}
    splitbank               = ${which:pycbc_splitbank}
    segment_query           = ${which:ligolw_segment_query_dqsegdb}
    segments_from_cats      = ${which:ligolw_segments_from_cats_dqsegdb}
    llwadd                  = ${which:ligolw_add}
    ligolw_combine_segments = ${which:ligolw_combine_segments}
    injections              = ${which:lalapps_inspinj}
    jitter_skyloc           = ${which:ligolw_cbc_jitter_skyloc}
    align_total_spin        = ${which:ligolw_cbc_align_total_spin}
    split_inspinj           = ${which:pycbc_split_inspinj}
    em_bright_filter        = ${which:pycbc_dark_vs_bright_injections}
    trig_combiner           = ${which:pylal_cbc_cohptf_trig_combiner}
    trig_cluster            = ${which:pylal_cbc_cohptf_trig_cluster}
    injfinder               = ${which:pylal_cbc_cohptf_injfinder}
    injcombiner             = ${Which:pylal_cbc_cohptf_injcombiner}
    sbv_plotter             = ${which:pylal_cbc_cohptf_sbv_plotter}
    efficiency              = ${which:pylal_cbc_cohptf_efficiency}
    inj_efficiency          = ${which:pylal_cbc_cohptf_efficiency}
    horizon_dist            = ${which:pylal_cbc_cohptf_inspiral_horizon}
    html_summary            = ${which:pycbc_make_grb_summary_page}

Here we are getting the executable paths from our environment for flexibility,
rather than supplying them as fixed paths.

The options to be given to every job run by an executable are then given
within a secion with the relevant name, for example our ``inspiral`` jobs (in
this case, lalapps_coh_PTF_inspiral) use the options in the following section::

    [inspiral]
    ligo-calibrated-data = real_8
    approximant = SpinTaylorT4
    order = threePointFivePN
    .
    .
    .

If the workflow were to contain multiple subclasses of ``inspiral`` jobs --
for example one for standard signal hunting and some for finding injected
signals -- options could be provided separately to these subclasses in tagged
sections. If the injection jobs are tagged in the workflow by the string
``coherent_injections``, then options specific to these jobs may be given in
the section::

    [inspiral-coherent_injections]
    inj-search-window = 1
    inj-mchirp-window = 0.05
    analyze-inj-segs-only =

Sections which share a common set of options may be given together::

    [inspiral&workflow-exttrig_segments]
    pad-data = 8

Here the ``workflow-exttrig_segments`` section and the ``inspiral`` executable
section are sharing a common option.


------------------------
Pegasus Profile Sections
------------------------

If, for example, we wished to ask condor to request nodes with 2000M of memory
for the ``trig_combiner`` executable jobs, we may do this via::

    [pegasus_profile-trig_combiner]
    condor|request_memory=2000M

This can be generalised to any executable or tagged jobs.

.. _howtorunpygrb:

==========
How to run
==========

Here we document the stages needed to run the triggered coherent GRB search.

Once PyCBC is installed, you should be able to run the following help command for the workflow generation
script::

    pycbc_make_offline_grb_workflow --help

This should produce a help message like the following

.. command-output:: pycbc_make_offline_grb_workflow --help

This outlines the command line arguments that may be passed to the executable.
The majority of options passed to the workflow will come from configuration
files, and these are known to the executable via the option
``--config-files``.

----------------------
Set up a run directory
----------------------

Navigate to the directory you wish to run in::

    RUN_DIR=/path/to/run/directory
    mkdir -p $RUN_DIR
    cd $RUN_DIR

Next gather together configuration files for your run.

---------------------------------------------------------------------------------
Configuration files - Are you running from production configuration (.ini) files?
---------------------------------------------------------------------------------
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Yes, I want to run in a standard production configuration
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

The option ``--config-files`` takes a space separated list of files locations.
These can be URLs to remote file locations. Production configuration files may
be found here_ (LIGO.ORG protected).

.. _here: https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-config

Therefore, an example run on a GRB from the 8th Advance LIGO engineering run
might use the following config files::

    pycbc_make_offline_grb_workflow \
    --config-files \
    https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-config/download/master/ER8/pygrb/analysis_er8.ini \
    https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-config/download/master/ER8/pygrb/injections_er8.ini \
    https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-config/download/master/ER8/pygrb/postprocessing_er8.ini \
    https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-config/download/master/ER8/pygrb/data_er8b.ini \
    https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-config/download/master/ER8/pygrb/offline_er8.ini \

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
No, I have my own configuration files
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

The option ``--config-files`` takes a space separated list of files locations.
For example, you could provide a pair of local files::

    pycbc_make_offline_grb_workflow \
    --config-files \
    /path/to/config_file_1.ini \
    /path/to/config_file_2.ini

Now go down to :ref:`pygrbgenerate`.

.. _pygrbgenerate:

=====================
Generate the workflow
=====================

When you are ready, you can generate the workflow. As this is a triggered
gravitational wave search, a number of key pieces of information will change
between one GRB and the next, such as the time of the GRB, or its position on
the sky. This may perhaps be most easily done by setting a number of variables
in your environment before launching the generation script.

First we need to set the trigger time, ie. the GPS Earth-crossing time of the
GRB signal. You should also set the GRB name. For example::

    GRB_TIME=1125614344
    GRB_NAME=150906B

We should next set the sky coordinates of the GRB in RA and Dec, in this
example::

    RA=159.239
    DEC=-25.603
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
    --config-files \
    https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-config/download/master/ER8/pygrb/analysis_er8.ini \
    https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-config/download/master/ER8/pygrb/injections_er8.ini \
    https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-config/download/master/ER8/pygrb/postprocessing_er8.ini \
    https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-config/download/master/ER8/pygrb/data_er8a.ini \
    https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-config/download/master/ER8/pygrb/offline_er8.ini \
    --config-overrides \
    workflow:ra:${RA} \
    workflow:dec:${DEC} \
    workflow:sky-error:${SKY_ERROR} \
    workflow:trigger-name:${GRB_NAME} \
    workflow:trigger-time:${GRB_TIME} \
    workflow:start-time:$(( GRB_TIME - 4096 )) \
    workflow:end-time:$(( GRB_TIME + 4096 )) \
    workflow:html-dir:${HTML_DIR}

.. _pygrbplan:

====================================
Planning and Submitting the Workflow
====================================

Change directory into the directory where the dax was generated::

    cd GRB${GRB_NAME}

From the directory where the dax was created, run the submission script::

    pycbc_submit_dax --dax pygrb_offline.dax --accounting-group <your.accounting.group.tag>

.. note::

    If running on the ARCCA cluster, please provide a suitable directory via
    the option --local-dir, ie. /var/tmp/${USER}

-------------------------------------------------------------------------------------------------------------------------------------------
Monitor and Debug the Workflow (`Detailed Pegasus Documentation <https://pegasus.isi.edu/wms/docs/latest/tutorial.php#idm78622034400>`_)
-------------------------------------------------------------------------------------------------------------------------------------------

To monitor the above workflow, one can run::

    pegasus-status -cl /path/to/analysis/run
    
To get debugging information in the case of failures.::

    pegasus-analyzer /path/to/analysis/run

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

======================================
Reuse of data from a previous workflow
======================================

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

