############################################################
``pycbc_make_coinc_workflow``: A CBC analysis workflow generator
############################################################

===============
Introduction
===============

The executable ``pycbc_make_coinc_workflow`` is a tool used to analyse data from multiple detectors independently and then perform a coincidence test and various signal-based veto cuts and data quality cuts to determine whether or not a compact binary coalescence is present in the given data.

The output is a webpage containing the plots that can be used to understand the results of the analysis

.. _howtorunworkflow:

=======================
How to run
=======================

Here we document the stages needed to run ``pycbc_make_coinc_workflow`` to generate an offline matched filtered CBC search.

---------------------------
Install lalsuite and pycbc
---------------------------

The first thing that is needed is to install lalsuite and pycbc. This is
described on the page here:

.. toctree::
   :maxdepth: 1

   ../install

----------------------------------------------------------------------------
The configuration file - Do you already have configuration (.ini) file(s)?
----------------------------------------------------------------------------
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Yes, I already have my own local config files
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

Great! Then copy the configuration files into your run directory::

    cp /path/to/config_file1.ini /path/to/config_file2.ini .

and set the names of these configuration files in your path. If you have more than one configuration file they must be space separated::

    LOCAL_CONFIG_FILES="config_file1.ini config_file2.ini"

Now go down to :ref:`coincworkflowgenerate`.

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Yes, I would like to use the unmodified preinstalled configuration files
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

For a full list of the preinstalled configuration files, see :ref:`configuration_files`.

Set the configurations files in your path and proceed to workflow generation::

    INSTALLED_CONFIG_FILES="example_pycbc.ini example_inj.ini example_pipedown.ini"

Now go down to :ref:`coincworkflowgenerate`.

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
No, I need to make a configuration file - Editing the example files
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

An example configuration file set is found in three parts::

    /src/dir/pycbc/workflow/ini_files/example_pycbc.ini
    /src/dir/pycbc/workflow/ini_files/example_pipedown.ini
    /src/dir/pycbc/workflow/ini_files/example_inj.ini

These files contain all the details needed to run ``pycbc_make_coinc_workflow``

.. note::

    If you are unfamiliar with pycbc workflows, look through these files.
    example_pipedown.ini will look familiar if you are used to ihope workflows.

* example_pycbc.ini contains options that are used when running the pycbc.workflow parts of the workflow
* example_pipedown.ini contains options that are used when running pipedown
* example_inj.ini contains the parameters used when generating simulation files

Alternatively, if you want to run with lalapps executables replace example_pycbc.ini with::

   example_lalapps.ini

If you want to run in this default configuration please jump down the "Generate the workflow".

If you want to run on non-S6 data, or want to analyse a different set of ifos, you will have to edit some additional options::

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
    segments-H1-science-name = H1:DMT-SCIENCE:4
    segments-L1-science-name = L1:DMT-SCIENCE:4
    segments-V1-science-name = V1:ITF_SCIENCEMODE
    segments-database-url = https://segdb.ligo.caltech.edu
    segments-veto-definer-url = https://www.lsc-group.phys.uwm.edu/ligovirgo/cbc/public/segments/S6/H1L1V1-S6_CBC_LOWMASS_B_OFFLINE-937473702-0.xml

    ALL the [tisi], [tisi-zerolag], [tisi-slides] sections (potentially)

To run through this

* The ``[workflow-ifos]`` section supplies which ifos will be analysed if data is found and available.
* The ``X1-channel-name`` options are the h(t) channel name in the frames
* The ``datafind-X1-frame-type`` is the type of the frames for use when calling gw_data_find
* The ``segments-X1-science-name`` is the flag used to store science times in the segment database
* The ``segments-database-url`` points to the segment database
* The ``segments-veto-definer-url`` points to the url where the veto-definer file can be found.
* The ``[tisi]`` sections give instructions to the workflow module on how to set up what time slides will be performed. See :ref:`workflowtimeslidesmod` for more details on how to supply this for other situations. Normally you will just need to add or remove detectors.

The remaining options affect how the jobs run, these should not be edited unless you know what you are doing ... but can freely be added if you do know what you are doing and want to change something. To find out more details about the possible options for any stage of the workflow, follow the links at :ref:`workflowhomepage`.

Now you have configuration files and can follow the same instructions as above. That is: 

Copy the configuration files into your run directory::

    cp /path/to/weekly_ahope.ini /path/to/inj.ini /path/to/pipedown.ini .

and set the names of these configuration files in your path. If you have more than one configuration file they must be space separated::

    LOCAL_CONFIG_FILES="weekly_ahope.ini inj.ini pipedown.ini"

.. _coincworkflowgenerate:

-----------------------
Generate the workflow
-----------------------

When you are ready, you can generate the workflow. First we need to choose a time span. Here is an example::

    export GPS_START_TIME=967593543
    export GPS_END_TIME=967679943

You also need to specify the directory for storing log files.

 * For CIT,LHO,LLO or SYR set::

    export LOGPATH=/usr1/${USER}/log
    export PIPEDOWNTMPSPACE=/usr1/${USER}
    mkdir -p $LOGPATH

 * For Atlas set::

    export LOGPATH=/local/user/${USER}/log/
    export PIPEDOWNTMPSPACE=/local/user/${USER}
    mkdir -p $LOGPATH 

 * For UWM set::

    export LOGPATH=/people/${USER}/log/
    export PIPEDOWNTMPSPACE=/localscratch/${USER}
    mkdir -p $LOGPATH

 * On the TACC XSEDE cluster, it is recommended to store your ihope directory under the work filesystem.
   For the TACC XSEDE cluster set::

    export LIGO_DATAFIND_SERVER=tacc.ligo.org:80
    export LOGPATH=${SCRATCH}/log
    export PIPEDOWNTMPSPACE=/tmp
    mkdir -p $LOGPATH

You also need to choose where the html page will be generated. For example::

    export HTMLDIR=/home/${USER}/public_html/ahope

If you are using locally editted or custom configuration files then you can
create the workflow using::

    pycbc_make_coinc_workflow --local-config-files ${LOCAL_CONFIG_FILES} \
                              --config-overrides workflow:start-time:${GPS_START_TIME} \
                                                 workflow:end-time:${GPS_END_TIME} \
                                                 workflow:workflow-html-basedir:${HTMLDIR} \
                                                 workflow:pipedown-log-path:${LOGPATH} \
                                                 workflow:pipedown-tmp-space:${PIPEDOWNTMPSPACE}
                                              
                                              
If you are using default installed configuration files then you can create the
workflow using::

    pycbc_make_coinc_workflow --installed-config-files ${INSTALLED_CONFIG_FILES} \
                              --config-overrides workflow:start-time:${GPS_START_TIME} \
                                                 workflow:end-time:${GPS_END_TIME} \
                                                 workflow:workflow-html-basedir:${HTMLDIR} \
                                                 workflow:pipedown-log-path:${LOGPATH} \
                                                 workflow:pipedown-tmp-space:${PIPEDOWNTMPSPACE}

.. _coincworkflowplan:

-----------------------------------------
Planning and Submitting the Workflow
-----------------------------------------
CD into the directory where the dax was generated::

    cd ${GPS_START_TIME}-${GPS_END_TIME}

From the directory where the dax was created, run the planning script::

    pycbc_basic_pegasus_plan weekly_ahope.dax $LOGPATH
    
Submit the workflow by following the instructions at the end of the script output, which looks something like 
the following.::

    ...
    10:49:15:INFO : Entering post-processing preperation stage.
    10:49:15:INFO : Leaving post-processing separation module.
    10:49:18:INFO : Finished.
    2014.03.26 10:49:28.676 EDT:   


    I have concretized your abstract workflow. The workflow has been entered 
    into the workflow database with a state of "planned". The next step is 
    to start or execute your workflow. The invocation required is


    pegasus-run  /usr1/ahnitz/log/ahnitz/pegasus/weekly_ahope/run0011

     
    2014.03.26 10:49:28.983 EDT:   Time taken to execute is 7.095 seconds 
    
In this case, the workflow would be submitted as follows.::

    pegasus-run  /usr1/ahnitz/log/ahnitz/pegasus/weekly_ahope/run0011

If the workflow runs successfully, you will find the output under your html directory some time later.

-------------------------------------------------------------------------------------------------------------------------------------------
Monitor and Debug the Workflow (`Detailed Pegasus Documentation <https://pegasus.isi.edu/wms/docs/latest/tutorial.php#idm78622034400>`_)
-------------------------------------------------------------------------------------------------------------------------------------------

To monitor the above workflow, one would run::

    pegasus-status /usr1/ahnitz/log/ahnitz/pegasus/weekly_ahope/run0011
    
To get debugging information in the case of failures.::

    pegasus-analyzer /usr1/ahnitz/log/ahnitz/pegasus/weekly_ahope/run0011

=======================
Post-processing
=======================

-----------------------------------------
Summary page
-----------------------------------------

A summary page will be created at the end of your weekly workflow. The directory is specified by the evironment varaible HTMLDIR that was set when you ran weekly_ahope.py to generate the workflow. For example::

    /home/${USER}/public_html/workflow

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Full data summary
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

A summary of the results are displayed in Section 3.

Section 3.1 displays the search sensitivity of the detectors over time. The inspiral horizion distance is the distance an optimally-oriented equal-mass system would give SNR equal to 8.

Section 3.2 to 3.4 contain cumulative histograms of coincident triggers against IFAR (inverse false-alarm rate). The blue triangles are coincident triggers and a table of the loudest events is provided below the plot. We use time slides to find background triggers and calculate the false-alarm rate as the number of triggers louder than a given trigger divided by the total background time. If a trigger is louder than all background triggers then we set its false-alarm rate to 0. A low false-alarm rate gives a high IFAR so we plot the trigger with an arrow pointing to the right. This indicates that its true IFAR is somewhere to the right.

Section 3.2 includes hardware injections, so note that a number of these signals are recovered with 0 combined false-alarm rate. These would be detection candidates if they were not hardware injections.

Section 3.4 removes hardware injections and times marked by CAT_3 vetoes.

Section 3.5 shows the recovery of the simulated signals that were added in this workflow.

=============================
Workflow visualization
=============================

-----------------------------
Pegasus Dashboard
-----------------------------

The `pegeasus dashboard <http://pegasus.isi.edu/wms/docs/latest/ch02s11.php>`_ is a visual and interactive way to get information about the progress, status, etc of your workflows. 

The software can be obtained from a seprate pegasus package here <https://github.com/pegasus-isi/pegasus-service>. 

Pegasus Dashboard is currently installed on sugar. To view your Pegasus Dashboard, in a browser go to::

    https://sugar.phy.syr.edu/pegasus/~albert.einstein

You will be greet with a page that has a table of all your workflows that were submitted from sugar. You can view the details of a workflow by clicking on the link in the Workflow Details column of the table.

Clicking on the Workflow Details link will take you to a webpage that gives a high-level overview of the workflow, telling you how many many jobs succeeded, fail, the submit directory, etc. There is a table with tabs at the bottom of the page. If you click the tabs Failed, Running, and Successful the page will generate a table that lists all the failed, running, and successful jobs in the workflow respectively. You also have the ability to search the table for a particular kind of job using the Search bar.

You can view the details of a job by clicking the link in the Job Name column. This will take you to a Job Details page. This page will tell you where to find stdout files, stderr files, how much wall clock time the job took to run, etc. There is a table at the bottom of the page with a Failed and Successful tab. If you click on the respective tab, it will list all invocations of that job. You can click on the link in the Invocations column for more information.

On the Invocation Details page there is information about the command line arguments, executable path, CPU time, wall clock time, etc.

-----------------------------
Pegasus Analyzer
-----------------------------

The `pegeasus analyzer <http://pegasus.isi.edu/wms/docs/trunk/cli-pegasus-analyzer.php>`_ is a command-line tool for reporting sucessful and failed jobs.

To run ``pegasus_analyzer`` on your workflow, type::

    pegasus-analyzer /usr1/ahnitz/log/ahnitz/pegasus/weekly_ahope/run0011

``pegasus_analyzer`` will display a summary of suceeded, failed, and unsubmitted jobs in the workflow. After the summary information, ``pegasus_analyzer`` will display information about each failed job. An example would be::

    ************************************Summary*************************************

    Submit Directory   : /usr1/cbiwer/log/H1L1V1-s6d_test-970012743-258000.9apn7X
    Total jobs         :     24 (100.00%)
    # jobs succeeded   :     19 (79.17%)
    # jobs failed      :      5 (20.83%)
    # jobs unsubmitted :      0 (0.00%)

    ******************************Failed jobs' details******************************

    =====================ligolw_cbc_hardware_inj_page_ID000020======================

    last state: POST_SCRIPT_FAILED
         site: local
    submit file: ligolw_cbc_hardware_inj_page_ID000020.sub
    output file: ligolw_cbc_hardware_inj_page_ID000020.out.001
    error file: ligolw_cbc_hardware_inj_page_ID000020.err.001

    -------------------------------Task #1 - Summary--------------------------------

    site        : local
    hostname    : avhe2010.sugar.phy.syr.edu
    executable  : /home/cbiwer/projects/test_workflow/970012743-970270743/executables/ligolw_cbc_hardware_inj_page
    arguments   : --source-xml hardware_injection_summary/H1L1V1-S6_CBC_HW_INJECTIONS-930493015-42111800.xml --outfile hardware_injection_summary/H1L1V1-HWINJ_SUMMARY_CAT_2-9
    70012743-258000.html ----segments-xml-glob ../segments/*-SCIENCE_SEGMENTS-*-*.xml --v1-injections ----vetos-xml-glob ../segments/*-COMBINED_CAT_2_VETO_SEGS-*-*.xml --gps-
    start-time 970012743 --segment-dir hardware_injection_summary --gps-end-time 970270743 --l1-injections --analyze-injections --cache-file full_data/H1L1V1-INSPIRAL_HIPE_FU
    LL_DATA_CAT_2_VETO-970012743-258000.cache --h1-injections --cache-pattern *SIRE_FIRST*
    exitcode    : 2
    working dir : /home/cbiwer/projects/test_workflow/970012743-970270743

    Task #1 - ligo-hwinjpagejob::ligolw_cbc_hardware_inj_page:1.0 - ID000020 - Kickstart stderr

    Usage:  ligolw_cbc_hardware_inj_page [options]
    Program to parse the inspiral injection log
    ligolw_cbc_hardware_inj_page: error: no such option: ----segments-xml-glob

The output provides you with the ``stderr``, the command line, and where the job was run.

If you have a subdax that failed, ``pegasus_analyzer`` will provide you with a command to recieve more information about the failed jobs in the subdax.

-----------------------------
Pegasus Plots
-----------------------------

Pegasus has a tool called pegasus-plan to visualize workflows. To generate these charts and create an summary html page with this information, one would run::

    export PPLOTSDIR=/home/ahnitz/public_html/workflow/pegasus_plots
    pegasus-plots --plotting-level all --output ${PPLOTSDIR} /usr1/ahnitz/log/ahnitz/pegasus/weekly_ahope/run0011

The Invocation Breakdown Chart section gives a snapshot of the workflow. You can click on the slices of the pie chart and it will report the number of failures, average runtime, and max/min runtime for that type of jobs in the workflow. The radio button labeled runtime will organize the pie chart by total runtime rather than the total number of jobs for each job type.

The Workflow Execution Gantt Chart section breaks down the workflow how long it took to run each job. You can click on a job in the gantt chart and it will report the job name and runtime.

The Host Over Time Chart section displays a gantt chart where you can see what jobs in the workflow ran on a given machine.

.. _weeklyahopereuse:

================================
Reuse of workflow file products
================================

One of the features of  Pegasus is to reuse the data products of prior runs.
This can be used to expand an analysis or recover a run with mistaken settings without
duplicating work.

-----------------------------------------
Generate the full workflow you want to do
-----------------------------------------

First generate the full workflow for the
run you would like to do as normal. Follow the instructions of this page from :ref:`howtorunworkflow`,
but stop before planning the workflow with plan.sh in :ref:`coincworkflowplan`.

-----------------------------------------------------
Select the files you want to reuse from the prior run
-----------------------------------------------------

Locate the directory of the run that you would like to reuse. There is a file
called ${GPS_START_TIME}-${GPS_END_TIME}/output.map. This file contains a 
listing of all of the data products of the prior workflow, and can be used to tell
pegasus to skip regenerating them.

Select the entries in the file that you would like to skip generating again and
place that into a new file. The example below selects all the inspiral and 
tmpltbank jobs and places their entries into a new listing called prior_data.map.::

    # Lets get the tmpltbank entries
    cat /path/to/old/run/${GPS_START_TIME}-${GPS_END_TIME}/output.map | grep 'TMPLTBANK' > prior_data.map
    
    # Add in the inspiral  files
    cat /path/to/old/run/${GPS_START_TIME}-${GPS_END_TIME}/output.map | grep 'INSPIRAL' >> prior_data.map

.. note::

    You can include files in the prior data listing that wouldn't be generated anyway by your new run. These are simply
    ignored.

---------------------------
Plan the workflow
---------------------------

From the directory where the dax was created, now plan the workflow with an additional argument as follows.::

    pycbc_basic_pegasus_plan weekly_ahope.dax $LOGPATH --cache /path/to/prior_data.map

Follow the remaining :ref:`coincworkflowplan` instructions to submit your reduced
workflow.
