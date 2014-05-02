########################################################
Weekly ahope: A CBC analysis workflow generator
########################################################

===============
Introduction
===============

Weekly ahope is a tool used to analyse data from multiple detectors independently and then perform a coincidence test and various signal-based veto cuts and data quality cuts to determine whether or not a compact binary coalescence is present in the given data.

The output is a webpage containing the plots that can be used to understand the results of the analysis

=======================
How to run weekly ahope
=======================

Here we document the stages needed to run weekly ahope.

---------------------------
Install lalsuite and pycbc
---------------------------

The first thing that is needed is to install lalsuite and pycbc. This is
described on the page here:

.. toctree::
   :maxdepth: 1

   ../install

----------------------
Find the run scripts
----------------------

The scripts to run weekly ahope currently reside within the pycbc source tree.
These will be moved to be installed executables at some point. For now this
can be found in::

    examples/ahope/weekly_ahope

cd to this directory::

    cd ${PYCBC_SRC_DIR}/examples/ahope/weekly_ahope

If you want to run in a different directory then you can copy the files to that directory::

    cp ${PYCBC_SRC_DIR}/examples/ahope/weekly_ahope/* /path/to/your/run/directory
    cd /path/to/your/run/directory

----------------------------------------------------------------------------
The configuration file - Do you already have configuration (.ini) file(s)?
----------------------------------------------------------------------------

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Yes, I already have configuration files
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

Great! Then copy the configuration files into your run directory (overwrite the example files if you need to)::

    cp /path/to/config_file1.ini /path/to/config_file2.ini .

and set the names of these configuration files in your path. If you have more than one configuration file they must be space separated::

    CONFIG_FILES="config_file1.ini config_file2.ini"

Now go down to :ref:`weeklyahopegenerate`.

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
No, I need to make a configuration file - Editing the example files
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

The default configuration file for weekly_ahope is found in three parts::

    weekly_ahope.ini
    pipedown.ini
    inj.ini

These files contain all the details needed to run weekly_ahope

.. note::

    If you are unfamiliar with ahope workflows, look through these files.
    pipedown.ini will look familiar if you are used to ihope workflows.

* weekly_ahope.ini contains options that are used when running the ahope parts of the workflow
* pipedown.ini contains options that are used when running pipedown
* inj.ini contains the parameters used when generating simulation files

Alternatively, if you want to run with pycbc executables replace weekly_ahope.ini with::

   weekly_ahope_pycbc.ini

The weekly_ahope.ini example is set up to run on S6 data and analysing only H1 and L1. The weekly_ahope_pycbc.ini example is set up to run on S6 data analysing H1, L1 and V1.

If you want to run in this default configuration please jump down the "Generate the workflow".

If you want to run on non-S6 data, or want to analyse a different set of ifos, you will have to edit some additional options::

    [ahope]
    h1-channel-name = H1:LDAS-STRAIN
    l1-channel-name = L1:LDAS-STRAIN

    [ahope-ifos]
    ; This is the list of ifos to analyse
    h1 =
    l1 =

    [ahope-datafind]
    datafind-h1-frame-type = H1_LDAS_C02_L2
    datafind-l1-frame-type = L1_LDAS_C02_L2

    [ahope-segments]
    segments-H1-science-name = H1:DMT-SCIENCE:4
    segments-L1-science-name = L1:DMT-SCIENCE:4
    segments-V1-science-name = V1:ITF_SCIENCEMODE
    segments-database-url = https://segdb.ligo.caltech.edu
    segments-veto-definer-url = https://www.lsc-group.phys.uwm.edu/ligovirgo/cbc/public/segments/S6/H1L1V1-S6_CBC_LOWMASS_B_OFFLINE-937473702-0.xml

    ALL the [tisi], [tisi-zerolag], [tisi-slides] sections (potentially)

To run through this

* The [ahope-ifos] section supplies which ifos will be analysed if data is found and available.
* The X1-channel-name options are the h(t) channel name in the frames
* The datafind-X1-frame-type is the type of the frames for use when calling gw_data_find
* The segments-X1-science-name is the flag used to store science times in the segment database
* segments-database-url points to the segment database
* segments-veto-definer-url points to the url where the veto-definer file can be found.
* The [tisi] sections give instructions to ahope on how to set up what time slides will be performed. See :ref:`ahopetimeslidesmod` for more details on how to supply this for other situations. Normally you will just need to add or remove detectors.

The remaining options affect how the jobs run, these should not be edited unless you know what you are doing ... but can freely be added if you do know what you are doing and want to change something. To find out more details about the possible options for any stage of the workflow, follow the links at :ref:`ahopehomepage`.

Now you have configuration files and can follow the same instructions as above. That is: 

Copy the configuration files into your run directory::

    cp /path/to/weekly_ahope.ini /path/to/inj.ini /path/to/pipedown.ini .

and set the names of these configuration files in your path. If you have more than one configuration file they must be space separated::

    CONFIG_FILES="weekly_ahope.ini inj.ini pipedown.ini"

.. _weeklyahopegenerate:

-----------------------
Generate the workflow
-----------------------

When you are ready, you can generate the workflow. First we need to choose a time span. Here is an example::

    export GPS_START_TIME=961585543
    export GPS_END_TIME=961671943

You also need to specify the directory in which pipedown  will store log files. Ahope does not need this, but pipedown does.

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

Then you can generate the workflow::

    python weekly_ahope.py --config-files ${CONFIG_FILES} \
                           --config-overrides ahope:start-time:${GPS_START_TIME} \
                                              ahope:end-time:${GPS_END_TIME} \
                                              ahope:ahope-html-basedir:${HTMLDIR} \
                                              ahope:pipedown-log-path:${LOGPATH} \
                                              ahope:pipedown-tmp-space:${PIPEDOWNTMPSPACE}

-----------------------------------------
Planning and Submitting the Worklfow
-----------------------------------------
First, copy the files needed for planning into the directory where the dax 
was generated.::

    cp plan.sh ${GPS_START_TIME}-${GPS_END_TIME}/
    cp site-local.xml ${GPS_START_TIME}-${GPS_END_TIME}/
    cp pegasus.conf ${GPS_START_TIME}-${GPS_END_TIME}/

Then CD into the directory where the dax was generated::

    cd ${GPS_START_TIME}-${GPS_END_TIME}

From the directory where the dax was created, run the planning script::

    sh plan.sh weekly_ahope.dax

.. note::    
   Note, that if you have changed the name of the ini file, or are using an alternative one,
   the dax file's name will change accordingly.
    
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

-----------------------------------------
Monitor and Debug the Workflow
-----------------------------------------

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

An ahope summary page will be created at the end of your weekly ahope workflow. The directory is specified by the evironment varaible HTMLDIR that was set when you ran weekly_ahope.py to generate the workflow. For example::

    /home/${USER}/public_html/ahope

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Full data summary
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

A summary of the results are displayed in Section 3.

Section 3.1 displays the search sensitivity of the detectors over time. The inspiral horizion distance is the distance an optimally-oriented equal-mass system would give SNR equal to 8.

Section 3.2 to 3.4 contain cumulative histograms of coincident triggers against IFAR (inverse false-alarm rate). The blue triangles are coincident triggers and a table of the loudest events is provided below the plot. We use time slides to find background triggers and calculate the false-alarm rate as the number of triggers louder than a given trigger divided by the total background time. If a trigger is louder than all background triggers then we set its false-alarm rate to 0. A low false-alarm rate gives a high IFAR so we plot the trigger with an arrow pointing to the right. This indicates that its true IFAR is somewhere to the right.

Section 3.2 includes hardware injections, so note that a number of these signals are recovered with 0 combined false-alarm rate. These would be detection candidates if they were not hardware injections.

Section 3.4 removes hardware injections and times marked by CAT_3 vetoes.

Section 3.5 shows the recovery of the simulated signals that were added in this workflow.

-----------------------------------------
Workflow visualization
-----------------------------------------

Pegasus has a tool called pegasus-plan to visualize workflows. To generate these charts and create an summary html page with this information, one would run::

    export PPLOTSDIR=/home/ahnitz/public_html/ahope/pegasus_plots
    pegasus-plots --plotting-level all --output ${PPLOTSDIR} /usr1/ahnitz/log/ahnitz/pegasus/weekly_ahope/run0011

The Invocation Breakdown Chart section gives a snapshot of the workflow. You can click on the slices of the pie chart and it will report the number of failures, average runtime, and max/min runtime for that type of jobs in the workflow. The radio button labeled runtime will organize the pie chart by total runtime rather than the total number of jobs for each job type.

The Workflow Execution Gantt Chart section breaks down the workflow how long it took to run each job. You can click on a job in the gantt chart and it will report the job name and runtime.

The Host Over Time Chart section displays a gantt chart where you can see what jobs in the workflow ran on a given machine.

