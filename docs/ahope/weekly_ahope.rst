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

-----------------------------
Edit the configuration file
-----------------------------

The configuration file for weekly_ahope is split into parts::

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

The remaining options affect how the jobs run, these should not be edited unless you know what you are doing ... but can freely be added if you do know what you are doing and want to change something.

-----------------------
Generate the workflow
-----------------------

When you are ready, you can generate the workflow. First we need to choose a time span::

    export GPS_START_TIME=961585543
    export GPS_END_TIME=961671943

You also need to specify the directory in which pipedown  will store log files. Ahope does not need this, but pipedown does.

 * For CIT,LHO,LLO or SYR set::

    export LOGPATH=/usr1/${USER}/log
    export PIPEDOWNLOG=/usr1/${USER}
    mkdir -p $LOGPATH

 * For Atlas set::

    export LOGPATH=/local/user/${USER}/log/
    export PIPEDOWNLOG=/local/user/${USER}
    mkdir -p $LOGPATH 

 * For UWM set::

    export LOGPATH=/people/${USER}/log/
    export PIPEDOWNLOG=/localscratch/${USER}
    mkdir -p $LOGPATH

You also need to choose where the html page will be generated. For example::

    export HTMLDIR=/home/${USER}/public_html/ahope

Then you can generate the workflow::

    python weekly_ahope.py --config-files weekly_ahope.ini pipedown.ini inj.ini \
                           --config-overrides ahope:start-time:${GPS_START_TIME} \
                                              ahope:end-time:${GPS_END_TIME} \
                                              ahope:ahope-html-basedir:${HTMLDIR} \
                                              ahope:pipedown-log-path:${LOGPATH} \
                                              ahope:pipedown-tmp-space:${PIPEDOWNLOG}

-----------------------------------------
Planning and Submitting the Worklfow
-----------------------------------------
First, copy the files needed for planning into the directory where the dax 
was generated.

    cp plan.sh ${GPS_START_TIME}-${GPS_END_TIME}/
    cp site-local.xml ${GPS_START_TIME}-${GPS_END_TIME}/
    cp pegasus.conf ${GPS_START_TIME}-${GPS_END_TIME}/

Then CD into the directory where the dax was generated::

    cd ${GPS_START_TIME}-${GPS_END_TIME}

From the directory where the dax was created, run the planning script::

    sh plan.sh weekly_ahope.dax
    
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

If the workflow runs successfully you will find the output under your html directory some time later.

-----------------------------------------
Monitor and Debug the Worklfow
-----------------------------------------

To monitor the above workflow, one would run::

    pegasus-status /usr1/ahnitz/log/ahnitz/pegasus/weekly_ahope/run0011
    
To get debugging information in the case of failures.::

    pegasus-analyzer /usr1/ahnitz/log/ahnitz/pegasus/weekly_ahope/run0011


