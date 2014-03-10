########################################################
Weekly ahope: A CBC analysis workflow generator
########################################################

===============
Introduction
===============

Weekly ahope is a tool used to analyse data from multiple detectors independently and then perform a coincidence test and various signal-based veto cuts and data quality cuts to determine whether or not a compact binary coalescence is present in the given data.

The output is a webpage containing the plots that can be used to understand the results of the analysis

=======================
How to run daily ahope
=======================

Here we document the stages needed to run daily ahope.

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

The scripts to run daily ahope currently reside within the pycbc source tree.
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

These files contain all the details needed to run daily_ahope::

    IF YOU ARE UNFAMILIAR WITH AHOPE WORKFLOWS, LOOK THROUGH THESE FILES. pipedown.ini WILL LOOK FAMILIAR IF YOU ARE USED TO IHOPE WORKFLOWS

* weekly_ahope.ini contains options that are used when running the ahope parts of the workflow
* pipedown.ini contains options that are used when running pipedown
* inj.ini contains the parameters used when generating simulation files

Some things that *will* need to be changed for each user in the weekly_ahope.ini file::

    [ahope]
    ahope-html-basedir = /home/spxiwh/public_html/ahope/development/weekly_ahope/test

    [executables]
    ; setup of condor universe and location of executables
    some_exe_name         = /path/to/executable
    some_other_exe        = /path/to/that/one/too

and some things that *will* need to be changed for each user in the pipedown.ini file::

    [executables]
    ; setup of condor universe and location of executables
    some_exe_name         = /path/to/executable
    some_other_exe        = /path/to/that/one/too

    [pipeline]
    node-tmp-dir = /usr/spxiwh


To run through this. 

 * The ahope-html-basedir is the directory in which you want the output html page to appear. The page will be in a subdirectory in this corresponding to the gps times.
 * Everything under [executables] in both .ini files points to the executables that will be used. These should be changed as appropriate
 * node-tmp-dir refers to a directory that must, a) exist, and b) point to a local disk on every node in the cluster you are running on.
   * /usr/${USER} is approriate for CIT, LLO, LHO and SYR.
   * /localscratch/${USER} is appropriate on UWM
   * /local/user/${USER} is appropriate on Atlas.

The example is also set up to run on S6 data and analysing only H1 and L1. If you are running on non-S6 data, or want to analyse a different set of ifos, you will have to edit some additional options::

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

To run through this::

 * The [ahope-ifos] section supplies which ifos will be analysed if data is found and available.
 * The X1-channel-name options are the h(t) channel name in the frames
 * The datafind-X1-frame-type is the type of the frames for use when calling gw_data_find
 * The segments-X1-science-name is the flag used to store science times in the segment database
 * segments-database-url points to the segment database
 * segments-veto-definer-url points to the url where the veto-definer file can be found.
 * The [tisi] sections give instructions on how to set up the input files used for time sliding. See :ref:`ahopetimeslidesmod` for more details on how to supply this for other situations. Some significant changes are needed to make ligolw_tisi and ahope work nicer with each other (and avoid this hassle). This is on the to-do list.

The remaining options affect how the jobs run, these should not be edited unless you know what you are doing ... but can freely be added if you do know what you are doing and want to change something.

-----------------------
Generate the workflow
-----------------------

When you are ready, you can generate the workflow. First we need to choose a time span::

    GPS_START_TIME=961585543
    GPS_END_TIME=961671943

You also need to specify the directory in which pipedown will store log files. Ahope does not need this, but pipedown does.

 * For CIT,LHO,LLO or SYR set::

    export LOGPATH=/usr1/${USER}/log
    mkdir -p $LOGPATH

 * For Atlas set::

    export LOGPATH=/local/user/${USER}/log/
    mkdir -p $LOGPATH 

 * For UWM set::

    export LOGPATH=/people/${USER}/log/
    mkdir -p $LOGPATH

Then you can generate the workflow::

    python weekly_ahope.py --config-files weekly_ahope.ini pipedown.ini inj.ini --config-overrides ahope:start-time:${GPS_START_TIME} ahope:end-time:${GPS_END_TIME} --pipedown-log-dir ${LOGPATH}

Then CD into the directory where the dag was generated::

    cd ${GPS_START_TIME}-${GPS_END_TIME}

where the directory naming is constructed from the year, month and day that is being analysed. Then submit the dag::

    condor_submit_dag weekly_ahope.dag

If the dag runs successfully you will find the output under your html directory some time later.

---------------------
Monitor the dagman
---------------------

One can follow the process of the dagman by running::

    tail -f weekly_ahope.dag.dagman.out

in the run directory to watch the progress of the dag. If jobs fail you should
look in the::

    SUBDIRECTORIES/logs

directory to see all the stderr and stdout files from each job. You can match these files with the condor process numbers given in the dagman.out to figure out which file corresponds to the failing jobs. You can also use::

    weekly_ahope.sh

to find the command line for each job in the ahope dag if you want to run by hand to debug any job.

We will soon be transitioning to pegasus which will make some of this easier!
