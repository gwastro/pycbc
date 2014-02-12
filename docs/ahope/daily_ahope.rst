########################################################
Daily ahope: A single-ifo detchar tool
########################################################

===============
Introduction
===============

Daily ahope is a tool used to analyse data from single detector(s), with no
coincidence stage. The output from this single detector matched-filter runs
can be used to identify times where the detectors are producing glitches
that have a large single-detector detection statistic for the CBC searches.
If the detectoralists can identify what causes these most egregious glitches
it can increase the overall search sensitivity.

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

    examples/ahope/er_daily_ahope

CD to this directory::

    cd ${SRC_DIR}/examples/ahope/er_daily_ahope

-----------------------------
Edit the configuration file
-----------------------------

The configuration file::

    daily_ahope.ini

contains all the details needed to run daily_ahope::

    LOOK THROUGH THIS FILE

Some things that *will* need to be changed for each user::

    [ahope]
    ahope-asset-dir = /home/spxiwh/lscsoft_git/src/lalsuite/lalapps/src/inspiral
    ahope-html-basedir = /home/spxiwh/public_html/ER4/test

    [ahope-segments]
    segments-veto-definer-file = /home/spxiwh/lscsoft_git/src/pycbc/pycbc/examples/ahope/er_daily_ahope/H1L1V1-ER3_CBC_OFFLINE-1011571215-0.xml

    [ahope-omega]
    omega-conf-file = /home/spxiwh/ERs/ER4/daily_ihope_test/old_conf_omega.txt

    [executables]
    ; setup of condor universe and location of executables
    tmpltbank         = /usr/bin/lalapps_tmpltbank
    inspiral          = /usr/bin/lalapps_inspiral
    splittable = /usr/bin/lalapps_splitbank
    segment_query = /usr/bin/ligolw_segment_query
    segments_from_cats = /home/spxiwh/lscsoft_git/executables_master/bin/ligolw_segments_from_cats
    ligolw_add = /usr/bin/ligolw_add
    siclustercoarse = /home/spxiwh/lscsoft_git/executables_master/bin/ligolw_sicluster
    siclusterfine = /home/spxiwh/lscsoft_git/executables_master/bin/ligolw_sicluster
    ihope_daily_page = /home/spxiwh/ERs/ER4/daily_ihope_test/lalapps_ihope_daily_page
    cbc_glitch_page = /home/spxiwh/lscsoft_git/executables_master/bin/ligolw_cbc_glitch_page
    cbc_hardware_inj_page = /home/spxiwh/lscsoft_git/executables_master/bin/ligolw_cbc_hardware_inj_page

To run through this. 

 * The ahope-asset-dir points to the location where the CSS files needed for the html page are stored. This is your lalsuite source directory and then /lalapps/src/inspiral. This will be removed in the future and read from a web-accessible location.
 * The ahope-html-basedir is the directory in which you want the output html pages to appear. Pages will be in subdirectories in this corresponding to unique days. So for me I might see output in directories like /home/spxiwh/public_html/ER4/test/201308/20130812/
 * The segments-veto-definer-file points to the veto-definer-file. This will be changed to take a url input not a filepath.
 * omega-conf-file points to the configuration for omega. Currently omega doesn't work within daily ahope, this will be fixed.
 * Everything under [executables] points to the executables that will be used. These should be changed as appropriate

The example is also set up to run on ER4 data. If you are running on non-ER4 data you may have to edit some additional options, for e.g.::

    [ahope]
    h1-channel-name = H1:FAKE-STRAIN
    l1-channel-name = L1:FAKE-STRAIN
    v1-channel-name = V1:FAKE_h_16384Hz_4R

    [ahope-datafind]
    datafind-h1-frame-type = H1_ER_C00_L1
    datafind-l1-frame-type = L1_ER_C00_L1
    datafind-v1-frame-type = V1Online

    [ahope-segments]
    segments-H1-science-name = H1:DMT-SCIENCE:1
    segments-L1-science-name = L1:DMT-SCIENCE:1
    segments-V1-science-name = V1:ITF_SCIENCEMODE
    segments-database-url = https://segdb-er.ligo.caltech.edu
    segments-veto-definer-url = https://www.lsc-group.phys.uwm.edu/ligovirgo/cbc/public/segments/ER4/H1L1V1-ER4_CBC_OFFLINE-1011571215-0.xml

    [ahope-omega]
    omega-frame-dir = /frames/ER4/L1_ER_C00_L1/L1/L-L1_ER_C00_L1-%%d/L-L1_ER_C00_L1-%%d

To run through this::

 * The X1-channel-name options are the h(t) channel name in the frames
 * The datafind-X1-frame-type is the type of the frames for use when calling gw_data_find
 * The segments-X1-science-name is the flag used to store science times in the segment database
 * segments-database-url points to the segment database
 * segments-veto-definer-url points to the url where the veto-definer file can be found.

The remaining options affect how the jobs run, these should not be edited unless you know what you are doing!

-----------------------
Generate the workflow
-----------------------

When you are ready, you can generate the workflow. First we need to choose a time::

    GPS_START_TIME=1059436816

This time should be a gps time during the *same day* that you want to analyse. Daily ahope will analyse from 00:00:00 to 23:59:59 of that day. Then you can generate the workflow::

    python daily_ahope.py -s ${GPS_START_TIME} -i daily_ahope.ini -d ${PWD}

Then CD into the directory where the dag was generated::

    cd 201308/20130802

where the directory naming is constructed from the year, month and day that is being analysed. Then submit the dag::

    condor_submit_dag daily_ahope.dag

If the dag runs successfully you will find the output under your html directory some time later.

---------------------
Monitor the dagman
---------------------

One can follow the process of the dagman by running::

    tail -f daily_ahope.dag.dagman.out

in the run directory to watch the progress of the dag. If jobs fail you should
look in the::

    logs

directory to see all the stderr and stdout files from each job. You can match these files with the condor process numbers given in the dagman.out to figure out which file corresponds to the failing jobs. You can also use::

    daily_ahope.sh

to find the command line for each job in the ahope dag if you want to run by hand to debug any job.
