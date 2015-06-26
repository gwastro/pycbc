########################################################
``pycbc_make_daily_workflow``: A single-ifo detchar tool
########################################################

===============
Introduction
===============

The executable ``pycbc_make_daily_workflow`` is a tool used to analyse data
from single detector(s), with no coincidence stage. The output from this single
detector matched-filter runs can be used to identify times where the detectors
are producing glitches that have a large single-detector detection statistic
for the CBC searches. If the detectoralists can identify what causes these most
egregious glitches it can increase the overall search sensitivity.

=========================================
How to run ``pycbc_make_daily_workflow``
=========================================

Here we document the stages needed to run ``pycbc_make_daily_workflow``.

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

Now go down to :ref:`pycbcdailygenerate`.

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Yes, I would like to use the unmodified preinstalled configuration files
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

For a full list of the preinstalled configuration files, see :ref:`configuration_files`.

Set the configurations files in your path and proceed to workflow generation::

    INSTALLED_CONFIG_FILES="example_daily_lalapps.ini"

Now go down to :ref:`pycbcdailygenerate`.

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
No, I need to make a configuration file - Editing the example files
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

An example configuration file set is found here::

    /src/dir/pycbc/workflow/ini_files/example_daily_lalapps.ini

This file contains all the details needed to run ``pycbc_make_daily_workflow``

.. note::

    If you are unfamiliar with pycbc workflows, look through these files.
    example_daily_lalapps.ini will look familiar if you are used to daily ihope workflows.

The ``example_daily_lalapps.ini`` example is set up to run on S6 data and analysing only H1 and L1.

If you want to run in this default configuration please jump down the "Generate the workflow".

If you want to run on non-ER5 data, or want to analyse a different set of ifos, you will have to edit some additional options::

    [workflow]
    h1-channel-name = H1:FAKE-STRAIN
    l1-channel-name = L1:FAKE-STRAIN

    [workflow-ifos]
    ; This is the list of ifos to analyse
    h1 =
    l1 =

    [workflow-datafind]
    datafind-h1-frame-type = H1_ER_C00_L1
    datafind-l1-frame-type = L1_ER_C00_L1

    [workflow-segments]
    segments-H1-science-name = H1:DMT-SCIENCE:1
    segments-L1-science-name = L1:DMT-SCIENCE:1
    segments-database-url = https://segdb-er.ligo.caltech.edu
    segments-veto-definer-url = https://www.lsc-group.phys.uwm.edu/ligovirgo/cbc/public/segments/ER4/H1L1V1-ER4_CBC_OFFLINE-1011571215-0.xml

To run through this

* The ``[workflow-ifos]`` section supplies which ifos will be analysed if data is found and available.
* The ``X1-channel-name`` options are the h(t) channel name in the frames.
* The ``datafind-X1-frame-type`` options are the type of the frames for use when calling ``gw_data_find``.
* The ``segments-X1-science-name`` options are the flag used to store science times in the segment database.
* The ``segments-database-url`` option points to the segment database.
* The ``segments-veto-definer-url`` option points to the url where the veto-definer file can be found.

The remaining options affect how the jobs run, these should not be edited unless you know what you are doing. To find out more details about the possible options for any stage of the workflow, follow the links at :ref:`workflowhomepage`.

Now you have configuration files and can follow the same instructions as above. That is: 

Copy the configuration files into your run directory::

    cp /path/to/example_daily_lalapps.ini .

and set the names of these configuration files in your path. If you have more than one configuration file they must be space separated::

    LOCAL_CONFIG_FILES="example_daily_lalapps.ini"

.. _pycbcdailygenerate:

-----------------------
Generate the workflow
-----------------------

When you are ready, you can generate the workflow. First we need to choose a start time. Here is an example::

    export DAY=2/23/2014
    export GPS_START_TIME=`lalapps_tconvert ${DAY}`
    export MONTHDIR=`lalapps_tconvert -f "%Y%m" ${GPSSTART}`
    export DAYDIR=`lalapps_tconvert -f "%Y%m%d" ${GPSSTART}`

The executable ``pycbc_make_daily_workflow`` assumes that the end time will be 24 hours from the start time.

You also need to specify the directory for storing log files.

 * For CIT,LHO,LLO or SYR set::

    export LOGPATH=/usr1/${USER}/log
    mkdir -p $LOGPATH

 * For Atlas set::

    export LOGPATH=/local/user/${USER}/log/
    mkdir -p $LOGPATH 

 * For UWM set::

    export LOGPATH=/people/${USER}/log/
    mkdir -p $LOGPATH

You also need to choose where the html page will be generated. For example::

    export HTMLDIR=/home/${USER}/public_html/daily_cbc_offline

You also need to tell the post-processing dag the path of the css file. For example::

    export ASSETDIR=/src/dir/lalsuite/lalapps/src/inspiral

If you are using locally edited or custom configuration files then you can
create the workflow using::

    pycbc_make_daily_workflow --local-config-files ${LOCAL_CONFIG_FILES} \
                              --start-time ${GPS_START_TIME} --ouput-dir ${LOGPATH} \
                           --config-overrides workflow:workflow-html-basedir:${HTMLDIR} \
                                              workflow:workflow-asset-dir:${ASSETDIR} \
                                              
If you are using default installed configuration files then you can create the
workflow using::

    pycbc_make_daily_workflow --installed-config-files ${INSTALLED_CONFIG_FILES} \
                              --start-time ${GPS_START_TIME} --output-dir ${LOGPATH} \
                           --config-overrides workflow:workflow-html-basedir:${HTMLDIR} \
                                              workflow:workflow-asset-dir:${ASSETDIR} \

.. _weeklyahopeplan:

-----------------------------------------
Planning and submitting the workflow
-----------------------------------------
CD into the directory where the dax was generated::

    cd ${MONTHDIR}/${DAYDIR}/

From the directory where the dax was created, run the planning script::

    pycbc_submit_dax --dax daily_ahope.dax
    
This will plan an submit your workflow to the cluster.

-----------------------------------------
Monitor and Debug the Workflow
-----------------------------------------

To monitor the above workflow, one would run::

    pegasus-status /usr1/${USER}/log/${USER}/pegasus/daily_ahope/run0011
    
To get debugging information in the case of failures::

    pegasus-analyzer /usr1/${USER}/log/${USER}/pegasus/daily_ahope/run0011


