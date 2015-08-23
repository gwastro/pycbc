########################################################
``pycbc_make_sngl_workflow``: A single-ifo detchar tool
########################################################

===============
Introduction
===============

The executable ``pycbc_make_sngl_workflow`` is a tool used to analyse data
from single detector(s), with no coincidence stage. The output from this single
detector matched-filter runs can be used to identify times where the detectors
are producing glitches that have a large single-detector detection statistic
for the CBC searches. If the detectoralists can identify what causes these most
egregious glitches it can increase the overall search sensitivity.

=========================================
How to run ``pycbc_make_sngl_workflow``
=========================================

Here we document the stages needed to run ``pycbc_make_sngl_workflow``.

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

    INSTALLED_CONFIG_FILES="sngl/example_sngl.ini"

Now go down to :ref:`pycbcdailygenerate`.

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
No, I need to make a configuration file - Editing the example files
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

An example configuration file set is found here::

    /src/dir/pycbc/workflow/ini_files/sngl/example_sngl.ini

This file contains all the details needed to run ``pycbc_make_sngl_workflow``

.. note::

    If you are unfamiliar with pycbc workflows, look through these files.
    example_daily_lalapps.ini will look familiar if you are used to daily ihope workflows.

The ``example_sngl.ini`` example is set up to run on S6 data and analysing only H1 and L1.

If you want to run in this default configuration please jump down the "Generate the workflow".

If you want to analyse a different set of ifos, you will have to edit some additional options::

    [workflow]
    h1-channel-name = H1:GDS-FAKE_STRAIN
    l1-channel-name = L1:OAF-CAL_DARM_DQ

    [workflow-ifos]
    ; This is the list of ifos to analyse
    h1 =
    l1 =

    [workflow-datafind]
    datafind-h1-frame-type = H1_ER_C00_AGG
    datafind-l1-frame-type = L1_R

    [workflow-segments]
    segments-h1-science-name = H1:DMT-SCIENCE:1
    segments-l1-science-name = L1:DMT-ANALYSIS_READY:1
    segments-database-url = https://dqsegdb5.phy.syr.edu
    segments-veto-definer-url = https://www.lsc-group.phys.uwm.edu/ligovirgo/cbc/public/segments/ER6/H1L1V1-ER6_CBC_OAF_CAL_DARM_DQ.xml

To run through this

* The ``[workflow-ifos]`` section supplies which ifos will be analysed if data is found and available.
* The ``x1-channel-name`` options are the h(t) channel name in the frames.
* The ``datafind-x1-frame-type`` options are the type of the frames for use when calling ``gw_data_find``.
* The ``segments-x1-science-name`` options are the flag used to store science times in the segment database.
* The ``segments-database-url`` option points to the segment database.
* The ``segments-veto-definer-url`` option points to the url where the veto-definer file can be found.

The remaining options affect how the jobs run, these should not be edited unless you know what you are doing. To find out more details about the possible options for any stage of the workflow, follow the links at :ref:`workflowhomepage`.

Now you have configuration files and can follow the same instructions as above. That is: 

Copy the configuration files into your run directory::

    cp /path/to/example_sngl.ini .

and set the names of these configuration files in your path. If you have more than one configuration file they must be space separated::

    LOCAL_CONFIG_FILES="example_sngl.ini"

.. _pycbcdailygenerate:

-----------------------
Generate the workflow
-----------------------

When you are ready, you can generate the workflow. First we need to choose a start time. Here is an example::

    export GPS_START_TIME=1102089616
    export GPS_END_TIME=1102100616

Choose a name for the workflow. For example use the GPS times::

    export WORKFLOW_NAME=${GPS_START_TIME}-${GPS_END_TIME}

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

    export HTMLDIR=/home/${USER}/public_html/pycbc_sngl/${WORKFLOW_NAME}

You also need to tell the html page generation code how to link the static pages. For example::

    export HTMLLINK=/~${USER}/pycbc_sngl/${WORKFLOW_NAME}

If you are using locally edited or custom configuration files then you can
create the workflow using::

   pycbc_make_sngl_workflow --name ${WORKFLOW_NAME} \
                            --local-config-files ${INSTALLED_CONFIG_FILES} \
                            --config-overrides workflow:start-time:${GPS_START_TIME} \
                                               workflow:end-time:${GPS_END_TIME} \
                                               workflow:workflow-html-basedir:${HTMLDIR} \
                                               workflow:workflow-html-link:${HTMLLINK} \
                                               html:analysis-subtitle:${GPS_START_TIME}-${GPS_END_TIME}

If you are using default installed configuration files then you can create the
workflow using::

   pycbc_make_sngl_workflow --name ${WORKFLOW_NAME} \
                            --installed-config-files ${INSTALLED_CONFIG_FILES} \
                            --config-overrides workflow:start-time:${GPS_START_TIME} \
                                               workflow:end-time:${GPS_END_TIME} \
                                               workflow:workflow-html-basedir:${HTMLDIR} \
                                               workflow:workflow-html-link:${HTMLLINK} \
                                               html:analysis-subtitle:${GPS_START_TIME}-${GPS_END_TIME}

.. _weeklyahopeplan:

-----------------------------------------
Planning and submitting the workflow
-----------------------------------------
CD into the directory where the dax was generated::

    cd ${MONTHDIR}/${DAYDIR}/

From the directory where the dax was created, run the planning script::

    pycbc_submit_dax --dax daily_ahope.dax 
    
This will plan an begin running the workflow.

-----------------------------------------
Monitor and Debug the Workflow
-----------------------------------------

To monitor the above workflow, one would run::

    pegasus-status /usr1/${USER}/log/${USER}/pegasus/daily_ahope/run0011
    
To get debugging information in the case of failures::

    pegasus-analyzer /usr1/${USER}/log/${USER}/pegasus/daily_ahope/run0011


