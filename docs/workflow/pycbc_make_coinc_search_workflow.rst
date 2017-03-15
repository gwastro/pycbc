####################################################################################
``pycbc_make_coinc_search_workflow``: A workflow to search for gravitational waves
####################################################################################

===============
Introduction
===============

The executable ``pycbc_make_coinc_search_workflow`` is a tool used to search
for gravitational waves using data from a network of gravitational-wave
detectors.  It performs all of the necessary steps in the workflow, including
data management, template bank construction (if required), matched filtering
and signal-based vetoes, multi-detector coincidence, data-quality cuts,
background estimation, and software injections to test the search. Its
ultimate task is to determine whether or not a compact binary coalescence is
present in the given data.  The output is a webpage containing the plots that
can be used to understand the results of the analysis

.. _configurationfiles:

==================
Configuration file
==================

The behavior of the workflow is controlled by a configuration file (also known as an ``ini`` file) that is made up of three types of sections: workflow, the pegasus profile and the executable options. The workflow sections control how different parts of the the workflow hang together. The pegasus profile sections are equivalent to lines you would have in a condor_submit file (e.g. requirements, storage size etc). Anything you would do in condor you would do here. The third section type maps the options to an executable.

::

  [pegasus_profile]
  condor|accounting_group=ligo.dev.o1.cbc.bns_spin.pycbcoffline

This is used for keeping an account of the pycbc usage

::

  [workflow]
  ; https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/workflow/initialization.html
  h1-channel-name = H1:GDS-FAKE_STRAIN
  l1-channel-name = L1:GDS-CALIB_STRAIN
  file-retention-level = all_triggers

You tend to put things in here which will be referred to later (but you can leave it empty). This is a nice way to keep options the same without the need to repeatedly define them. Here the L1/H1 data channel name are given. We also add an option to keep all the triggers/plots the analysis produces

::

  [workflow-ifos]
  l1 =
  h1 =
  
Set up which detectors you are going to run over. A blank space after an equals sign denotes True.

::

  [workflow-datafind]
  ; See https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/workflow/datafind.html
  datafind-method = AT_RUNTIME_SINGLE_FRAMES
  datafind-h1-frame-type = H1_ER_C00_AGG
  datafind-l1-frame-type = L1_ER_C01_L1
  datafind-check-segment-gaps = update_times
  datafind-check-frames-exist = raise_error
  datafind-check-segment-summary = no_test
  
This section defines which frames we are going to use and employs different levels of checks to see whether the data exists, there are gaps etc. 

- ``‘datafind-method’`` states how we are going to find the frames. The ‘AT_RUNTIME_SINGLE_FRAMES’ means the executable returns a list of single frame files. You can however provide a cache file, but then the options need to be changed. 
- ``‘datafind-h1-frame-type’`` refers to the frame type the H1 channel name will be found in for the time specified. Same for L1. 
- ``‘datafind-check-segment-gaps’`` option checks to see if there are gaps in the segments from the segment database and the option ‘update_times’ will change the analysis times to skip over these gaps. 
- ``‘datafind-check-frames-exist’`` checks to see if the frames you are looking at actually exists, and if they don’t the ‘raise_error’ option will stop the workflow. 
- ``‘datafind-check-segment-summary’`` Checks the segment summary table and makes sure that the frames exist for all times that the segments are known

::

  [workflow-segments]
  ; See https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/workflow/segments.html
  segments-method = AT_RUNTIME
  segments-h1-science-name = H1:DMT-SCIENCE:1
  segments-l1-science-name = L1:DMT-ANALYSIS_READY:1
  segments-database-url = https://dqsegdb5.phy.syr.edu
  segments-veto-definer-url = https://www.lsc-group.phys.uwm.edu/ligovirgo/cbc/public/segments/ER6/H1L1V1-ER6_GDS_CALIB_STRAIN.xml
  segments-science-veto = 1
  segments-veto-groups = 
  segments-final-veto-group = 1

This section does a series of checks to the segment database for the segments you need for your analysis. 

- ``‘segments-method’`` option should not change. 
- ``‘segments-h1-science-name’`` option specifies the segment name at LHO we consider to flag science time. The same is given for L1. 
- ``‘segments-data-url’`` specifies the url for the segment database we want to query. 
- ``‘segments-veto-definer-url’`` is the url for the veto definer file we want to use for the search. 
- ``‘segments-science-veto’`` species which category of veto you want to eliminate from your search before it is performed to consider the data science. In this instance, 1 denotes that all the times of Cat 1 vetoes. Time vetoed here is not used in any part of the analysis, and is treated as if it were not collected. 
- ``‘segments-veto-groups’`` is an option you can populate with different veto categories and diagnostic plots will be made after each veto is employed. 
- ``‘segments-final-veto-group’`` is an important option as the vetoes defined here will be used to remove triggers from the search before coincidence is performed. An option of 1 will remove all Cat 1 veto times from the analysis before it is performed. If you want to add cat 2 then the option is 12.

::

  [workflow-tmpltbank]
  ; See https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/ahope/template_bank.html
  tmpltbank-method=PREGENERATED_BANK
  tmpltbank-pregenerated-bank=/home/jveitch/projects/mdc/spin/tmpltbanks/nonspin/BNS_NonSpin_30Hz_earlyaLIGO.xml

This section specifies which template bank to use

- ``’tmpltbank-method’`` option specifies whether you want to use a regenerated bank or to make it on the fly. In O1 we will be us a pregenerated bank. 
- ``‘tmpltbank-pregnerated-bank’`` specifies the location of the xml with the pregenerated bank. Note that this exact location is only valid for SUGAR, and that in general one must provide their own template bank. 

::

  [workflow-splittable]
  splittable-method = IN_WORKFLOW
  splittable-num-banks = 2

This section sets the options for splitting the bank to help with computational costs.

- ``‘splittable-method’`` tells you the method by which to split the bank, in this instance it is IN_WORKFLOW. If you do not want to split the bank, change this option to NOOP
- ``‘splittable-num-banks’`` specifies how many banks to split the original bank into.

::

  [workflow-matchedfilter]
  ; See https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/ahope/matched_filter.html
  matchedfilter-method=WORKFLOW_INDEPENDENT_IFOS
  min-analysis-segments = 5
  max-analysis-segments = 5
  output-type = hdf

This section defines how the matched filter is going to be performed. Whether it is going to be independent for each detector, and also how the analysis is actually going to be separated in to chunks given the data available.

- ``‘matched-filter-method’`` defines where the data is going to be separated and searched over, in this instance the data for each IFO will be considered independently and in the workflow
- ``‘min-analysis-segments’`` defines the minimum number of overlapping chunks you separate the data in to to analyze. This is a proxy for segment length. In this instance 5 has been stated. Therefore if the data cannot be split in to 5 overlapping chunks the code skips over the data. To understand how much time this is you need to look in the [inspiral] options and consider the segment-length and padding options specified. ‘max-analysis-segments’ is the same but for the maximum number of overlapping chunks. Be aware if you lower/raise either of these numbers you will affect the psd estimation. 
- ``‘output-type’`` is the format of the output trigger files from the matched filter search

::

  [workflow-coincidence]
  ; See https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/ahope/coincidence.html
  parallelization-factor = 10

This part of the workflow looks for coincidence between templates between detectors. All coincidences are kept. If you have a large template bank you probably want make the ``‘parallelization-factor’`` large

::

  [workflow-injections]
  injections-method=IN_WORKFLOW

This section deals with software injections. Here you are specifying whether to use either pregenerated injections sets or ones made within the workflow itself. In this case, we will use one that is created within the workflow. 

::

  [executables]
  ; setup of condor universe and location of executables
  inspiral          = ${which:pycbc_inspiral}
  injections = ${which:lalapps_inspinj}
  splittable = ${which:pycbc_splitbank}
  segment_query = ${which:ligolw_segment_query_dqsegdb}
  segments_from_cats = ${which:ligolw_segments_from_cats_dqsegdb}
  llwadd = ${which:ligolw_add}
  ligolw_combine_segments = ${which:ligolw_combine_segments}
  bank2hdf = ${which:pycbc_coinc_bank2hdf}
  hdfinjfind = ${which:pycbc_coinc_hdfinjfind}
  coinc = ${which:pycbc_coinc_findtrigs}
  statmap = ${which:pycbc_coinc_statmap}
  statmap_inj = ${which:pycbc_coinc_statmap_inj}
  plot_sensitivity = ${which:pycbc_page_sensitivity}
  plot_foundmissed = ${which:pycbc_page_foundmissed}
  plot_snrifar = ${which:pycbc_page_snrifar}
  page_foreground = ${which:pycbc_page_foreground}
  page_injections = ${which:pycbc_page_injtable}
  hdf_trigger_merge = ${which:pycbc_coinc_mergetrigs}
  plot_snrchi = ${which:pycbc_page_snrchi}
  plot_coinc_snrchi = ${which:pycbc_page_coinc_snrchi}
  plot_segments = ${which:pycbc_page_segments}
  results_page = ${which:pycbc_make_html_page}

This section defines where each of the executables live; it tells the workflow which files to process. It might be worth checking you can find all of these paths before you set the code running. 

The following options are those associated to a given executable. 

::

  [llwadd]
  [datafind]
  urltype=file

This is the format for the return of the data find executable - you want a file.

::

  [segments_from_cats]

Some sections are left empty. That is fine, but you have to define each option otherwise the code will complain

::

  [ligolw_combine_segments]

  [splittable]
  ; options for splittable job
  random-sort =

This option randomly sorts the bank to be split up before processing

::

  [injections]
  waveform = SpinTaylorT4threePointFivePN
  
Define the waveforms you want to use for injections

::

  [injections-bnslininj]
  f-lower = 20
  min-distance = 1000
  max-distance = 150000
  d-distr = uniform
  l-distr = random
  i-distr = uniform
  min-mass1 = 1.0
  max-mass1 = 3.1
  min-mass2 = 1.0
  max-mass2 = 3.1
  m-distr = componentMass
  min-mtotal = 2.0
  max-mtotal = 6.2
  disable-spin =
  time-step = 89.155
  time-interval = 10
  seed = 1234
  
These are the injections parameters you want to define. Only defining ones which aren’t so obvious

- ``f-lower`` = low frequency cut off
- ``min-distance`` =  (kpc)
- ``max-distance`` = (kpc)
- ``d-distr`` = the distance distribution of the injections
- ``l-distr`` = the distribution of injections in the sky
- ``i-distr`` = inclination of the injection 
- ``time-step`` = time between injections. This can be whatever time you want, but remember if the injections are too close together you can screw up your psd estimation. ~90s seems ok. 
- ``time-interval`` = time interval to inject the signal. It will not always be exactly at time-step, but at a time of time-step +/- random_number(0,time-interval)
- ``seed`` = random seed, choose different numbers to get different realizations of the same background distribution

::

  [inspiral]
  ; inspiral analysis parameters -- added to all inspiral jobs
  chisq-bins = 256
  snr-threshold = 5.0
  approximant = SPAtmplt
  order = 7
  cluster-method = window
  cluster-window = 1.0
  segment-length = 512
  segment-start-pad = 64
  segment-end-pad = 16
  psd-estimation = median
  psd-segment-length = 16
  psd-segment-stride = 8
  psd-inverse-length = 16
  strain-high-pass = 30
  pad-data = 8
  processing-scheme = mkl
  sample-rate = 4096
  filter-inj-only =
  low-frequency-cutoff = 40
  
These are the parameters you want to define for the inspiral search

- ``chisq-bins`` = number of chisq bins for the standard Bruce Allen chisq
- ``snr-threshold`` = SNR threshold
- ``approximant`` = approximation you want to use. SPAtmplt is stationary phase approximation template which is a fast implementation of Taylor F2.
- ``order`` = PN order, the numbers are double the order. So 7=3.5PN
- ``cluster-method`` = method over which to identify the loudest trigger - in this case a window
- ``cluster-window`` = take a 1 second window around the loudest trigger
- ``segment-length`` = the length of a segment you want to analyze. Remember previously we mention we want 5 overlapping segments
- ``segment-start-pad`` = the amount of time we want to pad the start of the data by. In this instance we want to not use the first 64 seconds of data, as it will contain errors from filtering. This takes in to account the length of time we lose due to PSD corruption (16s) and the wrap around effect we have due to the template (48s) 
- ``segment-end-pad`` = the amount of time we want to pad the end of the data by. See above.
- ``psd-estimation`` = the method by which we want to estimate the psd
- ``psd-segment-length`` = length of time used in each psd calculation
- ``psd-segment-stride`` = time spacing between each psd calculation. 16s length with 8s stride implies a 50% overlap
- ``psd-inverse-length`` = time length used to truncate the inverse FFT (that is, the time domain realization) of the psd 
- ``strain-high-pass`` = high pass filter applied to strain data before psd estimation
- ``pad-data`` = 8 second padding added to beginning of data to account for filter corruption for resampling and high-pass before data is broken up into chunks
- ``processing-scheme`` = indicates which software to use for processing (MKL = math kernel library made by Intel)
- ``sample-rate`` = sample rate of data (will be down sampled in workflow)
- ``filter-inj-only`` = Use only segments with injections in them for matched filter
- ``low-frequency-cutoff`` = low frequency limit for the matched filter search

::

  [inspiral-h1]
  ; h1 specific inspiral parameters
  channel-name = ${workflow|h1-channel-name}

Specify the name of the channel you want to run the inspiral analysis over for H1. Here we are referring back to the name in the workflow module

::

  [inspiral-l1]
  ; l1 specific inspiral parameters
  channel-name = ${workflow|l1-channel-name}

  [bank2hdf]
  [trig2hdf]

  [coinc]
  coinc-threshold = 0.000

Here we are doing exact match coincidence. So we take the light travel time between detectors and look for triggers which are coincident within this time window. The threshold defines if you want to extend the window.

::

  [coinc-full]
  decimation-factor = 1000
  loudest-keep = 200
  timeslide-interval=1.1

This section concerns time slides without injections, and its purpose is to keep a small number of timesmlide triggers for background estimation. Time slides are done at all relative offsets that are multiple of the 'timeslide-interval', which is defined here to be 1.1 seconds. We don’t store all the coincident triggers due from time slides. We keep 200 of the loudest triggers from each template time slide, given by the second option, which gives a good estimation of the background at low FAR. The top option specifies for which timeslides we will keep all triggers, to get an overall estimation of background (not just the loudest). In this instance we would keep the triggers from 1000th, 2000th, 3000th timeslide. 

::

  [coinc-injfull&coinc-fullinj]
  timeslide-interval={coinc-full:timeslide-interval}
  loudest-keep-value = 8.5
  cluster-window = {statmap|cluster-window}

This section concerns time slides with injections in the data. We assume only one injection will be coincident with a timeslide (done every 1.1 seconds - see first option) trigger and we keep its coincidence if its ranking statistic (newSNR) > 8.5 as specified in the second option. This is to limit storage of unimpactful triggers only. 

::

  [coinc-injinj]

  [pegasus_profile-statmap&pegasus_profile-statmap_inj]
  condor|request_memory = 20GB

This is the amount of memory the jobs might take

::

  [statmap&statmap_inj]
  veto-window = 0.050
  cluster-window = 10.0

This controls the final clustering after all coincidence testing. The ``cluster-window`` indicates the time window used for clustering.
The ``veto-window`` is used to remove all coincident zero-lag triggers so that they aren't included in background estimation

::

  [hdfinjfind]
  injection-window = 1.0

The rest of the config file concerns plotting formats

::

  [page_foreground]
  [plot_snrifar]

  [plot_snrchi]
  [plot_coinc_snrchi]
  [plot_coinc_snrchi-inj]
  [plot_coinc_snrchi-bkg]
  background-front=
  [plot_coinc_snrchi-inj&plot_coinc_snrchi-bkg&plot_snrchi]
  newsnr-contours =  6 8 10

  [plot_sensitivity]
  sig-type = ifar
  sig-bins = 1 3 10 30 100 300 1000 3000 10000 30000 100000

  [plot_sensitivity-mchirp]
  bin-type =  mchirp 
  bins = 0.89 1.31 1.74 2.17 2.60 
  min-dist = 40 
  max-dist = 120 
  dist-bins = 50 

  [plot_sensitivity-mtotal]
  bin-type =  total_mass
  bins = 2 2.4 3.2 4 6 
  min-dist = 40 
  max-dist = 120 
  dist-bins = 50 

  [plot_sensitivity-spin]
  bin-type =  spin
  bins = -0.4 -0.2 0.2 0.4 
  min-dist = 40 
  max-dist = 120 
  dist-bins = 50 

  [plot_sensitivity-mchirp_binless]
  bin-type =  mchirp 
  bins = 0.89 1.31 1.74 2.17 2.60 
  min-dist = 40 
  max-dist = 120 

  [plot_sensitivity-mtotal_binless]
  bin-type =  total_mass
  bins = 2 2.4 3.2 4 6 
  min-dist = 40 
  max-dist = 120 

  [plot_sensitivity-spin_binless]
  bin-type =  spin
  bins = -0.4 -0.2 0.2 0.4 
  min-dist = 40 
  max-dist = 120  

  [plot_foundmissed]
  [plot_foundmissed-mchirp]
  axis-type=mchirp
  dynamic=
  [plot_foundmissed-chirpdistmchirp]
  axis-type=mchirp
  dynamic=
  distance-type=chirp_distance
  [plot_foundmissed-time]
  axis-type=time
  dynamic=

  [plot_foundmissed-mchirp_static]
  axis-type=mchirp
  log-distance=
  [plot_foundmissed-chirpdistmchirp_static]
  axis-type=mchirp
  distance-type=chirp_distance
  log-distance=
  [plot_foundmissed-time_static]
  axis-type=time
  log-distance=

  [hdf_trigger_merge]
  [pegasus_profile-hdf_trigger_merge]
  condor|request_memory = 10GB

  [page_injections]
  [plot_segments]

  [results_page]
  analysis-title="PyCBC Coincident Analysis"
  analysis-subtitle="..."
  

.. _coincworkflowgenerate:

=======================
Generating the workflow
=======================

The workflow is generated by running the script ``pycbc_make_coinc_search_workflow``. This program takes the command line arguments

.. command-output:: pycbc_make_coinc_search_workflow --help

The configuration files can either be passes as local files, or given as URLs
to specific configuration files managed for an analysis. For example, to
generate a workflow to search two weeks of S6D data and place the results in
your ``public_html`` directory, run the command::

    pycbc_make_coinc_search_workflow --workflow-name s6d_chunk3 --output-dir output \
      --config-files https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-config/download/master/S6/pipeline/s6_run_pycbc_er8_pre_release.ini \
      https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-config/download/master/S6/pipeline/executables.ini \
      https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-config/download/master/S6/pipeline/injections.ini \
      https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-config/download/master/S6/pipeline/data_S6.ini \
      https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-config/download/master/S6/pipeline/gps_times_s6d_big_dog_two_weeks.ini \
      --config-overrides "results_page:output-path:${HOME}/public_html/s6/s6d-big-dog-weeks"

The configuration ``results_page:output-path`` can be changed appropriately to
set the output web page location.

.. note::

   To use released executables for production analysis, you should specify
   the URL to an ``executables.ini`` file from the 
   `PyCBC Software repository <https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-software>`_.

.. _coincworkflowplan:

====================================
Planning and Submitting the Workflow
====================================

Pegasus is used to plan and submit the workflow. To involve Pegasus to plan a
PyCBC workflow, you use the command ``pycbc_submit_dax`` which takes the
command line arguments

.. command-output:: pycbc_submit_dax --help

Note that  you are running on a resource that mandates accounting, then you
will also need to add a valid tag with the ``--accounting-tag`` command line
argument. Please see
`the LDG accounting page <https://ldas-gridmon.ligo.caltech.edu/ldg_accounting/user>`_. to
determine the correct tags. These can be applied by adding the following line
to your submit invocation.

For example, to plan and submit the workflow in the example above, change to the directory that you specified with the ``--output``
command line option to ``pycbc_make_coinc_search_workflow`` and plan and submit
the workflow::

    cd output
    pycbc_submit_dax --accounting-group ligo.dev.o1.cbc.explore.test --dax s6d_chunk3.dax

.. note::

    The above example uses the accounting tag ``ligo.dev.o1.cbc.explore.test``
    which should not be used in practice.

You can monitor the status of the workflow with Pegasus Dashboard, or the
other Pegasus tools described below. 

If the workflow runs successfully, the output will be place under the
directory specified by ``results_page:output-path`` when the workflow is
complete.

-------------------------------------------------------------------------------------------------------------------------------------------
Monitor and Debug the Workflow (`Detailed Pegasus Documentation <https://pegasus.isi.edu/wms/docs/latest/tutorial.php#idm78622034400>`_)
-------------------------------------------------------------------------------------------------------------------------------------------

To monitor the above workflow, one would run::

    pegasus-status /usr1/ahnitz/log/ahnitz/pegasus/weekly_ahope/run0011
    
To get debugging information in the case of failures.::

    pegasus-analyzer /usr1/ahnitz/log/ahnitz/pegasus/weekly_ahope/run0011

-----------------------------
Pegasus Dashboard
-----------------------------

The `pegasus dashboard <http://pegasus.isi.edu/wms/docs/latest/ch02s11.php>`_ is a visual and interactive way to get information about the progress, status, etc of your workflows.

The software can be obtained from a separate pegasus package here <https://github.com/pegasus-isi/pegasus-service>.

Pegasus Dashboard is currently installed on sugar. To view your Pegasus Dashboard, in a browser go to::

    https://sugar.phy.syr.edu/pegasus/u/albert.einstein

This shows a page that has a table of all your workflows that were submitted from sugar. You can view the details of a workflow by clicking on the link in the Workflow Details column of the table.

Clicking on the Workflow Details link will take you to a webpage that gives a high-level overview of the workflow, telling you how many many jobs succeeded, fail, the submit directory, etc. There is a table with tabs at the bottom of the page. If you click the tabs Failed, Running, and Successful the page will generate a table that lists all the failed, running, and successful jobs in the workflow respectively. You also have the ability to search the table for a particular kind of job using the Search bar.

You can view the details of a job by clicking the link in the Job Name column. This will take you to a Job Details page. This page will tell you where to find stdout files, stderr files, how much wall clock time the job took to run, etc. There is a table at the bottom of the page with a Failed and Successful tab. If you click on the respective tab, it will list all invocations of that job. You can click on the link in the Invocations column for more information.

On the Invocation Details page there is information about the command line arguments, executable path, CPU time, wall clock time, etc.

In certain cases, the pegasus monitor daemon may crash and this could result in
invalid or nonsensical information on the dashboard (e.g. a cumulative
computing time of None). This problem can be solved by running
``pegasus-plots`` on the workflow directory: the command should tell you what
to do. Typically this will be running ``pegasus-monitord`` in replay mode (see
its man page).

-----------------------------
Pegasus Analyzer
-----------------------------

The `pegasus analyzer <http://pegasus.isi.edu/wms/docs/trunk/cli-pegasus-analyzer.php>`_ is a command-line tool for reporting sucessful and failed jobs.

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

.. _weeklyahopereuse:

======================================
Reuse of data from a previous workflow
======================================

One of the features of Pegasus is reuse the data products of prior runs.
This can be used to e.g. expand an analysis or recover a run with mistaken settings without
duplicating work. The steps below explain how to do this.

------------------------------------
Setting up a workflow for data reuse
------------------------------------

The first step is to generate a new workflow that performs the analysis that
you would like to do. This workflow should be generated in a new directory so that it does not overwrite data from your previous workflows.
Data reuse happens at the ``pycbc_submit_dax`` step, so
first run ``pycbc_make_coinc_search_workflow`` to build a new workflow,
following the instructions in the section :ref:`coincworkflowgenerate` of this
page.

**Stop** before you plan and submit the workflow with ``pycbc_submit_dax``.
You will pass an additional file to ``pycbc_submit_dax`` using the
``--cache-file`` option with a list of files that Pegasus can re-use from a
previous run.  The Pegasus Workflow Planner will reduce the workflow
using this cache file. Reduction works by deleting jobs from the workflow
whose output files have been found in some location in this cache file.

The key to data reuse is building the cache file passed to ``pycbc_submit_dax``. This file maps a file created in the workflow to a URL and a site where that URL can be found. The syntax of the cache file is plain ASCII with each line in the file giving the location of a file in the format::

    LOGICAL_FILE_NAME PHYSICAL_FILE_URL pool="SITE"

where ``LOGICAL_FILE_NAME`` is the name of the file as it appears in the
workflow. This should include any subdirectory path used by the workflow to organize files in the case of, e.g.,
``INSPIRAL`` files but it should not be the absolute path to the file. ``PHYSICAL_FILE_URL`` is a
full URL where the file can be found, and ``SITE`` is the site on which that URL
resides.

The URI in the ``PHYSICAL_FILE_URL`` can be any of the URIs that Pegasus
recognizes. The URIs ``file://``, ``gsiftp://``, and ``http://`` are likely
the most useful. Pegasus will take care of adding transfer jobs for
``gsiftp://`` and ``http://`` URIs, if the data is not available locally.

The string ``SITE`` is a hint that tells Pegasus on which site the
``PHYSICAL_FILE_URL`` can be found. The ``SITE`` string should be one of the
names used by ``pycbc_submit_dax`` to identify the cluster where jobs are run.
In practice there are only two execution sites used by PyCBC workflows:

1. ``local`` which is the regular Condor pool on the local cluster where the workflow is being run from. This is typically used when re-using data that exists on the filesystem of the local cluster.
2. ``osg`` which is the Open Science Grid pool, as described in :ref:`weeklyahopeosg` below. This is only used if the data to be re-used is accessible via the ``/cvmfs`` filesystem.

If the ``SITE`` string for a file matches the site where a job will be run,
then Pegasus assumes that the file can be accessed locally via the regular
file open commands. If the ``SITE`` string does not match the site where a job
will be run, then Pegasus adds transfer jobs to the workflow to move the file
to the site where it will be needed by a job.

To tell Pegasus that the file is neither accessible via file open on the
``local`` submit host nor on the ``osg`` pool, then the ``SITE`` string can be
set to ``remote``. This tells Pegasus that the file is neither on the
``local`` or the ``osg`` site and so Pegasus must add file transfer jobs to
fetch the file from some other site.  This ``SITE`` attribute is needed
beacuse a map between the job execution site and the location of the file
might not be obvious from the hostname in the ``PHYSICAL_FILE_URL``.  

The following rule should be helpful when chosing the ``SITE`` string:

* If you are re-using a file that is available locally with a ``file://`` URI in its ``PHYSICAL_FILE_URL`` (or has an implicit ``file://`` URI since the ``PHYSICAL_FILE_URL`` starts with a ``/``) then the string ``SITE`` should be set to ``local``.
* If you are re-using a file from another cluster, e.g. you are on the Syracuse cluster and want to re-use data from AEI Atlas cluster, then the string ``SITE`` should be set to ``remote`` for that file. In this case, the URI in ``PHYSICAL_FILE_URL`` will be either ``gsiftp://`` or ``http://`` depending on how the file can be accessed.

To illustrate this, an example of a simple cache file containing four files for re-use from the ``local`` site is::

    H1-VETOTIME_CAT3-1169107218-1066800.xml file://localhost/home/dbrown/projects/aligo/o2/analysis-4/o2-analysis-4/output/results/1._analysis_time/1.01_segment_data/H1-VETOTIME_CAT3-1169107218-1066800.xml pool="local"
    L1-VETOTIME_CAT3-1169107218-1066800.xml file://localhost/home/dbrown/projects/aligo/o2/analysis-4/o2-analysis-4/output/results/1._analysis_time/1.01_segment_data/L1-VETOTIME_CAT3-1169107218-1066800.xml pool="local"
    116912/H1-INSPIRAL_FULL_DATA_JOB0-1169120586-1662.hdf file://localhost/home/dbrown/projects/aligo/o2/analysis-4/o2-analysis-4/output/full_data/H1-INSPIRAL_FULL_DATA_JOB0-1169120586-1662.hdf pool="local"
    116912/H1-INSPIRAL_FULL_DATA_JOB1-1169120586-1662.hdf file://localhost/home/dbrown/projects/aligo/o2/analysis-4/o2-analysis-4/output/full_data/H1-INSPIRAL_FULL_DATA_JOB1-1169120586-1662.hdf pool="local"

Note that the ``LOGICAL_FILE_NAME`` for the veto files is just the name of the
file, but for the two inspiral files it contains the subdirectory that the
workflow uses to organize the files by GPS time. In the case of this file Pegasus will delete from the workflow the jobs that create the files ``H1-VETOTIME_CAT3-1169107218-1066800.xml``, ``L1-VETOTIME_CAT3-1169107218-1066800.xml``, ``116912/H1-INSPIRAL_FULL_DATA_JOB0-1169120586-1662.hdf``, and ``116912/H1-INSPIRAL_FULL_DATA_JOB1-1169120586-1662.hdf`` when it plans the workflow. Insted, the data will be re-used from the URLs specified in the cache. Since ``site="local"`` for these files, Pegasus expects that the files all exist on the host where the workflow is run from.

To re-use data from a remote cluster, the URLs must contain a file transfer
mechanism and the ``SITE`` should be set to ``remote``. For example, if the
files listed in the example above are available on
``sugwg-condor.phy.syr.edu`` and you want to re-use them in a workflow on the
AEI Atlas cluster, then the cache file would contain::

    H1-VETOTIME_CAT3-1169107218-1066800.xml gsiftp://sugwg-condor.phy.syr.edu/home/dbrown/projects/aligo/o2/analysis-4/o2-analysis-4/output/results/1._analysis_time/1.01_segment_data/H1-VETOTIME_CAT3-1169107218-1066800.xml pool="remote"
    L1-VETOTIME_CAT3-1169107218-1066800.xml gsiftp://sugwg-condor.phy.syr.edu/home/dbrown/projects/aligo/o2/analysis-4/o2-analysis-4/output/results/1._analysis_time/1.01_segment_data/L1-VETOTIME_CAT3-1169107218-1066800.xml pool="remote"
    116912/H1-INSPIRAL_FULL_DATA_JOB0-1169120586-1662.hdf gsiftp://sugwg-condor.phy.syr.edu/home/dbrown/projects/aligo/o2/analysis-4/o2-analysis-4/output/full_data/H1-INSPIRAL_FULL_DATA_JOB0-1169120586-1662.hdf pool="remote"
    116912/H1-INSPIRAL_FULL_DATA_JOB1-1169120586-1662.hdf gsiftp://sugwg-condor.phy.syr.edu/home/dbrown/projects/aligo/o2/analysis-4/o2-analysis-4/output/full_data/H1-INSPIRAL_FULL_DATA_JOB1-1169120586-1662.hdf pool="remote"

Note that the URL now contains ``gsiftp://sugwg-condor.phy.syr.edu`` rather
than ``file://localhost`` and the files are listes as ``pool="remote"`` rather
than ``pool="local"``. Pegasus will re-use these data files adding in
file transfer jobs to the workflow to move them into the appropriate
locations.

Once a cache file has been constructed, to enable data re-use, you follow the
standard instructions for planning and submitting the workflow in the section
:ref:`coincworkflowplan`, but add the ``--cache-file`` argument that points to
the cache file that you have created. For example:: 

    pycbc_submit_dax --cache-file /path/to/prior_data.map --accounting-group ligo.dev.o1.cbc.explore.test --dax s6d_chunk3.dax

will use the URLs from the file ``/path/to/prior_data.map`` to implement
data re-use and subsequent workflow reduction. If more than once cache file is
provided, pass the paths as a comma separated list to ``pycbc_submit_dax``::

    pycbc_submit_dax --cache-file /path/to/prior_data.map,/path/to/other.map --accounting-group ligo.dev.o1.cbc.explore.test --dax s6d_chunk3.dax

Which file URLs should be included in the reuse cache? There is no single
correct way of deciding this, as it depends on exactly what you are trying to do. The sections
below explain how to do this for a few common situations.

.. note::

    The ``[workflow]`` section of the ini configuration file contains an
    option ``file-retention-level``. This is commonly set to ``all_files`` or
    ``all_triggers``, in which case the data products re-used will be copied
    from the input locations and stored into the output location of the new
    workflow when the new workflow is run with data re-use. This can be
    wasteful of disk space, so you may want to set this option to either
    ``merged_triggers`` or ``results`` to store a smaller sub-set of the
    workflow's data products. These setting will allow the use of data from
    a previous run, but not make duplicate copies of intermediate data files.
    See the documentation under :ref:`workflowconfigparsermod` for more
    details of the ``file-retention-level`` configuration option.

.. note::

    At present you *cannot* re-use ``.dax`` and ``.map`` files from a previous
    run. A workflow using data reuse must regenerate and re-run any sub-daxes
    from scratch. If you re-use a ``.map`` file rather than re-generating it,
    then the new workflow will write results files in the location of the old
    workflow. All of the examples below use an ``egrep -v '(dax|map)'`` to
    filter out these files.

.. _workflow_rerun_extend:

-------------------------------------------------
Extending the GPS end time of a previous workflow
-------------------------------------------------

A common mode of data re-use is to extend the GPS end time of a previous
workflow to generate a new result page that e.g. extends the analysis by a few
days. This assumes that: 

* The previous workflow completed successfully.

* There are no changes to the workflow configuration file, other than incrementing the end time of the workflow.

In this case, first re-run ``pycbc_make_coinc_search_workflow`` to build the
new workflow. The normal file retention level will copy a lot of reused data
from the previous workflow directory into the new workflow directory. If you
do not want to do this, use a ``--config-override`` to change the value of
``workflow:file-retention-level`` as described on the page
:ref:`workflowconfigparsermod`.

Then create a cache file in the following way:

1. Locate the PyCBC result page for the workflow that you wish to extend.

2. In the menu under **Section 8: Workflow**, locate the **Output map** section (usually Section 8.06) and open that page.

3. This page will show three output cache files that contain the URLs of the data created by the workflow. Locate the file that ends ``main.map`` and download it by clicking on the **Link to file**. This file contains the main intermediate and output data products of the workflow.

4. Edit this file so that it only contains the output of the ``pycbc_inspiral`` jobs, i.e. delete all of the lines that do not match the pattern ``*INSPIRAL*hdf``. You can do this in a text editor, or with your favorite combination of UNIX ``grep``, ``sed``, ``awk``, or ``perl`` commands.
For example::

    egrep 'INSPIRAL.*hdf' /path/to/downloaded/workflow-main.map > inspiral_files.map

will pull out all cache file lines for the outputs of ``pycbc_inspiral`` files and write them to a new cache file called ``inspiral_files.map``.  

5. If the files in the new cache file exist locally on the cluster where you are submitting the workflow, then the cache file is complete. If they do not, you will need to modify the file to change the ``PHYSICAL_FILE_URL`` to a valid ``gsiftp://`` or ``http://`` URL on the remote cluster, and change ``pool="local"`` to ``pool="remote"``. Again, these changes can be made with a text editor or UNIX shell tools. For example, if the file URLs begin with ``/home/dbrown`` and they are on the Syracuse cluster, to run on Atlas you would use the following ``sed`` commands to change the ``SITE`` and the URI in the cache file::

    sed 's/pool="local"/pool="remote"/g' inspiral_files.map > inspiral_files.map.tmp
    sed 's+/home/dbrown+gsiftp://sugwg-condor.phy.syr.edu/home/dbrown+g' inspiral_files.map.tmp > inspiral_files.map
    rm inspiral_files.map.tmp

6. Finally, copy the file ``inspiral_files.map`` to your new workflow directory and then run ``pycbc_submit_dax`` as usual, giving the path to ``inspiral_files.map`` as the ``--cache-file`` argument.

---------------------------------------------------
Re-running a workflow using a new veto definer file
---------------------------------------------------

Data reuse can be used to re-running a workflow with a new veto definer file, assuming that:

* The previous workflow completed successfully.
* No changes to the configuration file are made, other than changing the ``segments-veto-definer-url`` in the ``[workflow-segments]`` section of the workflow configration file (although the GPS end time can also be extended at the same time, if necessary).

In this case, first re-run ``pycbc_make_coinc_search_workflow`` to build the
new workflow. The normal file retention level will copy a lot of reused data
from the previous workflow directory into the new workflow directory. If you
do not want to do this, use a ``--config-override`` to change the value of
``workflow:file-retention-level`` as described on the page
:ref:`workflowconfigparsermod`.

Then create the cache file as follows:

1. Locate the PyCBC result page for the workflow that you wish to extend.

2. In the menu under **Section 8: Workflow**, locate the **Output map** section (usually Section 8.06) and open that page.

3. This page will show three output cache files that contain the URLs of the data created by the workflow. Locate the file that ends ``main.map`` and download it by clicking on the **Link to file**. This file contains the main intermediate and output data products of the workflow.

4. If only category 2 and higher vetoes have change, remove the output files that match the following strings from the output map file: 

  * ``VETOTIME`` to remove the files containing the old veto segments.
  * ``LIGOLW_COMBINE_SEGMENTS`` to remove the files that combine the veto segments into categories.
  * ``CUMULATIVE_CAT_12H_VETO_SEGMENTS`` to remove the files that contain times to veto.
  * ``COINC`` to remove the output of the coincidence code.
  * ``FIT`` to remove the background bin statistic results.
  * ``STATMAP`` to remove the detection statistic ranking output.
  * ``INJFIND`` to remove the results of software-injection tests.
  * ``PAGE`` to remove the results make with the loudest events.
  * ``FOREGROUND_CENSOR`` to remove the veto files used to remove events from the closed box plots.
  * ``html`` to remove any output web pages genereated.
  * ``png`` to remove any output plots generated.
  * ``dax`` to remove any follow-up workflows generated.

This can be acomplished with the following command::

    egrep -v '(VETOTIME|LIGOLW_COMBINE_SEGMENTS|CUMULATIVE_CAT_12H_VETO_SEGMENTS|COINC|FIT|STATMAP|INJFIND|PAGE|FOREGROUND_CENSOR|html|png|dax)' /path/to/main.map > /path/to/reuse_cache.map

If category 1 vetoes have changed, you must also remove files matching ``PSD``, ``OPTIMAL``, and ``MERGE`` to remove the PSD estimation jobs, the jobs that compute the optimal SNR of injections, and the merged single-detector inspiral trigger files which may also change if the category 1 vetoes change.

6. Copy the file ``reuse_cache.map`` to your new workflow directory and then run ``pycbc_submit_dax`` as usual, giving the path to ``reuse_cache.map`` as the ``--cache-file`` argument.

----------------------------
Re-running a failed workflow
----------------------------

Occasionally it may be necessary to use data from a partially completed
workflow, e.g. if there a bug in an executable and you wish to re-run the
workflow with a new version of the executable. If the workflow failed, no
results web page will have been generated and the output data may not have
been copied to the locations in ``main.map``. To re-use data from a previous
failed workflow, you need to create a cache file containing the completed jobs
from the previous workflow. 

To do this, ``cd`` into the ``local-site-scratch/work`` directory of your
failed workflow. For example, if you used ``--output-dir output`` when
planning the workflow, and then run the command::

    cd /path/to/workflow/output/local-site-scratch/work

Once in this directory there should be a directory that ends with
``main_ID0000001`` (e.g. ``my-workflow-main_ID0000001``) Change into that
directory.

Once in the ``main_ID0000001`` directory, run the command::

    for pfn in `find . -type f | sed 's+^./++g'` ; do echo $pfn file://`pwd`/$pfn pool=\"local\" ; done | egrep -v '(dax|map)' > /path/to/partial_workflow.map

changing ``/path/to`` to a location where you want to save the cache.
 
Now you can than use the ``partial_workflow.map`` cache file as the ``--cache-file`` argument to ``pycbc_submit_dax``.

-----------------------------------------------
Using partial products from a previous workflow
-----------------------------------------------

If you are changing the configuration parameters of a workflow, then you can
build a cache file from a previous ``main.map`` file or the files under
``local-site-scratch``, but you will need to filter the cache file to remove
the files for jobs that have a changed configuration.  Here are a few
examples:

* If you are changing the configuration of ``pycbc_inspiral`` you must regenerate almost all the files in the workflow so it easier to start from scratch.

* If you are changing the injections, but want to re-use the ``FULL_DATA`` previous analysis, you can filter the ``main.map`` to keep the veto files, template bank files, full data inspiral files, and PSD files but filtering out any plots and result pages. For example::

    egrep '(VETO|BANK|INSPIRAL_FULL_DATA|MERGE_FULL_DATA|PSD)' /path/to/main.map | egrep -v '(png|html|dax)' > /path/to/reuse.map

* If you are changing the configuration of the coincident code, you can reuse all the injection files and inspiral files. For example::

    egrep '(VETO|BANK|FULL_DATA|PSD)' /path/to/main.map | egrep -v '(COINC|FIT|STATMAP|INJFIND|html|png|dax)' /path/to/main.map > /path/to/reuse.map

.. note::

    There is no rule for exactly which products can be reused as it depends on what you are changing in the workflow configuration. For partial reuse, it is best to consult an expert on how to build the cache file.

.. _weeklyahopeosg:

================================
Running on the Open Science Grid
================================

-------------
Prerequisites
-------------

There are a number of requirements on the machine on which the workflow will be started:

- Pegasus version 4.7.1 or later.

- The bundled executables available on the submit machine

- A gridftp server running on the submit machine

- Condor configured on the head node to connect to OSG as documented at::

    https://its-condor-blog.syr.edu/dokuwiki/doku.php?id=researchgroups:physics:sl7_cluster_setup

------------------------
Configuring the workflow
------------------------

.. note::

    Standard PyCBC installs build a version of ``pycbc_inspiral`` that uses weave
    to compile code at runtime. Many OSG machines do not have all of the
    compiler tools required to support weave compilation. In order to run on
    the OSG, the ``pycbc_inspiral`` executable must be built as a PyInstaller
    bundle that contains all of the weave-compiled code inside the bundle.
    This bundle must be built on a lowest-common denominator platform so that
    the shared libraries that it needs at runtime (e.g. glibc) are available.
    RHEL6 (or a similar derivative) is a suitable platform. The ``/cvmfs``
    filesystem contains ``pycbc_inspiral`` bundles that are built on the
    ``x86_64_rhel_6`` platform and are suitable for use on the OSG. For
    instructions on how to build PyInstaller bundled executables, see the page :ref:`building_bundled_executables`.


In order for ``pycbc_inspiral`` to be sent to worker nodes it must be
available via a remote protocol, either http, gsiftp, or CVMFS. Releases of
pycbc are installed in CVMFS and the LDG head nodes run a gridftp server that
can serve your own development copy.  Specify this path when you run
``pycbc_make_coinc_search_workflow``. To run from a released bundle in CVMFS 
give the following argument to the ``--config-overrides`` option (changing the
path to point to the release that you want to use)::

    'executables:inspiral:/cvmfs/oasis.opensciencegrid.org/ligo/sw/pycbc/x86_64_rhel_6/bundle/v1.6.6/pycbc_inspiral'

If you are running your own build of ``pycbc_inspiral``, you will need to give
a path to a gsiftp URL and tell Pegasus that the executable is not installed
on the OSG with the two ``--config-overrides`` options::

    'executables:inspiral:gsiftp://server.hostname/path/to/pycbc_inspiral'

Make sure this executable is build following the instructions on the page :ref:`building_bundled_executables`.

Add the following to the list of ``--config-overrides`` when running ``pycbc_make_coinc_search_workflow`` to tell Pegasus to run the inspiral code on the OSG::
     
    'pegasus_profile-inspiral:pycbc|site:osg'
    'pegasus_profile-inspiral:hints|execution.site:osg'
    'pegasus_profile-inspiral:pycbc|installed:False'
    'inspiral:fixed-weave-cache'

You also need a ``--config-overrides`` to ``pycbc_make_coinc_search_workflow`` that sets the staging site for the main workflow to the local site. To do this, add the following argument, replacing ``${WORKFLOW_NAME}`` with the string that is given as the argument to the option ``--workflow-name``::

    'workflow-${WORKFLOW_NAME}-main:staging-site:osg=local'

Optionally, you can add a configuration that will check that your grid proxy
is valid locally before submitting the job. This means that if your grid proxy
expires before the workflow is complete, the failure will be on the local site
before the job is actually submitted, and not on the remote site once the job
has been scheduled and matched::

    'pegasus_profile-inspiral:dagman|pre:/usr/bin/grid-proxy-info'

Another useful enhancement for OSG running is to add profiles to your inspiral
job that will tell Condor to put it on hold if it has been running for more
that 48 hours and terminate it after 5 failed attempts. To do this, add the
follwing lines to your ``executables.ini`` file::

    [pegasus_profile-inspiral]
    condor|periodic_hold = (JobStatus == 2) && ((CurrentTime - EnteredCurrentStatus) > (2 * 86400))
    condor|periodic_release = (JobStatus == 5) && (HoldReasonCode == 3) && (NumJobStarts < 5) && ((CurrentTime - EnteredCurrentStatus) > (300))
    condor|periodic_remove = (NumJobStarts >= 5)

--------------------
Running the workflow
--------------------

Add the following arguments to ``pycbc_submit_dax``::

    --no-create-proxy \
    --execution-sites osg \
    --append-pegasus-property 'pegasus.transfer.bypass.input.staging=true' \
    --remote-staging-server `hostname -f` \

``hostname -f`` will give the correct value if there is a gsiftp server running on the submit machine.  If not, change this as needed. The remote-staging-site is the intermediary computer than can pass files between the submitting computer and the computers doing the work.  ``hostname -f`` returns the full name of the computer. The full name of the computer that ``hostname -f`` has to be one that is accessible to both the submit machine and the workers. 

