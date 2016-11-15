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

    pycbc_make_hdf_coinc_workflow --workflow-name s6d_chunk3 --output-dir output \
      --config-files https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-config/download/master/S6/pipeline/s6_run_pycbc_er8_pre_release.ini \
      https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-config/download/master/S6/pipeline/executables.ini \
      https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-config/download/master/S6/pipeline/injections.ini \
      https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-config/download/master/S6/pipeline/data_S6.ini \
      https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-config/download/master/S6/pipeline/gps_times_s6d_big_dog_two_weeks.ini \
      --config-overrides "results_page:output-path:${HOME}/public_html/s6/s6d-big-dog-weeks"

The configuration ``results_page:output-path`` can be changed appropriately to
set the output web page location.

.. note::

   To use released exectutables for production analysis, you should specify
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
command line option to ``pycbc_make_hdf_coinc_workflow`` and plan and submit
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
called ``output.map`` in the directory that you specified with the
``--output`` argument to ``pycbc_make_coinc_search_workflow``. This file contains a 
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

    pycbc_submit_dax --accounting-group ligo.dev.o1.cbc.explore.test --dax s6d_chunk3.dax --cache /path/to/prior_data.map

Follow the remaining :ref:`coincworkflowplan` instructions to submit your reduced
workflow.


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

In order for ``pycbc_inspiral`` to be sent to worker nodes it must be available
via a remote protocol, either http or gsiftp. The LDG head nodes run a gridftp
server that should be able to serve this, or code can obtain the code from the
web server on ``code.pycbc.phy.syr.edu``. Edit your ``executables.ini`` to 
specify this path or give the path when you run
``pycbc_make_coinc_search_workflow``, for example, with the option::

    --config-overrides 'executables:inspiral:gsiftp://server.name/path/to/pycbc_inspiral'

Add the following to the list of ``--config-overrides`` when running ``pycbc_make_coinc_search_workflow``::
     
    'pegasus_profile-inspiral:pycbc|site:osg' \
    'pegasus_profile-inspiral:hints|execution.site:osg' \
    'pegasus_profile-inspiral:condor|request_memory:1920M' \

You also need a ``--config-overrides`` to ``pycbc_make_coinc_search_workflow`` that sets the staging site for the main workflow to the local site. To do this, add the following argument, replacing ``${WORKFLOW_NAME}`` with the string that is given as the argument to the option ``--workflow-name ``::

    'workflow-${WORKFLOW_NAME}-main:staging-site:osg=local' \

Optionally, you can add a configuration that will check that your grid proxy
is valid locally before submitting the job. This means that if your grid proxy
expires before the workflow is complete, the failure will be on the local site
before the job is actually submitted, and not on the remote site once the job
has been scheduled and matched::

    'pegasus_profile-inspiral:dagman|pre:/usr/bin/grid-proxy-info' \

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
    --cache [URL/location of osg cache file] \

.. note::
   Cache files for the C02 frames can be obtained from https://code.pycbc.phy.syr.edu/ligo-cbc/pycbc-config/tree/master/O1/osg-cache-frame-files. The cache file contains the location of the frame files on Sugar, XSede, and xroot.

``hostname -f`` will give the correct value if there is a gsiftp server running on the submit machine.  If not, change this as needed. The remote-staging-site is the intermediary computer than can pass files between the submitting computer and the computers doing the work.  ``hostname -f`` returns the full name of the computer. The full name of the computer that ``hostname -f`` has to be one that is accessible to both the submit machine and the workers. 

