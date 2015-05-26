=====
Instructions for running pycbc on sugar or atlas cluster
=====

Using the include environment install script :download:`example script <resources/install_pycbc_env>`

::

  chmod 755 install
  ./install_pycbc_env /abs/path/to/working/directory

After this, pycbc and its dependencies should be installed. You just need to source it

::

  source /path/to/dir/source

You may need to make sure you activate mkl (math kernel library) before you get going.

::

  source /opt/intel/bin/compilervars.sh intel64

Now if you have a script (which includes your start/end time, the parameters you need for pycbc_make_hdf_coinc_workflow) and a config file you are ready to set up the dax and submit it. See an example of both and explanations below:

::

  ./run.sh
  pycbc_submit_dax test.dax

run.sh
-----

``GPS_START_TIME=1102291216``
start time of your analysis - this does not not need to be a time when the detector is locked.
``GPS_END_TIME=1102489216``
end time of your analysis

``pycbc_make_hdf_coinc_workflow \``

executable

``--output-dir er6_test \``

Name of output directory - what ever you want

``--workflow-name er6_test \``

Name of the dax file - what ever you want

``--local-config-files er6_bns.ini \``

Name of the config file

``--config-overrides \``

Options to be added to the config file 

``workflow:start-time:${GPS_START_TIME} \``

Add start time of analysis to workflow module

``workflow:end-time:${GPS_END_TIME} \``

Add end time of the analysis to workflow module

``results_page:output-path:/home/laura.nuttall//public_html/cbc/pycbc_test/ER6``

Add output-path to the results_page module. This is where all the plots will be made


config file
-----

This config file is made up of 3 sections: the workflow, the pegasus profile and the options. The workflow controls how different parts of the the workflow hang together. The pegasus profile section are equivalent to lines you would have in a condor_submit file (e.g. requirements, storage size etc). Anything you would do in condor you would do here. The third section maps the options to the executables.

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
  ; PIPEDOWN demands we use AT_RUNTIME
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
- ``‘segments-science-veto’`` species which category of veto you want to eliminate from your search before it is performed to consider the data science. In this instance, 1 denotes that all the times of Cat 1 vetoes. 
- ``‘segments-veto-groups’`` is an option you can populate with different veto categories and diagnostic plots will be made after each veto is employed. 
- ``‘segments-final-veto-group’`` is an important option as the vetoes defined here will be removed from the search before it is performed. An option of 1 will remove all Cat 1 veto times from the analysis before it is performed. If you want to add cat 2 then the option is 12.

::

  [workflow-tmpltbank]
  ; See https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/ahope/template_bank.html
  tmpltbank-method=PREGENERATED_BANK
  tmpltbank-pregenerated-bank=/home/jveitch/projects/mdc/spin/tmpltbanks/nonspin/BNS_NonSpin_30Hz_earlyaLIGO.xml
  ; Remove the option below to disable linking with matchedfilter_utils

This section specifies which template bank to use

- ``’tmpltbank-method’`` option specifies whether you want to use a regenerated bank or to make it on the fly. In O1 we will be us a pregenerated bank. 
- ``‘tmpltbank-pregnerated-bank’`` specifies the location of the xml with the pregenerated bank

::

  [workflow-splittable]
  splittable-method = IN_WORKFLOW
  splittable-num-banks = 10

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

This section deals with software injections. Here you are specifying where to perform them. In this instance it is in the workflow

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
- ``l-distr`` = the longitudinal distribution of the injections (i.e. how spread out on the sky they are)
- ``i-distr`` = inclination of the injection 
- ``time-step`` = time between injections. This can be whatever time you want, but remember if the injections are too close together you can screw up your psd estimation. ~90s seems ok. 
- ``time-interval`` = time interval to inject the signal. It will not always be exactly at time-step, but at a time of time-step +/- random_number(0,time-interval)
- ``seed`` = random seed, choose whatever number you want so you

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

- ``chisq-bins`` = number of bins to chop the signal (in frequency) into. Each bin has equal weighting and leads to newSNR
- ``snr-threshold`` = newSNR threshold
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

This section concerns time slides without injections, and its purpose is to keep a small number of timesmlide triggers for background estimation. Time slides are done in a fixed window, which is defined by this bottom option of 1.1 seconds. We don’t store all the coincident triggers due to the time slides. We keep 200 of the loudest triggers for each template time slide, given by the second option, which gives a good estimation of the background at low FAR. The top option specifies which trigger number to keep from each template time slide to get an overall estimation of background (not just the loudest). In this instance we would keep the 1000th, 2000th, 3000th trigger etc.

::

  [coinc-injfull&coinc-fullinj]
  timeslide-interval=1.1
  loudest-keep-value = 8.5
  cluster-window = 10.0

This section concerns time slides with injections in the data. We assume only one injection will be coincident with a timeslide (done every 1.1 seconds - see first option)trigger and we keep it if its newSNR>8.5 as specified in the second option. 
Cluster window indicates the time window that we want to cluster over in the SNR time series of a template

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
