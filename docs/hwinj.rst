######################################
Hardware injection waveform generation
######################################

.. contents::

============
Introduction
============

This page describes how to generate waveforms and save them as single-column ASCII waveform files that can be used by ``awgstream`` to inject into the detector.

There are two executables that can be used to generate single-column ASCII files; they are ``pycbc_generate_hwinj`` and ``pycbc_generate_hwinj_from_xml``. Both executables use the PyCBC injection module (``pycbc.inject``) to inject the coherent waveform into a time series of zeroes.

The executable ``pycbc_generate_hwinj`` generates a waveform using parameters from the command line. The user inputs parameters such as ``--mass1``, ``--mass2``, etc. on the command line. This executable is useful for generating a specific coherent waveform for hardware injections.

The executable ``pycbc_generate_hwinj_from_xml`` generates all the waveforms in a LIGOLW ``sim_inspiral`` table. The output of ``lalapps_inspinj`` (an executable for generating a population of injections) is a LIGOLW ``sim_inspiral`` table. This executable is useful if you want to generate a population of coherent waveforms for hardware injections.

==============================================================
Generate waveform from command line (``pycbc_generate_hwinj``)
==============================================================

Here is a usage example for generating a CBC waveform using detector data with ``pycbc_generate_hwinj``.

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Select a time for the injection
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

First on the command line set a variable for the GPS geocentric end time of the coherent injection ::

  GEOCENT_END_TIME=1124381661

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Select data for PSD estimation
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

We need to set what data we will use to estimate the PSD. This should be a 2048 second interval. So in this example we will set ::

  GPS_START_TIME=1124380361
  GPS_END_TIME=1124382409

We will want to use data from when the detector was in science mode time to estimate the PSD. To check if the detector was in science mode you can query the segment database, check the section :ref:`howtoquerysegdb` for an example command.

Since we are using detector data we will need to read the data from frame files. We will need to specify the frame type and channel name. The frame type is a way to identify what list of channels is in the frame file. So on the command line we set ::

  FRAME_TYPE=H1_HOFT_C00
  CHANNEL_NAME=H1:GDS-CALIB_STRAIN

We can check that frames exist for this time by querying the LDR server, check the section :ref:`howtoqueryldr` for an example command.

**Do not assume that the frame type and channel name are the same as in this section. These values are correct for this example.**

&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Run ``pycbc_generate_hwinj``
&&&&&&&&&&&&&&&&&&&&&&&&&&&&

The ``--instruments`` option specifies the IFOs to include in the injection. The executable ``pycbc_generate_hwinj`` will write the waveform for all detectors listed with the ``--instruments`` option. There are a number of options, including ``--instruments`` that take a space-seperated list, example usage of this option is ::

    --instruments H1 L1


Now let's specify the parameters of the waveform. A full list of the command line options is available with ``pycbc_generate_hwinj --help``. Here we provide an example for generating a coherent 1.4-1.4 component mass binary using the ``EOBNRv2`` approximant. Here is an example of this command line case ::

    --approximant SEOBNRv2 --order pseudoFourPN --mass1 25.0 --mass2 25.0 --inclination 0.0 --polarization 0.0 --ra 0.0 --dec 0.0 

The sample rate of the output ASCII file is given by the ``--sample-rate`` option. The sample rate for each IFO must match, example usage of this option is ::

    --sample-rate H1:16384 L1:16384

We have the option to query the LDR server directly with ``pycbc_generate_hwinj``. So we only need to specify a ``--frame-type`` option to get the data. Alternatively you an pass a space-seperated list with ``--frame-files`` or a LAL frame cache with ``--frame-cache``. We also need to say what channel to use for the PSD estimation with the ``--channel-name`` option. Example usage of ``--frame-type`` would be ::

    --frame-type H1:${FRAME_TYPE} L1:${FRAME_TYPE} --channel-name H1:${CHANNEL_NAME} L1:${CHANNEL_NAME}

We do not want to inject a step-like response into the detector, therefore we taper the waveform at the beginning. The ``EOBNRv2`` has a ringdown at the end so we do not need to taper the end. Example usage of this option is ::

    --taper TAPER_START

We specify the network SNR we want the coherent injection to have on the command line. The network SNR calculation includes all IFOs specified in the ``--insturments`` option, example usage of this option is ::

    --network-snr 28

In calculating the network SNR the executable ``pycbc_generate_hwinj`` will generate a PSD and calculate an SNR for the waveform. The options ``--psd-low-frequency-cutoff`` and ``--psd-high-frequency-cutoff`` set the min and max frequency for the SNR calculation. The waveform used in the SNR calculation is also generated at this low-frequency cutoff, note the waveform is not written to disk with this low-frequency cutoff. Example usage of the PSD options is ::

    --psd-low-frequency-cutoff 40.0 --psd-high-frequency-cutoff 1000.0 --psd-estimation median --psd-segment-length 16 --psd-segment-stride 8 --pad-data 8

The additional PSD options dictate how the PSD will be calculated, ie. how many seconds per FFT and how much overlap. The ``--pad-data`` option is how much data to disgard at the edges of our time series used in PSD estimation to avoid data corruption.

The ``--waveform-low-frequency-cutoff`` option is the frequency that ``pycbc_generate_hwinj`` will begin generating the waveform that is written to file.

Here is a full example command for generating an injection in only H1 ::

  pycbc_generate_hwinj --psd-high-frequency-cutoff 1000.0 --geocentric-end-time ${GEOCENT_END_TIME} --gps-start-time H1:${GPS_START_TIME} --gps-end-time H1:${GPS_END_TIME} --frame-type H1:${FRAME_TYPE} --channel-name H1:${CHANNEL_NAME} --approximant SEOBNRv2 --order pseudoFourPN --mass1 25.0 --mass2 25.0 --inclination 0.0 --polarization 0.0 --ra 0.0 --dec 0.0 --taper TAPER_START --network-snr 28 --waveform-low-frequency-cutoff 10.0 --psd-low-frequency-cutoff 40.0 --sample-rate 16384 --pad-data 8 --strain-high-pass 30.0 --psd-estimation median --psd-segment-length 16 --psd-segment-stride 8 --instruments H1

This will generate a single-column ASCII files that contains the h(t) time series for each detector and a LIGOLW XML file with the waveform parameters. The output filenames are not specified on the command line, they are determined internally by ``pycbc_generate_hwinj``. In this example the ASCII file with the waveform will be named ``L1-HWINJ_CBC-${START}-${DURATION}.txt`` where ``${START}`` is the start time stamp of the time series and ``${DURATION}`` is the length in seconds of the ASCII waveform file. The LIGOLW XML file will be named ``H1L1-HWINJ_CBC-${START}-${DURATION}.xml.gz``.

The LIGOLW XML file contains a ``process_params`` table that saves the command line that was used to generate the waveform for future reference. It also includes a ``sim_inspiral`` table and a ``sngl_inspiral`` table. The ``sim_inspiral`` entry allows us to use the parameters of the waveform as a software injection in the PyCBC matched filtering executable ``pycbc_inspiral``. The ``sngl_inspiral`` entry allows us to use the parameters of the waveform as the filter in ``pycbc_inspiral``.

The user should inspect the waveforms. For a waveform plotting executable see section :ref:`runpycbcplothwinj`.

=====================================================================================
Generate waveform from ``lalapps_inspinj`` output (``pycbc_generate_hwinj_from_xml``)
=====================================================================================

Here is a usage case for generating a population of waveforms with ``lalapps_inspinj``. This example generates an injection every Tuesday for three months.

&&&&&&&&&&&&&&&&&&&&&&&
Run ``lalapps_inspinj``
&&&&&&&&&&&&&&&&&&&&&&&

Here we show an example on how to use ``lalapps_inspinj`` to generate a population of injections.

In this example we will select distributions for time, distance, inclination, mass, spin, and sky location. Below is a explaination of the command line options in this example. A full list of command line options can be found with ``lalapps_inspinj --help``.

Our time distribution will be a fixed time step to perform an injection every Tuesday. We will allow the injections to be anytime during the day (86400 seconds) and have a minimum one week (604800 seconds) between injections. The command line options will be ::

  --time-interval 86400 --time-step 604800 --gps-start-time 1126368017 --gps-end-time 1130371217

Our distance distribution will be uniformly distributed in volume. We set the minimum and maximum chirp distance in units of kpc. The command line options will be ::

  --d-distr volume --min-distance 10000 --max-distance 40000

Our mass distribution will be uniform in total mass. We can select the minimum and maximum component masses in units of solar masses. The command line options will be ::

  --m-distr totalMass --min-mass1 1.0 --max-mass1 2.0 --min-mass2 1.0 --max-mass2 2.0

We can select the minimum and maximum component spins. The command line options will be ::

 --enable-spin --min-spin1 0.0 --max-spin1 0.04 --min-spin2 0.0 --max-spin2 0.04

Our inclination distribution will be uniform. The command line option will be ::

  --i-distr uniform

Our source distribution will be random. The command line option will be ::
  
  --l-distr random
  
We select to use the ``SpinTaylorT4`` approximant and begin the waveforms at 10.0Hz. Here we taper the injection at the start and end of the injection. The command line options will be ::

  --waveform SpinTaylorT4threePointFivePN --f-lower 10 --taper-injection startend --band-pass-injection

Now we can combine all the options above and run ``lalapps_inspinj`` as ::

  lalapps_inspinj --time-interval 86400 --time-step 604800 --gps-start-time 1126368017 --gps-end-time 1130371217 --d-distr volume --min-distance 10000 --max-distance 40000 --m-distr totalMass --min-mass1 1.0 --max-mass1 2.0 --min-mass2 1.0 --max-mass2 2.0 --enable-spin --min-spin1 0.0 --max-spin1 0.04 --min-spin2 0.0 --max-spin2 0.04 --i-distr uniform --l-distr random --waveform SpinTaylorT4threePointFivePN --f-lower 10 --taper-injection startend --band-pass-injection

In this example ``lalapps_inspinj`` will write a LIGOLW XML file called ``HL-INJECTIONS_1-1126368017-4003200.xml`` that has a ``sim_inspiral`` table with the population of injections.

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Run ``pycbc_generate_hwinj_from_xml``
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

Running ``lalapps_inspinj`` has written a LIGOLW XML file with a ``sim_inspiral`` table. Now we can run ``pycbc_generate_hwinj_from_xml`` to write single-column ASCII waveform files for the population of injections.

There are just two command line options ``--injection-file`` (path to the LIGOLW XML file that ``lalapps_inspinj`` had written) and ``--sample-rate`` (the sample rate of the ASCII waveform files).

In this example we set the sample rate to 16384Hz so on the command line do ::

  pycbc_generate_hwinj_from_xml --injection-file HL-INJECTIONS_1-1126368017-4003200.xml --sample-rate 16384

As this command runs it will generate a H1 and L1 ASCII waveform file for each row in the ``sim_inspiral`` table.

The ASCII waveform files will be named ``${IFO}-HWINJ_CBC_SIMULATION_ID_${SIMID}-${START}-${DURATION}.txt`` where where ``${SIMID}`` is the ``simulation_id`` number for the ``sim_inspiral`` row, ``${START}`` is the GPS start time of the ASCII waveform file, and ``${DURATION}`` is the duration of the file in seconds.

The user should inspect the waveforms. For a waveform plotting executable see section :ref:`runpycbcplothwinj`.

========================================
Checks for the hardware injection output
========================================

Here are some follow-up checks the user can do.

.. _runpycbcplothwinj:

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Plot ASCII waveform files with ``pycbc_plot_hwinj``
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

You can plot the ASCII waveform files with an X11 connection. It's strongly recommended to use the X11 connection instead of saving a static image of the entire waveform. The X11 connection allows the user to zoom in and inspect the waveform more closely. A basic inspection would include checking the amplitude, the tapering, and the ringdown  of the waveforms are reasonable. For the ``pycbc_generate_hwinj`` example above one would do ::

  pycbc_plot_hwinj L1-HWINJ_CBC-${START}-${DURATION}.txt

If you are using ``ssh`` or ``gsissh`` to log into a cluster, you can provide the ``-Y`` option to open an X11 connection. For example ::

  gsissh -Y ldas-pcdev1.ligo.caltech.edu

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Recover software injection with ``pycbc_inspiral``
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

The executable ``pycbc_generate_hwinj`` will create an XML file with both a ``sim_inspiral`` and ``sngl_inspiral`` table. Therefore we can inject the exact waveform parameters and recover them with the exact template.

The analogous software injection command for the example above would be ::

  TMPLTBANK_FILE=H1-HWINJ_CBC-${START}-${DURATION}.xml.gz
  INSPIRAL_FILE=H1-INSPIRAL_PYCBC-${GPS_START_TIME}-$((${GPS_END_TIME}-${GPS_START_TIME})).xml.gz
  pycbc_inspiral --segment-end-pad 64  --segment-length 256 --segment-start-pad 64 --psd-estimation median --psd-segment-length 16 --psd-segment-stride 8 --psd-inverse-length 16 --pad-data 8 --sample-rate 4096 --low-frequency-cutoff 40 --strain-high-pass 30 --filter-inj-only --processing-scheme cpu --cluster-method template --approximant SEOBNRv2 --order 8 --snr-threshold 5.5 --chisq-bins 16 --channel-name ${CHANNEL_NAME} --gps-start-time ${GPS_START_TIME} --gps-end-time ${GPS_END_TIME} --trig-start-time $(($GEOCENT_END_TIME - 2)) --trig-end-time $(($GEOCENT_END_TIME + 2)) --frame-type ${FRAME_TYPE} --injection-file ${TMPLTBANK_FILE}  --bank-file ${TMPLTBANK_FILE} --output ${INSPIRAL_FILE} --verbose

Where ``${START}`` is the start of the injection and ``${DURATION}`` is the length of the injection. We kept the same PSD options (eg. ``--psd-segment-length``, etc.), data, high-pass filter, and low-frequency-cutoff.

You can print out the recovered SNR and other parameters with ``lwtprint``, for example ::

  lwtprint -t sngl_inspiral -c end_time,snr ${INSPIRAL_FILE}

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Recover ASCII file injection with ``pycbc_inspiral``
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

There is an executable ``pycbc_insert_frame_hwinj`` that will read the single-column ASCII file and insert it into frame data. An example command is here ::

  HWINJ_FILE=H1-HWINJ_CBC-${START}-${DURATION}.txt
  pycbc_insert_frame_hwinj --frame-type ${FRAME_TYPE} --channel-name H1:${CHANNEL_NAME} --gps-start-time $((${GPS_START_TIME} - 16)) --gps-end-time $((${GPS_END_TIME} + 16)) --pad-data 8 --strain-high-pass 30.0 --sample-rate 16384 --hwinj-file ${HWINJ_FILE} --hwinj-start-time ${START} --ifo H1 --output-file H1-HWINJ.gwf

Where ``${START}`` is the start of the injection and ``${DURATION}`` is the length of the injection.

Then you can run pycbc on the output frame file ``H1-HWINJ.gwf``.

.. _howtoquerysegdb:

---------------------------------
How to query the segment database
---------------------------------

Here is an example on how to check if the detector was in science mode for a GPS time interval. To do this we query the segment database. A command line tool to do check the ``pycbc_generate_hwinj`` example above is ::

  ligolw_segment_query_dqsegdb --query-segments --segment-url https://dqsegdb5.phy.syr.edu --gps-start-time 1124380361 --gps-end-time 1124382409 --include-segments L1:DMT-ANALYSIS_READY:1 --output-file L1-SEGMENTS.xml
  ligolw_segment_query_dqsegdb --query-segments --segment-url https://dqsegdb5.phy.syr.edu --gps-start-time 1124380361 --gps-end-time 1124382409 --include-segments H1:DMT-ANALYSIS_READY:1 --output-file H1-SEGMENTS.xml

This should write two XML files ``L1-SEGMENTS.xml`` and another ``H1-SEGMENTS.xml``. You can check the ``segment`` table to see if the detector was in science mode for this time. A command line tool that helps is ::

  ligolw_print --table segment --column start_time --column end_time L1-SEGMENTS.xml
  ligolw_print --table segment --column start_time --column end_time H1-SEGMENTS.xml

The output should be ``1124380361,1124382409`` for both ``ligolw_print`` commands. This tells us that the detector was in science mode for the entire time since there is one segment that is the equal to the interval of ``--gps-start-time`` to ``--gps-end-time``.

**Do not assume that the segment databse URL and science-mode segment names are the same as in this section. These values are correct for this example.**

.. _howtoqueryldr:

---------------------------
How to query the LDR server
---------------------------

Here is an example on how to check if frame files exist for a GPS time interval. To do this we query the LDR server. A command line tool to do check the ``pycbc_generate_hwinj`` example above is ::

  gw_data_find --observatory L --type L1_RDS --gps-start-time 1124380361  --gps-end-time 1124382409 --url-type file --gaps
  gw_data_find --observatory H --type H1_RDS --gps-start-time 1124380361  --gps-end-time 1124382409 --url-type file --gaps

A list of frame files will be print to your terminal if they are accessible. If ``Missing segments`` is printed, then you will not be able to access all the frame files.

**Do not assume that the frame type and channel name are the same as in this section. These values are correct for this example.**
