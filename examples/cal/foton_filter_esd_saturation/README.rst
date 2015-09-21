#################################################
ESD Saturation Example
#################################################

.. contents::

=================================================
Introduction
=================================================

Here we walkthrough an example of using the foton filtering module in pycbc. This example is contained in the following example directory::

  examples/cal/foton_filter_esd_saturation

We will generate h(t) for a CBC waveform and filter it will several filterbanks to get ETMY DAC counts.

This example should be run at either the LLO or LHO cluster.

=================================================
Add ROOT and foton to pycbc installation
=================================================

Do the following to add both the ROOT and foton python packages to the virtual environment ::

  NAME=/path/to/virtualenv
  cd ${NAME}/lib64/python2.6/site-packages
  ln -s /usr/lib64/python2.6/site-packages/libPyROOT.so
  ln -s /usr/lib64/python2.6/site-packages/ROOT.py
  ln -s /usr/lib64/python2.6/site-packages/ROOTwriter.py
  cd ${NAME}/lib/python2.6/site-packages
  ln -s /usr/lib/python2.6/site-packages/foton.py

=================================================
Set variables
=================================================

First change directory into the examples directory inside the cloned git repository ::

  cd examples/cal/foton_filter_esd_saturation

We select a time for generating the waveform (``${GEOCENT_END_TIME}``) and the time to use for PSD estimation when we generate the waveform. Here we select our reference time to be ::

  GEOCENT_END_TIME=1126257416
  GPS_START_TIME=1126256416
  GPS_END_TIME=$((${GPS_START_TIME} + 2048))

We will use the following frame type and channel to estimate the PSD when we generate the CBC waveform ::

  IFO=H1
  FRAME_TYPE=H1_HOFT_C00
  CHANNEL_NAME=GDS-CALIB_STRAIN

We will generate the CBC waveform at 16384Hz ::

  SAMPLE_RATE=16384

=================================================
Checkout the foton filter file SVN
=================================================

The foton filter files are checked into an SVN repository. You can check out the SVN repositories with the following commands ::

  svn co --username=albert.einstein https://daqsvn.ligo-la.caltech.edu/svn/h1_filter_files/h1_archive
  svn co --username=albert.einstein https://daqsvn.ligo-la.caltech.edu/svn/l1_filter_files/l1_archive

We will use filterbanks from the ``CALCS`` and ``SUSETMY`` models for H1. So set ::

  CALCS_FILTER_FILE=h1_archive/H1CALCS.txt
  SUS_FILTER_FILE=h1_archive/H1SUSETMY.txt

=================================================
Generate a CBC waveform
=================================================

Now we generate a CBC waveform using the hardware injection executable ::

  pycbc_generate_hwinj --instruments ${IFO} --waveform-low-frequency-cutoff 30 --geocentric-end-time ${GEOCENT_END_TIME} --gps-start-time ${GPS_START_TIME} --gps-end-time ${GPS_END_TIME} --frame-type ${IFO}:${FRAME_TYPE} --channel-name ${IFO}:${CHANNEL_NAME} --approximant SEOBNRv2 --order pseudoFourPN --mass1 26.6637001 --mass2 23.2229004 --inclination 1.04719755 --polarization 0.0 --ra 0.0 --dec 0.0 --taper TAPER_START --network-snr 18.424 --spin1z -0.963 --spin2z  -0.988 --psd-low-frequency-cutoff 40.0 --sample-rate ${IFO}:${SAMPLE_RATE} --pad-data 8 --strain-high-pass 30.0 --psd-estimation median --psd-segment-length 16 --psd-segment-stride 8
  
There are a number of command line options you can change. See the hardware injection documentation for more details.

=================================================
Run the ESD saturation script
=================================================

In the example directory (``examples/cal/foton_filter_esd_saturation``) there is a bash script called ``pycbc_check_esd_saturation.sh``. This script will filter the waveform will several filterbanks and a matrix element to get from h(t) to ETMY DAC counts

#. ``CAL-INJ_HARDWARE`` filterbank
#. ``SUS-ETMY_L3_LOCK_L`` filterbank
#. ``ETMY_L3_DRIVEALIGN_L2L`` filterbank
#. ``ETMY_L3_EUL2ESD_2_1`` matrix element
#. ``ETMY_L3_ESDOUTF_LL`` filterbank

The input to the code is the single-column ASCII file (``--data-file``), the frame type for reading ``SWSTAT`` and ``GAIN`` channels (``--frame-type``), the GPS times to check the ``SWSTAT`` and ``GAIN`` channels, and the foton filter files.

The GPS options should be recent times to get the current gains of the filterbank. For example ::

 GPS_START_TIME=$((`lalapps_tconvert` - 2000))
 GPS_END_TIME=$((${GPS_START_TIME}+1))

We will need to change the frame type to the frames that contains the ``SWSTAT`` and ``GAIN`` channels ::

  FRAME_TYPE=H1_R



An example command is ::

  sh pycbc_check_esd_saturation.sh --data-file ${PWD}/${IFO}-HWINJ_CBC-*-*.txt --gps-start-time ${GPS_START_TIME} --gps-end-time $((${GPS_START_TIME} + 1)) --ifo ${IFO} --frame-type ${FRAME_TYPE} --sus-filter-file ${PWD}/${SUS_FILTER_FILE} --calcs-filter-file ${PWD}/${CALCS_FILTER_FILE}

The output files will be written to ``${PWD}/esd_output/``. It will contain a single-column ASCII file after each stage of filtering. The final time series will contain ``*ETMY_L3_ESDOUTF_LL*``.

=================================================
Plotting the output
=================================================

In the example directory (``examples/cal/foton_filter_esd_saturation``) there is a script to use gwpy to make some plots of the output; the script is called ``gwpy_plot_hwinj``. The script ``gwpy_plot_hwinj`` will plot a timeseries and a spectrogram.

The inputs are arguments for x-axis minimum, x-axis maximum, time series y-axis minimum, time series y-axis maximum, spectrogram colorbar minimum, and colorbar maximum respectively.

In a new terminal source the gwpy environment ::

  source /home/detchar/opt/gwpysoft/etc/gwpy-user-env.sh

Make the output directory ::

  HTMLDIR=/home/${USER}/public_html/esd_test/
  mkdir -p ${HTMLDIR}

To plot the h(t) CBC waveform do ::

  INPUT_FILE=`ls ${IFO}-HWINJ_CBC-*-*.txt`
  TIMESERIES_FILE=${HTMLDIR}/${IFO}-TIMESERIES_HWINJ_CBC.png
  SPECTROGRAM_FILE=${HTMLDIR}/${IFO}-SPECTROGRAM_HWINJ_CBC.png
  python gwpy_plot_hwinj ${INPUT_FILE} ${TIMESERIES_FILE} ${SPECTROGRAM_FILE} 0 20 -2e-21 2e-21 1e-21 1e-27

To plot the ETMY DAC counts time series do ::

  INPUT_FILE=`ls esd_output/${IFO}-FILTER_ETMY_L3_ESDOUTF_LL-*.txt`
  TIMESERIES_FILE=${HTMLDIR}/${IFO}-TIMESERIES_ETMY_L3_ESDOUTF_LL.png
  SPECTROGRAM_FILE=${HTMLDIR}/${IFO}-SPECTROGRAM_ETMY_L3_ESDOUTF_LL.png
  python gwpy_plot_hwinj ${INPUT_FILE} ${TIMESERIES_FILE} ${SPECTROGRAM_FILE} 0 20 -3e4 3e4 1e1 1e-7

To plot the ETMY DAC counts at the merger do ::

  INPUT_FILE=`ls esd_output/${IFO}-FILTER_ETMY_L3_ESDOUTF_LL-*.txt`
  TIMESERIES_FILE=${HTMLDIR}/${IFO}-TIMESERIES_MERGER_ETMY_L3_ESDOUTF_LL.png
  SPECTROGRAM_FILE=${HTMLDIR}/${IFO}-SPECTROGRAM_MERGER_ETMY_L3_ESDOUTF_LL.png
  python gwpy_plot_hwinj ${INPUT_FILE} ${TIMESERIES_FILE} ${SPECTROGRAM_FILE} 5.9 6.1 -30000 30000 1e+1 1e-7
  
If you have an X11 session open then you can use the interactive hardware injection plotting code called ``pycbc_plot_hwinj``. To use this do ::

  INPUT_FILE=`ls esd_output/${IFO}-FILTER_ETMY_L3_ESDOUTF_LL-*.txt`
  pycbc_plot_hwinj ${INPUT_FILE}
