##################################################
``pycbc_generate_mock_data``: Generating mock data
##################################################

``pycbc_generate_mock_data`` is a command-line tool which can be used for
generating mock data for a network of detectors. It also has a functionality
to choose from a population model, merger rate, and glitch distribution.
The basic idea of the tool is to do the following operations:

1. Generate noise realizations for a given network of detectors.
2. Inject signal population in the given noise realizations.
3. Provide zero noise realizations of signal and glitch population.

The following sections show concrete examples of common use cases.  For a
detailed list of all options, and option-specific help, see the output of::

    pycbc_generate_mock_data --help

=========================================================
Generate mock data realization for ground based detectors
=========================================================

The following command generate 16 second of data for a given noise realization 
without any signal or glitch::

    pycbc_generate_mock_data \
        --ifo-list H1 L1 V1 \
        --low-frequency-cutoff H1:20 L1:20 V1:20 \
        --psd-model H1:aLIGOZeroDetHighPower \
                    L1:aLIGOZeroDetHighPower \
                    V1:aLIGOZeroDetHighPower \
        --fake-strain-seed H1:1234 L1:2345 V1:3456 \
        --gps-start-time 1257294808 \
        --gps-end-time 1257294908 \
        --sample-rate 2048 \
        --channel-name H1:SIMULATED_STRAIN \
                       L1:SIMULATED_STRAIN \
                       V1:SIMULATED_STRAIN \
        --tag HLV_NOISE

This will produce simulated strain data files for each interferometer, tagged
with ``HLV_NOISE`` and stored in the default output directory. The GPS start
and end times determine the duration of the data segment, while the
``--fake-strain-seed`` option ensures that each detector's noise is drawn from
an independent random realization.

To generate time series data for a detector network, that contain the signal(s)
and noise realization, one need to provide an ``hdf`` injection file containing
the list of injection parameters. The following command generate the projected
signals in each detector of the network, without the noise::

     pycbc_generate_mock_data \
        --ifo-list H1 L1 V1 \
        --low-frequency-cutoff H1:20 L1:20 V1:20 \
        --psd-model H1:aLIGOZeroDetHighPower \
                    L1:aLIGOZeroDetHighPower \
                    V1:aLIGOZeroDetHighPower \
        --fake-strain-seed H1:1234 L1:2345 V1:3456 \
        --gps-start-time 1257294808 \
        --gps-end-time 1257294908 \
        --sample-rate 2048 \
        --channel-name H1:SIMULATED_STRAIN \
                       L1:SIMULATED_STRAIN \
                       V1:SIMULATED_STRAIN \
        --tag HLV_NOISE \
        --injection-file injection_10inj.hdf

The following command generate the projected signals in each detector of the 
network, without the noise::

     pycbc_generate_mock_data \
        --ifo-list H1 L1 V1 \
        --low-frequency-cutoff H1:20 L1:20 V1:20 \
        --psd-model H1:zeroNoise \
                    L1:zeroNoise \
                    V1:zeroNoise \
        --fake-strain-seed H1:1234 L1:2345 V1:3456 \
        --gps-start-time 1257294808 \
        --gps-end-time 1257294908 \
        --sample-rate 2048 \
        --channel-name H1:SIMULATED_STRAIN \
                       L1:SIMULATED_STRAIN \
                       V1:SIMULATED_STRAIN \
        --tag HLV_NOISE \
        --injection-file injection_10inj.hdf

=================================
Example Generating LISA mock data
=================================
In order to generate the LISA mock data, one need to provide additional options
related to LISA. The following command generate the LISA noise (without 
signals) for the given PSD::

    pycbc_generate_mock_data \
        --ifo-list LISA_A LISA_E LISA_T \
        --low-frequency-cutoff LISA_A:1e-4 LISA_E:1e-4 LISA_T:1e-4 \
        --psd-model LISA_A:analytical_psd_lisa_tdi_AE \
                    LISA_E:analytical_psd_lisa_tdi_AE \
                    LISA_T:analytical_psd_lisa_tdi_T \
        --fake-strain-seed LISA_A:100 LISA_E:150 LISA_T:200 \
        --gps-start-time 2121224.274911478 \
        --gps-end-time 4799924.274911478 \
        --sample-rate 0.2 \
        --len-arm 2.5e9 \
        --acc-noise-level 2.4e-15 \
        --oms-noise-level 7.9e-12 \
        --tdi '2.0' \
        --channel-name LISA_A:SIMULATED_STRAIN \
                    LISA_E:SIMULATED_STRAIN \
                    LISA_T:SIMULATED_STRAIN \
        --tag LISA_NOISE

To generate data containing signal+noise, one need to provide an injection 
``hdf`` file. This file may contain one or more injections. The following 
command generate the LISA mock data with signals and noise::

    pcbc_generate_mock_data \\
        --ifo-list LISA_A LISA_E LISA_T \
        --low-frequency-cutoff LISA_A:1e-4 LISA_E:1e-4 LISA_T:1e-4 \
        --psd-model LISA_A:analytical_psd_lisa_tdi_AE \
                    LISA_E:analytical_psd_lisa_tdi_AE \
                    LISA_T:analytical_psd_lisa_tdi_T \
        --fake-strain-seed LISA_A:100 LISA_E:150 LISA_T:200 \
        --gps-start-time 2121224.274911478 \
        --gps-end-time 4799924.274911478 \
        --sample-rate 0.2 \
        --len-arm 2.5e9 \
        --acc-noise-level 2.4e-15 \
        --oms-noise-level 7.9e-12 \
        --tdi '2.0' \
        --channel-name LISA_A:SIMULATED_STRAIN \
                       LISA_E:SIMULATED_STRAIN \
                       LISA_T:SIMULATED_STRAIN \
        --injection-file injections_lisa.hdf \
        --tag LISA_NOISE_PLUS_SIGNAL


The following command generate the LISA signals without noise::

    pcbc_generate_mock_data \
        --ifo-list LISA_A LISA_E LISA_T \
        --low-frequency-cutoff LISA_A:1e-4 LISA_E:1e-4 LISA_T:1e-4 \
        --psd-model LISA_A:zeroNoise LISA_E:zeroNoise LISA_T:zeroNoise \
        --fake-strain-seed LISA_A:100 LISA_E:150 LISA_T:200 \
        --gps-start-time 2121224.274911478 \
        --gps-end-time 4799924.274911478 \
        --sample-rate 0.2 \
        --len-arm 2.5e9 \
        --acc-noise-level 2.4e-15 \
        --oms-noise-level 7.9e-12 \
        --tdi '2.0' \
        --channel-name LISA_A:SIMULATED_STRAIN \
                       LISA_E:SIMULATED_STRAIN \
                       LISA_T:SIMULATED_STRAIN \
        --injection-file injections_lisa.hdf \
        --tag LISA_ZERONOISE_SIGNAL 

