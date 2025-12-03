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

==========================================================
Generate Mock Data Realizations for Ground-Based Detectors
==========================================================

This section demonstrates how to generate simulated strain data for a network
of ground-based gravitational-wave detectors using
:command:`pycbc_generate_mock_data`.

Mock data can include **noise only**, **signals only**, or **signals embedded
in noise**, depending on the options provided.

---------------------------------------------------------
Noise-only mock data
---------------------------------------------------------

The following command generates a **100-second noise realization** for the H1,
L1, and V1 interferometers, without any injected signals or glitches::

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

---------------------------------------------------------
Signals embedded in noise
---------------------------------------------------------

To generate time series data for a detector network, that contain the signal(s)
and noise realization, one need to provide an ``hdf`` injection file containing
the list of injection parameters. The following command generate the projected
signals in each detector of the network, embeded in the noise described by PSD
of each detectors::

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
        --tag HLV_NOISE_SIGNAL \
        --injection-file injection_hlv_10_bbh.hdf

---------------------------------------------------------
Signal-only (zero-noise) mock data
---------------------------------------------------------
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
        --tag HLV_SIGNAL \
        --injection-file injection_hlv_10_bbh.hdf

=================================
Example Generating LISA mock data
=================================

Generating mock data for LISA requires providing several LISA-specific options,
including the spacecraft arm length, optical metrology (OMS) noise, acceleration
noise, and the Time-Delay Interferometry (TDI) generation. The examples below
demonstrate how to generate noise-only datasets, signal+noise datasets, and
signal-only datasets for the LISA TDI channels A, E, and T.

---------------------------------------------------------
Noise-only LISA mock data
---------------------------------------------------------
The following command generates **noise-only** TDI data for the LISA A, E, and T
channels, using analytical PSD models for the corresponding TDI variables::

    pycbc_generate_mock_data \
        --ifo-list LISA_A LISA_E LISA_T \
        --low-frequency-cutoff LISA_A:1e-4 LISA_E:1e-4 LISA_T:1e-4 \
        --psd-model LISA_A:analytical_psd_lisa_tdi_AE \
                    LISA_E:analytical_psd_lisa_tdi_AE \
                    LISA_T:analytical_psd_lisa_tdi_T \
        --fake-strain-seed LISA_A:100 LISA_E:150 LISA_T:200 \
        --gps-start-time 0 \
        --gps-end-time 31536000 \
        --sample-rate 0.2 \
        --len-arm 2.5e9 \
        --acc-noise-level 2.4e-15 \
        --oms-noise-level 7.9e-12 \
        --tdi '2.0' \
        --channel-name LISA_A:SIMULATED_STRAIN \
                       LISA_E:SIMULATED_STRAIN \
                       LISA_T:SIMULATED_STRAIN \
        --tag LISA_NOISE \
        --fake-strain-filter-duration 31536000

---------------------------------------------------------
Signal + noise LISA mock data
---------------------------------------------------------

To generate data containing signal+noise, one need to provide an injection 
``hdf`` file. This file may contain one or more injections. The following 
command generate the LISA mock data with signals and noise::

    pycbc_generate_mock_data \
        --ifo-list LISA_A LISA_E LISA_T \
        --low-frequency-cutoff LISA_A:1e-4 LISA_E:1e-4 LISA_T:1e-4 \
        --psd-model LISA_A:analytical_psd_lisa_tdi_AE \
                    LISA_E:analytical_psd_lisa_tdi_AE \
                    LISA_T:analytical_psd_lisa_tdi_T \
        --fake-strain-seed LISA_A:100 LISA_E:150 LISA_T:200 \
        --gps-start-time 0 \
        --gps-end-time 31536000 \
        --sample-rate 0.2 \
        --len-arm 2.5e9 \
        --acc-noise-level 2.4e-15 \
        --oms-noise-level 7.9e-12 \
        --tdi '2.0' \
        --channel-name LISA_A:SIMULATED_STRAIN \
                       LISA_E:SIMULATED_STRAIN \
                       LISA_T:SIMULATED_STRAIN \
        --injection-file injection_lisa_1_smbhb.hdf \
        --tag LISA_NOISE_SIGNAL \
        --fake-strain-filter-duration 31536000


The command :command:`pycbc_create_injections` can be used to generate such an
injection configuration file. The example below shows how to create a file 
containing one LISA SMBHB injection:

.. code-block:: bash

    pycbc_create_injections --verbose \
        --config-files injection_lisa_1_smbhb.ini \
        --ninjections 1 \
        --seed 10 \
        --output-file injection_lisa_1_smbhb.hdf \
        --variable-params-section variable_params \
        --static-params-section static_params \
        --dist-section prior

The example ``injection_lisa_1_smbhb.ini`` file can be downloaded from
`examples/mdc_examples/injection_lisa.ini <examples/mdc_generation/injection_lisa_1_smbhb.ini>`_.

.. important::

   To generate LISA mock data, you must install the following packages in order
   inside your virtual environment:

   1. **BBHx**  
      https://github.com/mikekatz04/BBHx/tree/v1.0.5

   2. **BBHx waveform plugin**  
      https://github.com/gwastro/BBHX-waveform-model

---------------------------------------------------------
Signal-only (zero-noise) LISA mock data
---------------------------------------------------------
To generate **signal-only** mock data (no instrument noise), simply use the
``zeroNoise`` PSD::

    pycbc_generate_mock_data \
        --ifo-list LISA_A LISA_E LISA_T \
        --low-frequency-cutoff LISA_A:1e-4 LISA_E:1e-4 LISA_T:1e-4 \
        --psd-model LISA_A:zeroNoise LISA_E:zeroNoise LISA_T:zeroNoise \
        --fake-strain-seed LISA_A:100 LISA_E:150 LISA_T:200 \
        --gps-start-time 0 \
        --gps-end-time 31536000 \
        --sample-rate 0.2 \
        --len-arm 2.5e9 \
        --acc-noise-level 2.4e-15 \
        --oms-noise-level 7.9e-12 \
        --tdi '2.0' \
        --channel-name LISA_A:SIMULATED_STRAIN \
                       LISA_E:SIMULATED_STRAIN \
                       LISA_T:SIMULATED_STRAIN \
        --injection-file injection_lisa_1_smbhb.hdf \
        --tag LISA_ZERONOISE_SIGNAL \
        --fake-strain-filter-duration 31536000
