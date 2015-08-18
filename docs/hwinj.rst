################################################################
Hardware injection waveform generation (``pycbc_generate_hwinj``)
################################################################

=====================================
Introduction
=====================================

The executable ``pycbc_generate_hwinj`` is for generating waveforms to be used for hardware injections.

The procedure for generating a waveform is:

 * Generate the waveform using the parameters given on the command line.
 * Generate a PSD that will be used for calculating sigma squared.
 * Calculate the antenna pattern for each detector.
 * Calculate sigma squared for each detector using the waveform and the PSD. We do this to determine the network SNR.
 * Rescale the waveform to the desired network SNR.
 * Write an XML file with the waveform parameters and single-column ASCII files that contain the h(t) timeseries.

This executable uses the ``pycbc.waveform`` module to generate the waveforms. The ``pycbc.waveform`` module is calling lalsimulation routines to inject h(t).

The waveform-specific command line options to ``pycbc_generate_hwinj`` are:

 * ``--approximant`` the waveform approximant to use,
 * ``--order`` the post-Newtonian order of the waveform,
 * ``--inclination`` the inclination of the orbit,
 * ``--polarization`` the polarization,
 * ``--ra`` and ``--dec`` the RA and DEC of the waveform in radians,
 * ``--taper`` tells when to taper the waveform,
 * ``--mass1`` and ``--mass2`` the component masses of the binary,
 * ``--geocentric-end-time`` the geocentric end time of the hardware injection,
 * ``--network-snr`` the desired network SNR of the hardware injection,
 * ``--sample-rate`` the sample rate to generate the ASCII waveforms, and
 * ``--h1`` and ``--l1`` select what IFOs to generate the hardware injection.

==================================
How to generate a single waveform
==================================

Here is a usage example for generating a CBC waveform using a design noise curve. ::

  GPS_START_TIME=1117982000
  GPS_END_TIME=$(($GPS_START_TIME + 2048))

  pycbc_generate_hwinj --approximant EOBNRv2 --order pseudoFourPN --mass1 1.4 --mass2 1.4 --inclination 0.0 --polarization 0.0 --ra 1.0 --dec 1.0 --taper TAPER_START --network-snr 28 --geocentric-end-time 1117982241 --low-frequency-cutoff 15.0 --gps-start-time ${GPS_START_TIME} --gps-end-time ${GPS_END_TIME} --h1 --l1 --fake-strain aLIGOZeroDetHighPower --psd-model aLIGOZeroDetHighPower --sample-rate 16384


You can generate a CBC waveform using interferometer data to estimate the PSD instead of a model. Using the pycbc strain and psd options.
