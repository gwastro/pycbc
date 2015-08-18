################################################################
Hardware injection waveform generation (``pycbc_generate_hwinj``)
################################################################

=====================================
Introduction
=====================================

The executable ``pycbc_generate_hwinj`` is for generating waveforms to be
used for hardware injections.

The procedure for generating a waveform is

 * Generate the waveform using the parameters given on the command line.
 * Generate a PSD.
 * Calculate the time delay and antenna pattern for each detector.
 * Calculate sigma for each detector using the waveform and the PSD. Determine the network SNR.
 * Rescale the waveform to the desired network SNR.
 * Write an XML file with the waveform parameters and an ASCII file with the h(t) timeseries.

This executable uses the ``pycbc.waveform`` module to generate the waveforms that uses the lalsimulation ``lalsimulation.SimInspiralChooseTDWaveform`` routine to return the two h(t) polarizations.

The command line options to ``pycbc_generate_hwinj`` are

 * ``--approximant`` the waveform approximant to use
 * ``--mass1`` and ``--mass2`` the component masses of the binary
 * ``--gps-end-time`` the geocentric end time of the hardware injection
 * ``--network-snr`` the network SNR of the hardware injection, network SNR is 
 * ``--ifos`` a list of the interferometers to include in the hardware injection
 * ``--output-waveform-file`` a list of ASCII filenames that will have a single-column with the h(t) timeseries, the ordering of the IFOs must match `--ifos`
 * ``--output-xml-file`` an XML filename with a ``sim_inspiral`` and ``sngl_inspiral`` table for the injection
 * ``--psd-model`` or ``--psd-estimate`` if you want to use a model noise curve or estimate the noise curve from data, if you estimate the PSD then you will need to include the additional options found in the PSD option group (see ``pycbc_generate_hwinj --help``)

==================================
How to generate a single waveform
==================================

Here is a usage example for generating a CBC waveform that could
be used as a hardware injection. ::

    pycbc_generate_hwinj --approximant EOBNRv2 \
                         --mass1 1.4 \
                         --mass2 1.4 \
                         --output-xml-file injection.xml.gz \
                         --output-waveform-file h1_injection.txt l1_injection.txt \
                         --gps-end-time 1117982141 \
                         --ifos H1 L1 \
                         --psd-model aLIGOZeroDetHighPower \
                         --psd-output injection_psd.txt \
                         --inclination 0.0 \
                         --network-snr 8

You can generate a CBC waveform using data to estimate the PSD instead of a model.
Here is a usage example. ::

    pycbc_generate_hwinj --approximant EOBNRv2 \
                         --mass1 1.4 \
                         --mass2 1.4 \
                         --output-xml-file injection.xml.gz \
                         --output-waveform-file h1_injection.txt l1_injection.txt \
                         --gps-end-time 1117982141 \
                         --ifos H1 L1 \
                         --psd-model aLIGOZeroDetHighPower \
                         --psd-output injection_psd.txt \
                         --inclination 0.0 \
                         --network-snr 8
