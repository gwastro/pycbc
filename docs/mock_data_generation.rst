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

=====================================================
Generate noise realization for a network of detectors
=====================================================

The following command generate 16 second of data for a given noise realiztion::

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


