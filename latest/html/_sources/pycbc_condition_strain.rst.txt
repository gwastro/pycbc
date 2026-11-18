#############################################################
``pycbc_condition_strain``: operations with strain data files
#############################################################

``pycbc_condition_strain`` is a command-line tool which can be used for a
variety of operations related to gravitational-wave strain data files.
The basic idea of the tool is to combine the following operations:

1. Obtain a segment of gravitational-wave strain data, either by simulating
   detector noise, or reading an existing data file, or downloading public
   data from `GWOSC`_.
2. Optionally, simulate one or more signals and add them to the data.
3. Optionally, apply a number of operations that a matched-filter search would
   apply before matched-filtering, such as low- or high-passing, resampling or
   glitch gating. These operations are sometimes referred to as "data
   conditioning", hence the name of the tool.
4. Write the resulting data back to one or more files. Different output file
   formats are supported.

The following sections show concrete examples of common use cases.  For a
detailed list of all options, and option-specific help, see the output of::

    pycbc_condition_strain --help

================================================
Splitting existing data into shorter frame files
================================================

The following command downloads a 64 s segment of gravitational-wave data from
`GWOSC`_ and writes the data as four frame files of 16 s duration::

    pycbc_condition_strain \
        --gps-start-time 1242442818 \
        --gps-end-time 1242442882 \
        --frame-type GWOSC \
        --channel-name H1:GWOSC-16KHZ_R1_STRAIN \
        --frame-duration 16 \
        --output-strain-file 'H1-STRAIN-{start}-{duration}.gwf'

Note that the strings ``{start}`` and ``{duration}`` in the output file name
should be literally given like that. They are automatically replaced with the
right values by ``pycbc_condition_strain`` when the ``--frame-duration``
option is used.

==========================
Conditioning existing data
==========================

This command downloads 64 s of gravitational-wave data from `GWOSC`_, applies
a high-pass filter with corner frequency of 15 Hz, downsamples the data to
2048 Hz, and writes it to a single frame file::

    pycbc_condition_strain \
        --gps-start-time 1242442818 \
        --gps-end-time 1242442882 \
        --frame-type GWOSC \
        --channel-name H1:GWOSC-16KHZ_R1_STRAIN \
        --strain-high-pass 15 \
        --sample-rate 2048 \
        --output-strain-file H1-STRAIN_CONDITIONED-1242442818-64.gwf

==============================================
Injecting simulated signals into existing data
==============================================

Read 64 s of existing gravitational-wave data from `GWOSC`_, add one or more
simulated signals into it (specified by the ``injections.hdf`` file), and
write the result to four frame files of 16 s duration::

    pycbc_condition_strain \
        --gps-start-time 1242442818 \
        --gps-end-time 1242442882 \
        --frame-type GWOSC \
        --channel-name H1:GWOSC-16KHZ_R1_STRAIN \
        --injection-file injections.hdf \
        --frame-duration 16 \
        --output-strain-file 'H1-STRAIN_WITH_INJECTIONS-{start}-{duration}.gwf'

==================================
Simulating gravitational-wave data
==================================

The following command simulates 64 s of detector noise with a given power
spectral density (in this case, a model for the noise expected for Virgo
in the O4 run), adds one or more simulated signals to it (specified by the
``injections.hdf`` file), and writes the result to a single frame file::

    pycbc_condition_strain \
        --fake-strain aLIGOAdVO4T1800545 \
        --fake-strain-seed 1234 \
        --fake-strain-flow 10 \
        --sample-rate 16384 \
        --gps-start-time 1242442818 \
        --gps-end-time 1242442882 \
        --channel-name V1:SIMULATED_STRAIN \
        --injection-file injections.hdf \
        --output-strain-file V1-SIMULATED_STRAIN-1242442818-64.gwf

=======
Caveats
=======

An error will be raised if data is not available for the requested
time/detector combination.  See :ref:`example-valid-data` for an example of how
to determine when data is valid for a particular detector.

When applying filters to the data, care should be taken to use the
``--pad-data`` option appropriately, in order to remove the portion of the data
corrupted by the finite duration of the filter's impulse response.
``--pad-data`` defaults to several seconds, which is usually safe.

Depending on the value of the ``--pad-data`` option, ``pycbc_condition_strain``
might read more data than what specified by the ``--gps-start-time`` and
``--gps-end-time`` options. Errors will arise if ``--pad-data`` is set to a
value that causes ``pycbc_condition_strain`` to request data outside the range
of availability.

.. _GWOSC: https://www.gwosc.org/about/
