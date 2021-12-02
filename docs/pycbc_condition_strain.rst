########################################################
Manipulating frame files with ``pycbc_condition_strain``
########################################################

``pycbc_condition_strain`` is a command-line tool which can be used for a
variety of operations related to gravitational-wave strain data files.
The basic idea of the tool is to combine the following operations:

1. Obtain a segment of gravitational-wave strain data, either by simulating
detector noise, or reading an existing data file, or downloading the data from
a server like `GWOSC <https://www.gw-openscience.org/about/>`_.
2. Optionally, simulate one or more signals and add them to the data.
3. Optionally, apply a number of operations that a matched-filter search would
apply before matched-filtering, such as low- or high-passing, resampling or
glitch gating. These operations are sometimes referred to as "data
conditioning", hence the name of the tool.
4. Write the resulting data back to one or more files. Different output file
formats are supported.

The following sections show concrete examples of common use cases.  For a
detailed list of all options, and option-specific help, see
``pycbc_condition_strain --help``.

================================================
Splitting existing data into shorter frame files
================================================

Obtain a 64 s segment of gravitational-wave data from `GWOSC
<https://www.gw-openscience.org/about/>`_ and write it to multiple frame files
of 16 s duration:::

    pycbc_condition_strain \
        --gps-start-time 1242442818 \
        --gps-end-time 1242442882 \
        --frame-type LOSC \
        --channel-name H1:GWOSC-16KHZ_R1_STRAIN \
        --frame-duration 16 \
        --output-strain-file 'H1-STRAIN-{start}-{duration}.gwf'

Note that, when requesting GWOSC data, an error will be raised if the requested
time/detector combination does not have usable data.

==========================
Conditioning existing data
==========================

Obtain 64 s of gravitational-wave data from `GWOSC
<https://www.gw-openscience.org/about/>`_, apply a high-pass filter with corner
frequency of 15 Hz, downsample the data to 2048 Hz, and write it to a single
frame file:::

    pycbc_condition_strain \
        --gps-start-time 1242442818 \
        --gps-end-time 1242442882 \
        --frame-type LOSC \
        --channel-name H1:GWOSC-16KHZ_R1_STRAIN \
        --strain-high-pass 15 \
        --sample-rate 2048 \
        --output-strain-file 'H1-STRAIN_CONDITIONED-1242442818-64.gwf'

==============================================
Injecting simulated signals into existing data
==============================================

Read 64 s of existing gravitational-wave data, add one or more simulated
signals into it (specified by the ``injections.hdf`` file), and write the
result to frame files of 16 s duration:::

    pycbc_condition_strain \
        --gps-start-time 1242442818 \
        --gps-end-time 1242442882 \
        --frame-type LOSC \
        --channel-name H1:GWOSC-16KHZ_R1_STRAIN \
        --injection-file injections.hdf \
        --frame-duration 16 \
        --output-strain-file 'H1-STRAIN_WITH_INJECTIONS-{start}-{duration}.gwf'

==================================
Simulating gravitational-wave data
==================================

Simulate 64 s of detector noise with a given power spectral density, add one or more
simulated signals to it (specified by the ``injections.hdf`` file), and write
the result to a single frame file:::

    pycbc_condition_strain \
        --fake-strain aLIGOAdVO4T1800545 \
        --fake-strain-seed 1234 \
        --fake-strain-flow 10 \
        --fake-strain-sample-rate 16384 \
        --gps-start-time 1242442818 \
        --gps-end-time 1242442882 \
        --channel-name V1:SIMULATED_STRAIN \
        --injection-file injections.hdf \
        --output-strain-file 'V1-SIMULATED_STRAIN-1242442818-64.gwf'
