########################################################
Manipulating frame files with ``pycbc_condition_strain``
########################################################

``pycbc_condition_strain`` is a command-line utility which, despite the obscure
name, can be used for a variety of generic and common operations related to
frame files. The basic idea of the utility is to obtain a segment of
gravitational-wave strain data, apply a number of operations that a
matched-filter search would apply before matched-filtering, and then write the
transformed data back to frame files instead of running the matched filtering.

The following sections show some concrete examples.

================================================
Splitting existing data into shorter frame files
================================================

Read existing gravitational-wave data (either from pre-existing frame files
or from a data server such as GWOSC) and write it back split into frame files
of a given duration.

==========================
Conditioning existing data
==========================

Read existing gravitational-wave data, apply various conditioning procedures
to it (e.g. applying a high-pass or low-pass filter, resampling, gating, or
reducing the precision) and write the result to one or more frame files.

``pycbc_condition_strain`` was initially created with this specific task in
mind, hence the name.

==============================================
Injecting simulated signals into existing data
==============================================

Read existing gravitational-wave data, add one or more simulated signals
into it, and write the result to one or more frame files.

==================================
Simulating gravitational-wave data
==================================

Simulate detector noise with a given power spectral density, optionally
add one or more simulated signals to it, and write the result to one or more
frame files.
