###################################################
Reading Gravitational-wave Frames
###################################################

============
Introduction
============

All data generated and recorded by the current generation of ground-based laser-interferometer gravitational-wave detectors are recorded in gravitational-wave frame (GWF) files. These files typically contain data from a number of sources bundled into a single, time-stamped set, along with the metadata for each channel.

=====================
Querying a LDR server
=====================

The LIGO Data Replicator (LDR) is a tool for replicating data sets to the different data grids. If you have access to a LDR server you can read GWF files using ``pycbc.frame`` module as follows::

    >>> from pycbc import frame
    >>> tseries = frame.query_and_read_frame("G1_RDS_C01_L3", "G1:DER_DATA_H", 1049587200, 1049587200 + 60)

This returns a ``TimeSeries`` instance of the data. Note if you do not have access to frames through an LDR server then you will need to copy the frames to your run location.

Alternatively, if you just want the location of the frame files, you can do::

    >>> from pycbc import frame
    >>> frame_files = frame.frame_paths("G1_RDS_C01_L3", 1049587200, 1049587200 + 60)

This will return a ``list`` of the frame files' paths.

For public GWOSC strain, request 4 kHz frames with ``sample_rate=4096``::

    >>> strain = frame.query_and_read_frame('GWOSC_STRAIN', 'H1:GWOSC-4KHZ_R1_STRAIN', 1238166018, 1238166022, sample_rate=4096)

The default remains 4 kHz for S5, S6, and O1, and 16 kHz from O2 onward.
For O1 at 16 kHz, the strain reader selects the distinct
``GWOSC-16KHZ_R1_STRAIN`` channel automatically.
GWOSC run metadata uses PyCBC's existing download path and CI mirror; the GWF
files use that path with caching. Not every run has both rates; requesting an
unpublished rate raises ``ValueError``. In particular, the separate H1
background release around GW170608 is available only at 16 kHz through this
run-data interface.

=====================
Reading a frame file
=====================

The ``pycbc.frame`` module provides methods for reading these files into ``TimeSeries`` objects as follows::

    >>> from pycbc import frame
    >>> data = frame.read_frame('G-G1_RDS_C01_L3-1049587200-60.gwf', 'G1:DER_DATA_H', 1049587200, 1049587200 + 60)

Here the first argument is the path to the frame file of interest, while the second lists the `data channel` of interest whose data exist within the file.

=====================
Writing a frame file
=====================

The ``pycbc.frame`` modules provides a method for writing ``TimeSeries`` instances to GWF as follows::

    >>> from pycbc import frame
    >>> from pycbc import types
    >>> data = types.TimeSeries(numpy.ones(16384 * 16), delta_t=1.0 / 16384)
    >>> frame.write_frame("./test.gwf", "H1:TEST_UNITY", data)

Here the first argument is the path to the output frame file, the second is the name of the channel, and the last is the ``TimeSeries`` instance to be written to the frame.

====================
Method documentation
====================

.. automodule:: pycbc.frame
    :noindex:
    :members:
