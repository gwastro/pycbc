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

=====================
Reading a frame file
=====================

The ``pycbc.frame`` module provides methods for reading these files into ``TimeSeries`` objects as follows::

    >>> from pycbc import frame
    >>> data = frame.read_frame('G-G1_RDS_C01_L3-1049587200-60.gwf', 'G1:DER_DATA_H')

Here the first argument is the path to the framefile of interest, while the second lists the `data channel` of interest whose data exist within the file.

====================
Method documentation
====================

.. automodule:: pycbc.frame
    :noindex:
    :members:
