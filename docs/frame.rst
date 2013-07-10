###################################################
Reading Gravitational-wave Frames (``pycbc.frame``)
###################################################

============
Introduction
============

All data generated and recorded by the current generation of ground-based laser-interferometer gravitational-wave detectors are recorded in gravitational-wave frame (GWF) files. These files typically contain data from a number of sources bundled into a single, time-stamped set, along with the metadata for each channel.

The ``pycbc.frame`` module provides methods for reading these files into ``TimeSeries`` objects as follows::

    >>> from pycbc import frame
    >>> data = frame.read_frame('G-G1_RDS_C01_L3-1049587200-60.gwf', 'G1:DER_DATA_H')

Here the first argument is the path to the framefile of interest, while the second lists the `data channel` of interest whose data exist within the file.

====================
Method documentation
====================

.. autofunction:: pycbc.frame.read_frame
