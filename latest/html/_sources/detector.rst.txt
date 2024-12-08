###################################################
Gravitational-wave Detectors
###################################################

The pycbc.detector module provides the :py:mod:`pycbc.detector.Detector` class
to access information about gravitational wave detectors and key information
about how their orientation and position affects their view of a source

=====================================
Detector Locations
=====================================

.. literalinclude:: ../examples/detector/loc.py
.. command-output:: python ../examples/detector/loc.py

=====================================
Light travel time between detectors
=====================================

.. literalinclude:: ../examples/detector/travel.py
.. command-output:: python ../examples/detector/travel.py

======================================================
Time source gravitational-wave passes through detector
======================================================

.. literalinclude:: ../examples/detector/delay.py
.. command-output:: python ../examples/detector/delay.py

================================================================
Antenna Patterns and Projecting a Signal into the Detector Frame
================================================================

.. literalinclude:: ../examples/detector/ant.py
.. command-output:: python ../examples/detector/ant.py
   
==============================================================
Adding a custom detector / overriding existing ones
==============================================================
PyCBC supports observatories with arbitrary locations. For the study
of possible new observatories you can add them explicitly within a script
or by means of a config file to make the detectors visible to all codes
that use the PyCBC detector interfaces.

An example of adding a detector directly within a script.

.. plot:: ../examples/detector/custom.py
   :include-source:


The following demonstrates a config file which similarly can provide
custom observatory information. The options are the same as for the
direct function calls. To tell PyCBC the location of the config file, 
set the PYCBC_DETECTOR_CONFIG variable to the location of the file e.g.
PYCBC_DETECTOR_CONFIG=/some/path/to/detectors.ini. The following would
provide new detectors 'f1' and 'f2'.

.. literalinclude:: ../examples/detector/custom.ini
