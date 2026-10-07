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

==============================================================
Spacecraft orbit interface
==============================================================

Native spacecraft orbit providers inherit from
:class:`pycbc.detector.orbits.BaseOrbit`. The interface separates constellation
kinematics from coordinate transformations and from the detector/TDI response.
It does not prescribe a particular mission, orbit model or file format.

Providers expose ``t_interp`` (sample times, or ``None`` for an analytic orbit)
and ``num_sc``, and implement ``compute_position(t, sc=None)``,
``compute_velocity(t, sc=None)`` and ``compute_acceleration(t, sc=None)``.
Spacecraft labels start at 1. Scalar inputs keep both time and spacecraft
axes: requesting one spacecraft at one time returns shape ``(1, 1, 3)``.
Vector requests preserve input order, including repeated times or labels.
``sc=None`` returns all spacecraft.

All three methods use fixed J2000 ecliptic axes about the solar-system
barycentre, in metres, metres/second and metres/second squared, respectively.
A provider must document its coordinate-time epoch, time scale, supported
interval and extrapolation policy. The base class performs no epoch or frame
conversion. Evaluation times and ``t_interp`` use the same time convention.

The base constructor validates interpolation times and spacecraft count.
Subclasses can reuse ``_prepare_times`` and ``_sc_indices`` for evaluation
inputs; they remain responsible for returning the documented shapes and units.
In particular, interpolation times are increasing, whereas evaluation times
may be unordered or repeated (as needed for retarded-time calculations).
Interpolation, file readers and derivative algorithms belong to concrete
providers, not to the ABC. Existing third-party providers can still be used
through their compatible methods without inheriting from this base class.

.. autoclass:: pycbc.detector.orbits.BaseOrbit
   :members:
