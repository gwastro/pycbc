# Copyright (C) 2012  Alex Nitz
#
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""
"""
import lal
import numpy

# get the cached detector data
_detectors = {}

for _detector in lal.lalCachedDetectors:
    _detectors[_detector.frDetector.name]=_detector

def detectors():
    return _detectors.keys()

def antenna_pattern(detector_name, right_ascension, declination, polarization, gmst):
    """
    """
    return tuple(lal.ComputeDetAMResponse(_detectors[detector_name].response, right_ascension, 
                                    declination, polarization, gmst))

def fiducial_antenna_pattern(right_ascension, declination, polarization):
    """
    """
    no_detector = numpy.array([[1, 0, 0], [0, 1, 0], [0, 0, 0]], 
                                             dtype=numpy.float32)
    return tuple(lal.ComputeDetAMResponse(no_detector, right_ascension, 
                                    declination, polarization, 0))
    








