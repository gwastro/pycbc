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
"""This module provides utilities for calculating detector responses.
"""
import lal
import numpy
from math import cos, sin

# get the cached detector data
_detectors = {}

for _detector in lal.lalCachedDetectors:
    _detectors[_detector.frDetector.name]=_detector

def detectors():
    """Return a list of detectors. 
    """
    return _detectors.keys()

def get_detector(detector_name):
    """
    """
    return _detectors[detector_name]

def detector_antenna_pattern(detector_name, right_ascension, declination, 
                                 polarization, gmst):
    """Return the detector response.
    """
    return tuple(lal.ComputeDetAMResponse(_detectors[detector_name].response, 
                 right_ascension, declination, polarization, gmst))

def overhead_antenna_pattern(right_ascension, declination, polarization):
    """Return the detector response. 
    """
    f_plus = - (1.0/2.0) * (1.0 + cos(declination)*cos(declination)) * cos (2.0 * right_ascension) * cos (2.0 * polarization) - cos(declination) * sin(2.0*right_ascension) * sin (2.0 * polarization)
    f_cross=  (1.0/2.0) * (1.0 + cos(declination)*cos(declination)) * cos (2.0 * right_ascension) * sin (2.0* polarization) - cos (declination) * sin(2.0*right_ascension) * cos (2.0 * polarization)
    return f_plus, f_cross








