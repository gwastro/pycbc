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
import lalsimulation
import numpy as np
import lal
from numpy import cos, sin
from pycbc.types import TimeSeries


class Detector(object):
    """A gravitaional wave detector
    """
    def __init__(self, detector_name):
        self.name = str(detector_name)
        self.frDetector =  lalsimulation.DetectorPrefixToLALDetector(self.name)
        self.response = self.frDetector.response
        self.location = self.frDetector.location
        self.latitude = self.frDetector.frDetector.vertexLatitudeRadians
        self.longitude = self.frDetector.frDetector.vertexLongitudeRadians

    def light_travel_time_to_detector(self, det):
        """ Return the light travel time from this detector

        Parameters
        ----------
        detector: Detector
            The other detector to determine the light travel time to.

        Returns
        -------
        time: float
            The light travel time in seconds
        """
        return lal.LightTravelTime(self.frDetector, det.frDetector) * 1e-9

    def antenna_pattern(self, right_ascension, declination, polarization, t_gps):
        """Return the detector response.
        """
        gmst = lal.GreenwichMeanSiderealTime(t_gps)
        return tuple(lal.ComputeDetAMResponse(self.response,
                     right_ascension, declination, polarization, gmst))

    def time_delay_from_earth_center(self, right_ascension, declination, t_gps):
        """Return the time delay from the earth center
        """
        return lal.TimeDelayFromEarthCenter(self.location,
                      float(right_ascension), float(declination), float(t_gps))

    def time_delay_from_detector(self, other_detector, right_ascension,
                                 declination, t_gps):
        """Return the time delay from the given to detector for a signal with
        the given sky location; i.e. return `t1 - t2` where `t1` is the
        arrival time in this detector and `t2` is the arrival time in the
        other detector. Note that this would return the same value as
        `time_delay_from_earth_center` if `other_detector` was geocentric.

        Parameters
        ----------
        other_detector : detector.Detector
            A detector instance.
        right_ascension : float
            The right ascension (in rad) of the signal.
        declination : float
            The declination (in rad) of the signal.
        t_gps : float
            The GPS time (in s) of the signal.

        Returns
        -------
        float
            The arrival time difference between the detectors.
        """
        return lal.ArrivalTimeDiff(self.location, other_detector.location,
                                   float(right_ascension), float(declination),
                                   float(t_gps))

    def project_wave(self, hp, hc, longitude, latitude, polarization):
        """Return the strain of a wave with given amplitudes and angles as
        measured by the detector.
        """
        h_lal = lalsimulation.SimDetectorStrainREAL8TimeSeries(
                hp.astype(np.float64).lal(), hc.astype(np.float64).lal(),
                longitude, latitude, polarization, self.frDetector)
        return TimeSeries(
                h_lal.data.data, delta_t=h_lal.deltaT, epoch=h_lal.epoch,
                dtype=np.float64, copy=False)

    def optimal_orientation(self, t_gps):
        """Return the optimal orientation in right ascension and declination
           for a given GPS time.
        """
        ra = self.longitude + (lal.GreenwichMeanSiderealTime(t_gps) % (2*np.pi))
        dec = self.latitude
        return ra, dec

def overhead_antenna_pattern(right_ascension, declination, polarization):
    """Return the detector response where (0, 0) indicates an overhead source. 
    This functions uses coordinates such that the detector can be thought to
    be on the north pole.

    Parameters
    ----------
    right_ascention: float
    declination: float
    polarization: float

    Returns
    -------
    f_plus: float
    f_cros: float   
    """
    # convert from declination coordinate to polar (angle dropped from north axis)
    theta = np.pi / 2.0 - declination

    f_plus  = - (1.0/2.0) * (1.0 + cos(theta)*cos(theta)) * \
                cos (2.0 * right_ascension) * cos (2.0 * polarization) - \
                cos(theta) * sin(2.0*right_ascension) * sin (2.0 * polarization)

    f_cross =   (1.0/2.0) * (1.0 + cos(theta)*cos(theta)) * \
                cos (2.0 * right_ascension) * sin (2.0* polarization) - \
                cos(theta) * sin(2.0*right_ascension) * cos (2.0 * polarization)

    return f_plus, f_cross

def effective_distance(distance, inclination, f_plus, f_cross):
    return distance / np.sqrt( ( 1 + np.cos( inclination )**2 )**2 / 4 * f_plus**2 + np.cos( inclination )**2 * f_cross**2 )

