# -*- coding: UTF-8 -*-

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
"""This module provides utilities for calculating detector responses and timing
between observatories.
"""
import lalsimulation
import numpy as np
import lal
from pycbc.types import TimeSeries
from astropy.time import Time
from astropy import constants
from astropy.units.si import sday
from numpy import cos, sin

# Response functions are modelled after those in lalsuite and as also
# presented in https://arxiv.org/pdf/gr-qc/0008066.pdf

def gmst_accurate(gps_time):
    gmst = Time(gps_time, format='gps', scale='utc',
                location=(0, 0)).sidereal_time('mean').rad
    return gmst


def get_available_detectors():
    """Return list of detectors known in the currently sourced lalsuite.
    This function will query lalsuite about which detectors are known to
    lalsuite. Detectors are identified by a two character string e.g. 'K1',
    but also by a longer, and clearer name, e.g. KAGRA. This function returns
    both. As LAL doesn't really expose this functionality we have to make some
    assumptions about how this information is stored in LAL. Therefore while
    we hope this function will work correctly, it's possible it will need
    updating in the future. Better if lal would expose this information
    properly.
    """
    ld = lal.__dict__
    known_lal_names = [j for j in ld.keys() if "DETECTOR_PREFIX" in j]
    known_prefixes = [ld[k] for k in known_lal_names]
    known_names = [ld[k.replace('PREFIX', 'NAME')] for k in known_lal_names]
    return list(zip(known_prefixes, known_names))


class Detector(object):
    """A gravitational wave detector
    """
    def __init__(self, detector_name, reference_time=1126259462.0):
        """ Create class representing a gravitational-wave detector
        Parameters
        ----------
        detector_name: str
            The two-character detector string, i.e. H1, L1, V1, K1, I1
        reference_time: float
            Default is time of GW150914. In this case, the earth's rotation
        will be estimated from a reference time. If 'None', we will
        calculate the time for each gps time requested explicitly
        using a slower but higher precision method.
        """
        self.name = str(detector_name)
        self.frDetector = lalsimulation.DetectorPrefixToLALDetector(self.name)
        self.response = self.frDetector.response
        self.location = self.frDetector.location
        self.latitude = self.frDetector.frDetector.vertexLatitudeRadians
        self.longitude = self.frDetector.frDetector.vertexLongitudeRadians

        self.reference_time = reference_time
        self.sday = None
        self.gmst_reference = None

    def set_gmst_reference(self):
        if self.reference_time is not None:
            self.sday = float(sday.si.scale)
            self.gmst_reference = gmst_accurate(self.reference_time)
        else:
            raise RuntimeError("Can't get accurate sidereal time without GPS "
                               "reference time!")

    def gmst_estimate(self, gps_time):
        if self.reference_time is None:
            return gmst_accurate(gps_time)

        if self.gmst_reference is None:
            self.set_gmst_reference()
        dphase = (gps_time - self.reference_time) / self.sday * (2.0 * np.pi)
        gmst = (self.gmst_reference + dphase) % (2.0 * np.pi)
        return gmst

    def light_travel_time_to_detector(self, det):
        """ Return the light travel time from this detector
        Parameters
        ----------
        det: Detector
            The other detector to determine the light travel time to.
        Returns
        -------
        time: float
            The light travel time in seconds
        """
        d = self.location - det.location
        return float(d.dot(d)**0.5 / constants.c.value)

    def antenna_pattern(self, right_ascension, declination, polarization, t_gps):
        """Return the detector response.
        Parameters
        ----------
        right_ascension: float or numpy.ndarray
            The right ascension of the source
        declination: float or numpy.ndarray
            The declination of the source
        polarization: float or numpy.ndarray
            The polarization angle of the source
        Returns
        -------
        fplus: float or numpy.ndarray
            The plus polarization factor for this sky location / orientation
        fcross: float or numpy.ndarray
            The cross polarization factor for this sky location / orientation
        """
        if isinstance(t_gps, lal.LIGOTimeGPS):
            t_gps = float(t_gps)
        gha = self.gmst_estimate(t_gps) - right_ascension

        cosgha = cos(gha)
        singha = sin(gha)
        cosdec = cos(declination)
        sindec = sin(declination)
        cospsi = cos(polarization)
        sinpsi = sin(polarization)

        x0 = -cospsi * singha - sinpsi * cosgha * sindec
        x1 = -cospsi * cosgha + sinpsi * singha * sindec
        x2 =  sinpsi * cosdec
        x = np.array([x0, x1, x2])

        dx = self.response.dot(x)

        y0 =  sinpsi * singha - cospsi * cosgha * sindec
        y1 =  sinpsi * cosgha + cospsi * singha * sindec
        y2 =  cospsi * cosdec
        y = np.array([y0, y1, y2])
        dy = self.response.dot(y)

        if hasattr(dx, 'shape'):
            fplus = (x * dx - y * dy).sum(axis=0)
            fcross = (x * dy + y * dx).sum(axis=0)
        else:
            fplus = (x * dx - y * dy).sum()
            fcross = (x * dy + y * dx).sum()

        return fplus, fcross

    def time_delay_from_earth_center(self, right_ascension, declination, t_gps):
        """Return the time delay from the earth center
        """
        return self.time_delay_from_location(np.array([0, 0, 0]),
                                             right_ascension,
                                             declination,
                                             t_gps)

    def time_delay_from_location(self, other_location, right_ascension,
                                 declination, t_gps):
        """Return the time delay from the given location to detector for
        a signal with the given sky location
        In other words return `t1 - t2` where `t1` is the
        arrival time in this detector and `t2` is the arrival time in the
        other location.
        Parameters
        ----------
        other_location : numpy.ndarray of coordinates
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
        ra_angle = self.gmst_estimate(t_gps) - right_ascension
        cosd = cos(declination)

        e0 = cosd * cos(ra_angle)
        e1 = cosd * -sin(ra_angle)
        e2 = sin(declination)

        ehat = np.array([e0, e1, e2])
        dx = other_location - self.location
        return dx.dot(ehat) / constants.c.value

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
        return self.time_delay_from_location(other_detector.location,
                                             right_ascension,
                                             declination,
                                             t_gps)

    def project_wave(self, hp, hc, longitude, latitude, polarization):
        """Return the strain of a waveform as measured by the detector.
        Apply the time shift for the given detector relative to the assumed
        geocentric frame and apply the antenna patterns to the plus and cross
        polarizations.
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
        Parameters
        ----------
        t_gps: float
            Time in gps seconds
        Returns
        -------
        ra: float
            Right ascension that is optimally oriented for the detector
        dec: float
            Declination that is optimally oriented for the detector
        """
        ra = self.longitude + (self.gmst_estimate(t_gps) % (2.0*np.pi))
        dec = self.latitude
        return ra, dec

def overhead_antenna_pattern(right_ascension, declination, polarization):
    """Return the antenna pattern factors F+ and Fx as a function of sky
    location and polarization angle for a hypothetical interferometer located
    at the north pole. Angles are in radians. Declinations of ±π/2 correspond
    to the normal to the detector plane (i.e. overhead and underneath) while
    the point with zero right ascension and declination is the direction
    of one of the interferometer arms.
    Parameters
    ----------
    right_ascension: float
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

""" LISA class  """


class LISA(object):
    def __init__(self,kappa,_lambda_,reference_time=1126259462.0):
        self.reference_time=reference_time
        self.kappa=kappa
        self._lambda_=_lambda_

    def get_pos(self, t_gps):
        if t_gps is None:
            t_gps = Time(val = self.reference_time, format = 'gps',
                  scale = 'utc').to_datetime(timezone = None)
        elif isinstance(t_gps, np.ScalarType):
            t_gps = Time(val = t_gps, format = 'gps',
                  scale = 'utc').to_datetime(timezone = None)

        t_gps = np.sum(np.array([t_gps.year - 2034, t_gps.month/12, t_gps.day/(12*365),
                             t_gps.hour/(12*365*24), 
                             t_gps.minute/(12*365*24*60),
                             t_gps.second/(12*365*24*60*60),
                             t_gps.microsecond/(12*365*24*60*60*1e-6)]), axis=0)
        
        n = np.array(range(1, 4))
        kappa, _lambda_ = 0, 0
        alpha = 2. * np.pi * t_gps/1 + kappa
        beta_n = (n - 1) + (2. * np.pi/3) + _lambda_
        a, L = 1., .1   # units are in AU
        e = L/(2. * a * np.sqrt(3))

        x = a*cos(alpha) + a*e*(sin(alpha)*cos(alpha)*sin(beta_n) - (1 + sin(alpha)**2)*cos(beta_n))
        y = a*sin(alpha) + a*e*(sin(alpha)*cos(alpha)*sin(beta_n) - (1 + cos(alpha)**2)*sin(beta_n))
        z = -np.sqrt(3)*a*e*cos(alpha - beta_n)
        self.location = np.array([x,y,z])

        return self.location

    def plot_orbit(self):
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt
        dec_center = []
        for i in range(2000):
          dec_center.append(self.get_pos(self.reference_time + i).mean(axis = 1))
        """ values don't change much with GPS time"""
        t = Time(val = self.reference_time, format = 'gps', scale ='utc')
        sun = coordinates.get_sun(t).transform_to('icrs')
        earth = coordinates.get_body('earth', t, location = None).transform_to('icrs')
        sun.representation_type, earth.representation_type ='cartesian', 'cartesian'

        fig = plt.figure()
        ax = plt.axes(projection = "3d")
        ax.scatter(np.float32(earth.x), np.float32(earth.y), np.float32(earth.z), marker = 'o')
        ax.scatter(np.float32(sun.x), np.float32(sun.y), np.float32(sun.z), marker = '+')
        ax.scatter(dec_center[0], dec_center[1] ,dec_center[2], marker = '*')
        ax.set_xlabel('X axis (AU)')
        ax.set_ylabel('Y axis (AU)')
        ax.set_zlabel('Z axis (AU)')

    def time_delay_from_location(self, other_location, right_ascension,
                                 declination, t_gps):
        dec_loc = self.get_pos(t_gps)
        """signal = coordinates.SkyCoord(ra = right_ascension, dec = declination,
                                          unit = u.rad, frame = 'gcrs').transform_to('icrs')"""

        dx = np.array([other_location[0] - self.location[0],
                       other_location[1] - self.location[1],
                       other_location[2] - self.location[2]])

        """ra_angle = self.gmst_estimate(t_gps) - right_ascension"""
        cosd = cos(declination)
        e0 = cosd * cos(right_ascension)
        e1 = cosd * -sin(right_ascension)
        e2 = sin(declination)
        ehat = np.array([e0, e1, e2])
        return dx.dot(ehat) / constants.c.value

    def time_delay_from_detector(self, other_detector, right_ascension,
                                 declination, t_gps):
        return self.time_delay_from_location(other_detector.location,
                                             right_ascension,
                                             declination,
                                             t_gps)
    def time_delay_from_earth_center(self, right_ascension, declination, t_gps):
        if t_gps is None:
          t_gps = Time(val = self.reference_time, format = 'gps', scale ='utc')
        else:
          t_gps = Time(val = t_gps, format = 'gps', scale ='utc')
        earth = coordinates.get_body('earth', t, location = None).transform_to('icrs')
        return self.time_delay_from_location(
            np.array([np.float32(earth.x), np.float32(earth.y), np.float32(earth.z)]),
            right_ascension, declination, t_gps)        
