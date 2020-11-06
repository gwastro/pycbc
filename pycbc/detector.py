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
from astropy import constants, coordinates, units
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

    def antenna_pattern(self, right_ascension, declination, polarization, t_gps, polarization_type='tensor'):
        """Return the detector response.

        Parameters
        ----------
        right_ascension: float or numpy.ndarray
            The right ascension of the source
        declination: float or numpy.ndarray
            The declination of the source
        polarization: float or numpy.ndarray
            The polarization angle of the source
        polarization_type: string flag: Tensor, Vector or Scalar
            The gravitational wave polarizations. Default: 'Tensor'

        Returns
        -------
        fplus(default) or fx or fb : float or numpy.ndarray
            The plus or vector-x or breathing polarization factor for this sky location / orientation
        fcross(default) or fy or fl : float or numpy.ndarray
            The cross or vector-y or longitudnal polarization factor for this sky location / orientation
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

        z0 = -cosdec * cosgha
        z1 = cosdec * singha
        z2 = -sindec
        z = np.array([z0, z1, z2])

        dz = self.response.dot(z)

        if polarization_type == 'tensor':
            if hasattr(dx, 'shape'):
                fplus = (x * dx - y * dy).sum(axis=0)
                fcross = (x * dy + y * dx).sum(axis=0)
            else:
                fplus = (x * dx - y * dy).sum()
                fcross = (x * dy + y * dx).sum()               
            return fplus, fcross

        elif polarization_type == 'vector':
            if hasattr(dx, 'shape'):
                fx = (z * dx + x * dz).sum(axis=0)
                fy = (z * dy + y * dz).sum(axis=0)
            else:
                fx = (z * dx + x * dz).sum()
                fy = (z * dy + y * dz).sum()
            return fx, fy

        elif polarization_type == 'scalar':
            if hasattr(dx, 'shape'):
                fb = (x * dx + y * dy).sum(axis=0)
                fl = (z * dz).sum(axis=0)
            else:
                fb = (x * dx + y * dy).sum()
                fl = (z * dz).sum()
            return fb, fl

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

    def get_icrs_pos(self):
        """ Transforms GCRS frame to ICRS frame
        Returns
        ----------
        loc: numpy.ndarray shape (3,1) units: AU
             ICRS coordinates in cartesian system
        """
        loc = self.location
        loc = coordinates.SkyCoord(x=loc[0], y=loc[1], z=loc[2], unit=units.m,
                frame='gcrs', representation_type='cartesian').transform_to('icrs')
        loc.representation_type = 'cartesian'
        conv = np.float32(((loc.x.unit/units.AU).decompose()).to_string())
        loc = np.array([np.float32(loc.x), np.float32(loc.y),
                        np.float32(loc.z)])*conv
        return loc

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


"""     LISA class      """


class LISA(object):
    """For LISA detector
    """
    def __init__(self):
        None

    def get_pos(self, ref_time):
        """Return the position of LISA detector for a given reference time
        Parameters
        ----------
        ref_time : numpy.ScalarType
        Returns
        -------
        location : numpy.ndarray of shape (3,3)
               Returns the position of all 3 sattelites with each row
               correspoding to a single axis.
        """
        ref_time = Time(val=ref_time, format='gps', scale='utc').jyear
        n = np.array(range(1, 4))
        kappa, _lambda_ = 0, 0
        alpha = 2. * np.pi * ref_time/1 + kappa
        beta_n = (n - 1) * 2.0 * pi / 3.0 + _lambda_
        a, L = 1., 0.03342293561
        e = L/(2. * a * np.sqrt(3))

        x = a*cos(alpha) + a*e*(sin(alpha)*cos(alpha)*sin(beta_n) - (1 + sin(alpha)**2)*cos(beta_n))
        y = a*sin(alpha) + a*e*(sin(alpha)*cos(alpha)*cos(beta_n) - (1 + cos(alpha)**2)*sin(beta_n))
        z = -np.sqrt(3)*a*e*cos(alpha - beta_n)
        self.location = np.array([x, y, z])

        return self.location

    def get_gcrs_pos(self, location):
        """ Transforms ICRS frame to GCRS frame
        Parameters
        ----------
        loc : numpy.ndarray shape (3,1) units: AU
              Cartesian Coordinates of the location
              in ICRS frame
        Returns
        ----------
        loc : numpy.ndarray shape (3,1) units: meters
              GCRS coordinates in cartesian system
        """
        loc = location
        loc = coordinates.SkyCoord(x=loc[0], y=loc[1], z=loc[2], unit=units.AU,
                frame='icrs', representation_type='cartesian').transform_to('gcrs')
        loc.representation_type = 'cartesian'
        conv = np.float32(((loc.x.unit/units.m).decompose()).to_string())
        loc = np.array([np.float32(loc.x), np.float32(loc.y),
                        np.float32(loc.z)])*conv
        return loc

    def time_delay_from_location(self, other_location, right_ascension,
                                 declination, t_gps):
        """Return the time delay from the LISA detector to detector for
        a signal with the given sky location. In other words return
        `t1 - t2` where `t1` is the arrival time in this detector and
        `t2` is the arrival time in the other location. Units(AU)
        Parameters
        ----------
        other_location : numpy.ndarray of coordinates in ICRS frame
            A detector instance.
        right_ascension : float
            The right ascension (in rad) of the signal.
        declination : float
            The declination (in rad) of the signal.
        t_gps : float
            The GPS time (in s) of the signal.
        Returns
        -------
        numpy.ndarray
            The arrival time difference between the detectors.
        """
        dx = np.array([self.location[0] - other_location[0],
                       self.location[1] - other_location[1],
                       self.location[2] - other_location[2]])
        cosd = cos(declination)
        e0 = cosd * cos(right_ascension)
        e1 = cosd * -sin(right_ascension)
        e2 = sin(declination)
        ehat = np.array([e0, e1, e2])
        return dx.dot(ehat) / constants.c.value

    def time_delay_from_detector(self, det, right_ascension,
                                 declination, t_gps):
        """Return the time delay from the LISA detector for a signal with
        the given sky location in ICRS frame; i.e. return `t1 - t2` where
        `t1` is the arrival time in this detector and `t2` is the arrival
        time in the other detector.
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
        numpy.ndarray
            The arrival time difference between the detectors.
        """
        loc = Detector(det, t_gps).get_icrs_pos()
        return self.time_delay_from_location(loc, right_ascension,
                                             declination, t_gps)

    def time_delay_from_earth_center(self, right_ascension, declination, t_gps):
        """Return the time delay from the earth center in ICRS frame
        """
        t_gps = Time(val=t_gps, format='gps', scale='utc')
        earth = coordinates.get_body('earth', t_gps,
                                     location=None).transform_to('icrs')
        earth.representation_type = 'cartesian'
        return self.time_delay_from_location(
            np.array([np.float32(earth.x), np.float32(earth.y),
                      np.float32(earth.z)]), right_ascension,
            declination, t_gps)
