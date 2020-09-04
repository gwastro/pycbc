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
from scipy import interpolate

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
        beta_n = (n - 1) * 2.0 * np.pi / 3.0 + _lambda_
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

    def interpolate_waveform(self, hp, hc):
        """Returns the interpolated waveform
        Parameters
        ----------
        hp, hc : pycbc.types.TimeSeries
            plus and cross polarizations in order
        Returns
        -------
        numpy.ndarray
            interpolated waveform
        """
        fhp = scipy.interpolate.interpolate.interp1d(hp.sample_times, hp)
        fhc = scipy.interpolate.interpolate.interp1d(hc.sample_times, hc)
        return fhp(hp.sample_times), fhc(hc.sample_times)

    def eval_k(self, beta, lam):
        """Returns the unit propogation vector
        Parameters
        ----------
        beta, lam : numpy.Scalar
            position of the source in the sky in spherical coordinates
        Returns
        -------
        numpy.ndarray shape (1,3)
        """
        k = np.array([-cos(beta)*cos(lam), -cos(beta)*sin(lam), -sin(beta)])
        return k/(k[0]**2 + k[1]**2 + k[2]**2)**.5

    def eval_pol_tensor(self, beta, lam):
        """Evaluates the polarization basis tensors
        Parameters
        ----------
        beta, lam : numpy.Scalar
            position of the source in the sky in spherical coordinates
        Returns
        -------
        epij, ecij : two numpy.ndarray shape (3,3)
            polarization basis tensors in order
        """
        u = np.array([-sin(beta)*cos(lam), -sin(beta)*sin(lam), cos(beta)])
        v = np.array([sin(lam), -cos(lam), 0])
        return np.outer(u, u) - np.outer(v, v), np.outer(u, v) + np.outer(v, u)

    def eval_hij(self, waveform, beta, lam, psi, t):
        """Returns the values of hij
        Parameters
        ----------
        waveform : array of pycbc.types.TimeSeries
            plus and cross polarizations in order
        beta, lam : numpy.Scalar
            position of the source in the sky in spherical coordinates
        t : numpy.ndarray
            GPS time
        Returns
        -------
        numpy.ndarray shape (3,1)
            The values of hij at the required time
        """
        hp, hc = waveform
        epij, ecij = self.eval_pol_tensor(beta, lam)
        fhp, fhc = self.interpolate_waveform(hp, hc)
        h = []
        for i in range(len(t)):
            t0 = np.abs(hp.sample_times.data - t[i]).argmin()
            hij = np.zeros((3, 3))
            hij[0, 0:2] = (fhp[t0]*cos(2*psi) - fhc[t0]*sin(2*psi))*epij[0][0:2] + \
                          (fhp[t0]*sin(2*psi) - fhc[t0]*cos(2*psi))*ecij[0][0:2]
            hij[1, 0:2] = (fhp[t0]*cos(2*psi) - fhc[t0]*sin(2*psi))*epij[1][1:3] + \
                          (fhp[t0]*sin(2*psi) - fhc[t0]*cos(2*psi))*ecij[1][1:3]
            h.append(hij)
        return np.array(h)

    def GWresponse(self, waveform, beta, lam, psi, ref_time):
        """Returns the response to the GW signals
        Parameters
        ----------
        waveform : array of pycbc.types.TimeSeries
            plus and cross polarizations in order
        beta, lam : numpy.Scalar
            position of the source in the sky in spherical coordinates
        ref_time : numpy.Scalar
            GPS time
        Returns
        -------
        numpy.ndarray
            The values of response to GW signals.
        """
        guide_cen = np.transpose(np.mean(self.get_pos(ref_time), axis=1))
        guide_cen.shape = (3, 1)
        R = self.get_pos(ref_time) - guide_cen
        L_l = np.array([R[:, 1] - R[:, 0], R[:, 2] - R[:, 1], R[:, 0] - R[:, 2]])
        dist_L = np.array([np.sqrt(np.sum(L_l[0]**2)),
                           np.sqrt(np.sum(L_l[1]**2)),
                           np.sqrt(np.sum(L_l[2]**2))])
        dist_L.shape = (3, 1)
        n_hat = L_l/dist_L
        k = self.eval_k(beta, lam)
        kdotR = np.dot(k, R.transpose())
        kdotn = np.dot(k, n_hat.transpose())
        dist_L.shape = (3)

        phi_t1, phi_t2, phi_t_2 = [], [], []
        for i in range(3):
            for j in range(4):
                t1 = ref_time - kdotR - dist_L[0] - j*dist_L[0]
                t2 = ref_time - kdotR - j*dist_L[0]

                phi_t1.append(np.dot(np.dot(n_hat[i],
                              np.transpose(self.eval_hij(waveform, beta, lam, psi, t1))),
                              n_hat[i].transpose()))

                phi_t2.append(np.dot(np.dot(n_hat[i],
                              np.transpose(self.eval_hij(waveform, beta, lam, psi, t2))),
                              n_hat[i].transpose()))

                phi_t_2.append(np.dot(np.dot(n_hat[i],
                               np.transpose(self.eval_hij(waveform, beta, lam, psi, t2 + 2*kdotR))),
                               n_hat[i].transpose()))

        val = 2*(1 - kdotn)
        return np.array([phi_t1, phi_t2, phi_t_2]), val, kdotn


    def eval_XYZ(self, waveform, beta, lam, psi, ref_time):
        """Returns the values of X, Y, Z
        Parameters
        ----------
        waveform : array of pycbc.types.TimeSeries
            plus and cross polarizations in order
        beta, lam : numpy.Scalar
            position of the source in the sky in spherical coordinates
        ref_time : numpy.Scalar
            GPS time
        -------
        numpy.ndarray shape (3,1)
            The values of X, Y, Z at the ref_time.
        """
        data = self.GWresponse(waveform, beta, lam, psi, ref_time)
        Phi, den, kdotn = data[0], data[1], data[2]

        cycle = [([0, 1, 2] * 2)[x:x+3] for i in range(3) for x in [i % len([0, 1, 2])]]
        xyz = np.zeros(3)
        for i in range(3):
            a,b,c = cycle[i]
            xyz[i] = np.array([(Phi[0][3][a] - Phi[1][3][b])/(den + 4*kdotn)[c] +
                               (Phi[0][2][b] - Phi[2][2][a])/den[c] +
                               (Phi[0][1][b] - Phi[2][1][a])/den[c] +
                               (Phi[0][0][c] - Phi[1][0][a])/(den + 4*kdotn)[b] -
                               (Phi[0][3][a] - Phi[2][3][c])/den[b] -
                               (Phi[0][2][c] - Phi[1][2][a])/(den + 4*kdotn)[b] -
                               (Phi[0][1][a] - Phi[1][1][b])/(den + 4*kdotn)[c] -
                               (Phi[0][0][b] - Phi[2][0][a])/den[c]])
        return xyz 

    def eval_AET_from_XYZ(self, xyz):
        """Evaluates the values of E A T
        Parameters
        ----------
        xyz : numpy.ndarray shape (1,3)
            X, Y, Z in order
        Returns
        -------
        E A T : numpy.ndarray shape (1,3)
        """
        X, Y, Z = xyz[0], xyz[1], xyz[2]
        E = (X - 2*Y + Z) / 6**.5
        A = (Z - X) / 2**.5
        T = (X + Y + Z) / 3**.5
        return E, A, T
