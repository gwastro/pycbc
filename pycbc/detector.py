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
import os
import logging
import numpy as np
from numpy import cos, sin, pi

import lal
from astropy.time import Time
from astropy import constants, coordinates, units
from astropy.coordinates.matrix_utilities import rotation_matrix
from astropy.units.si import sday, meter

import pycbc.libutils
from pycbc.types import TimeSeries
from pycbc.types.config import InterpolatingConfigParser

logger = logging.getLogger('pycbc.detector')

# Response functions are modelled after those in lalsuite and as also
# presented in https://arxiv.org/pdf/gr-qc/0008066.pdf

def gmst_accurate(gps_time):
    gmst = Time(gps_time, format='gps', scale='utc',
                location=(0, 0)).sidereal_time('mean').rad
    return gmst

def get_available_detectors():
    """ List the available detectors """
    dets = list(_ground_detectors.keys())
    return dets

def get_available_lal_detectors():
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

_ground_detectors = {}

def add_detector_on_earth(name, longitude, latitude,
                          yangle=0, xangle=None, height=0,
                          xlength=4000, ylength=4000,
                          xaltitude=0, yaltitude=0):
    """ Add a new detector on the earth

    Parameters
    ----------

    name: str
        two-letter name to identify the detector
    longitude: float
        Longitude in radians using geodetic coordinates of the detector
    latitude: float
        Latitude in radians using geodetic coordinates of the detector
    yangle: float
        Azimuthal angle of the y-arm (angle drawn from pointing north)
    xangle: float
        Azimuthal angle of the x-arm (angle drawn from point north). If not set
        we assume a right angle detector following the right-hand rule.
    xaltitude: float
        The altitude angle of the x-arm measured from the local horizon.
    yaltitude: float
        The altitude angle of the y-arm measured from the local horizon.
    height: float
        The height in meters of the detector above the standard
        reference ellipsoidal earth
    """
    if xangle is None:
        # assume right angle detector if no separate xarm direction given
        xangle = yangle + np.pi / 2.0

    # baseline response of a single arm pointed in the -X direction
    resp = np.array([[-1, 0, 0], [0, 0, 0], [0, 0, 0]])
    rm2 = rotation_matrix(-longitude * units.rad, 'z')
    rm1 = rotation_matrix(-1.0 * (np.pi / 2.0 - latitude) * units.rad, 'y')
    
    # Calculate response in earth centered coordinates
    # by rotation of response in coordinates aligned
    # with the detector arms
    resps = []
    vecs = []
    for angle, azi in [(yangle, yaltitude), (xangle, xaltitude)]:
        rm0 = rotation_matrix(angle * units.rad, 'z')
        rmN = rotation_matrix(-azi *  units.rad, 'y')
        rm = rm2 @ rm1 @ rm0 @ rmN
        # apply rotation
        resps.append(rm @ resp @ rm.T / 2.0)
        vecs.append(rm @ np.array([-1, 0, 0]))

    full_resp = (resps[0] - resps[1])
    loc = coordinates.EarthLocation.from_geodetic(longitude * units.rad,
                                                  latitude * units.rad,
                                                  height=height*units.meter)
    loc = np.array([loc.x.value, loc.y.value, loc.z.value])
    _ground_detectors[name] = {'location': loc,
                               'response': full_resp,
                               'xresp': resps[1],
                               'yresp': resps[0],
                               'xvec': vecs[1],
                               'yvec': vecs[0],
                               'yangle': yangle,
                               'xangle': xangle,
                               'height': height,
                               'xaltitude': xaltitude,
                               'yaltitude': yaltitude,
                               'ylength': ylength,
                               'xlength': xlength,
                              }

# Notation matches
# Eq 4 of https://link.aps.org/accepted/10.1103/PhysRevD.96.084004
def single_arm_frequency_response(f, n, arm_length):
    """ The relative amplitude factor of the arm response due to
    signal delay. This is relevant where the long-wavelength
    approximation no longer applies)
    """
    n = np.clip(n, -0.999, 0.999)
    phase = arm_length / constants.c.value * 2.0j * np.pi * f
    a = 1.0 / 4.0 / phase
    b = (1 - np.exp(-phase * (1 - n))) / (1 - n)
    c = np.exp(-2.0 * phase) * (1 - np.exp(phase * (1 + n))) / (1 + n)
    return a * (b - c) * 2.0  # We'll make this relative to the static resp

def load_detector_config(config_files):
    """ Add custom detectors from a configuration file

    Parameters
    ----------
    config_files: str or list of strs
        The config file(s) which specify new detectors
    """
    methods = {'earth_normal': (add_detector_on_earth,
                                ['longitude', 'latitude'])}
    conf = InterpolatingConfigParser(config_files)
    dets = conf.get_subsections('detector')
    for det in dets:
        kwds = dict(conf.items('detector-{}'.format(det)))
        try:
            method, arg_names = methods[kwds.pop('method')]
        except KeyError:
            raise ValueError("Missing or unkown method, "
                             "options are {}".format(methods.keys()))
        for k in kwds:
            kwds[k] = float(kwds[k])
        try:
            args = [kwds.pop(arg) for arg in arg_names]
        except KeyError as e:
            raise ValueError("missing required detector argument"
                             " {} are required".format(arg_names))
        method(det.upper(), *args, **kwds)


# prepopulate using detectors hardcoded into lalsuite
for pref, name in get_available_lal_detectors():
    lalsim = pycbc.libutils.import_optional('lalsimulation')
    lal_det = lalsim.DetectorPrefixToLALDetector(pref).frDetector
    add_detector_on_earth(pref,
                          lal_det.vertexLongitudeRadians,
                          lal_det.vertexLatitudeRadians,
                          height=lal_det.vertexElevation,
                          xangle=lal_det.xArmAzimuthRadians,
                          yangle=lal_det.yArmAzimuthRadians,
                          xlength=lal_det.xArmMidpoint * 2,
                          ylength=lal_det.yArmMidpoint * 2,
                          xaltitude=lal_det.xArmAltitudeRadians,
                          yaltitude=lal_det.yArmAltitudeRadians,
                          )

# autoload detector config files
if 'PYCBC_DETECTOR_CONFIG' in os.environ:
    load_detector_config(os.environ['PYCBC_DETECTOR_CONFIG'].split(':'))


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
        
        lal_detectors = [pfx for pfx, name in get_available_lal_detectors()]
        if detector_name in _ground_detectors:
            self.info = _ground_detectors[detector_name]
            self.response = self.info['response']
            self.location = self.info['location']
        else:
            raise ValueError("Unkown detector {}".format(detector_name))

        loc = coordinates.EarthLocation(self.location[0],
                                        self.location[1],
                                        self.location[2],
                                        unit=meter)
        self.latitude = loc.lat.rad
        self.longitude = loc.lon.rad

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

    def lal(self):
        """ Return lal data type detector instance """
        import lal
        d = lal.FrDetector()
        d.vertexLongitudeRadians = self.longitude
        d.vertexLatitudeRadians = self.latitude
        d.vertexElevation = self.info['height']
        d.xArmAzimuthRadians = self.info['xangle']
        d.yArmAzimuthRadians = self.info['yangle']
        d.xArmAltitudeRadians = self.info['xaltitude']
        d.yArmAltitudeRadians = self.info['yaltitude']

        # This is somewhat abused by lalsimulation at the moment
        # to determine a filter kernel size. We set this only so that
        # value gets a similar number of samples as other detectors
        # it is used for nothing else
        d.yArmMidpoint = self.info['ylength'] / 2.0
        d.xArmMidpoint = self.info['xlength'] / 2.0

        x = lal.Detector()
        r = lal.CreateDetector(x, d, lal.LALDETECTORTYPE_IFODIFF)
        self._lal = r
        return r

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

    def antenna_pattern(self, right_ascension, declination, polarization, t_gps,
                        frequency=0,
                        polarization_type='tensor'):
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

        if frequency:
            e0 = cosdec * cosgha
            e1 = cosdec * -singha
            e2 = sin(declination)
            nhat = np.array([e0, e1, e2], dtype=object)

            nx = nhat.dot(self.info['xvec'])
            ny = nhat.dot(self.info['yvec'])

            rx = single_arm_frequency_response(frequency, nx,
                                               self.info['xlength'])
            ry = single_arm_frequency_response(frequency, ny,
                                               self.info['ylength'])
            resp = ry * self.info['yresp'] -  rx * self.info['xresp']
            ttype = np.complex128
        else:
            resp = self.response
            ttype = np.float64

        x0 = -cospsi * singha - sinpsi * cosgha * sindec
        x1 = -cospsi * cosgha + sinpsi * singha * sindec
        x2 =  sinpsi * cosdec

        x = np.array([x0, x1, x2], dtype=object)
        dx = resp.dot(x)

        y0 =  sinpsi * singha - cospsi * cosgha * sindec
        y1 =  sinpsi * cosgha + cospsi * singha * sindec
        y2 =  cospsi * cosdec

        y = np.array([y0, y1, y2], dtype=object)
        dy = resp.dot(y)

        if polarization_type != 'tensor':
            z0 = -cosdec * cosgha
            z1 = cosdec * singha
            z2 = -sindec
            z = np.array([z0, z1, z2], dtype=object)
            dz = resp.dot(z)

        if polarization_type == 'tensor':
            if hasattr(dx, 'shape'):
                fplus = (x * dx - y * dy).sum(axis=0).astype(ttype)
                fcross = (x * dy + y * dx).sum(axis=0).astype(ttype)
            else:
                fplus = (x * dx - y * dy).sum()
                fcross = (x * dy + y * dx).sum()
            return fplus, fcross

        elif polarization_type == 'vector':
            if hasattr(dx, 'shape'):
                fx = (z * dx + x * dz).sum(axis=0).astype(ttype)
                fy = (z * dy + y * dz).sum(axis=0).astype(ttype)
            else:
                fx = (z * dx + x * dz).sum()
                fy = (z * dy + y * dz).sum()

            return fx, fy

        elif polarization_type == 'scalar':
            if hasattr(dx, 'shape'):
                fb = (x * dx + y * dy).sum(axis=0).astype(ttype)
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

        ehat = np.array([e0, e1, e2], dtype=object)
        dx = other_location - self.location
        return dx.dot(ehat).astype(np.float64) / constants.c.value

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

    def project_wave(self, hp, hc, ra, dec, polarization,
                     method='lal',
                     reference_time=None):
        """Return the strain of a waveform as measured by the detector.
        Apply the time shift for the given detector relative to the assumed
        geocentric frame and apply the antenna patterns to the plus and cross
        polarizations.

        Parameters
        ----------
        hp: pycbc.types.TimeSeries
            Plus polarization of the GW
        hc: pycbc.types.TimeSeries
            Cross polarization of the GW
        ra: float
            Right ascension of source location
        dec: float
            Declination of source location
        polarization: float
            Polarization angle of the source
        method: {'lal', 'constant', 'vary_polarization'}
            The method to use for projecting the polarizations into the
            detector frame. Default is 'lal'.
        reference_time: float, Optional
            The time to use as, a reference for some methods of projection.
            Used by 'constant' and 'vary_polarization' methods. Uses average
            time if not provided.
        """
        # The robust and most fefature rich method which includes
        # time changing antenna patterns and doppler shifts due to the
        # earth rotation and orbit
        if method == 'lal':
            import lalsimulation
            h_lal = lalsimulation.SimDetectorStrainREAL8TimeSeries(
                    hp.astype(np.float64).lal(), hc.astype(np.float64).lal(),
                    ra, dec, polarization, self.lal())
            ts = TimeSeries(
                    h_lal.data.data, delta_t=h_lal.deltaT, epoch=h_lal.epoch,
                    dtype=np.float64, copy=False)

        # 'constant' assume fixed orientation relative to source over the
        # duration of the signal, accurate for short duration signals
        # 'fixed_polarization' applies only time changing orientation
        # but no doppler corrections
        elif method in ['constant', 'vary_polarization']:
            if reference_time is not None:
                rtime = reference_time
            else:
                # In many cases, one should set the reference time if using
                # this method as we don't know where the signal is within
                # the given time series. If not provided, we'll choose
                # the midpoint time.
                rtime = (float(hp.end_time) + float(hp.start_time)) / 2.0

            if method == 'constant':
                time = rtime
            elif method == 'vary_polarization':
                if (not isinstance(hp, TimeSeries) or
                   not isinstance(hc, TimeSeries)):
                    raise TypeError('Waveform polarizations must be given'
                                    ' as time series for this method')

                # this is more granular than needed, may be optimized later
                # assume earth rotation in ~30 ms needed for earth ceneter
                # to detector is completely negligible.
                time = hp.sample_times.numpy()

            fp, fc = self.antenna_pattern(ra, dec, polarization, time)
            dt = self.time_delay_from_earth_center(ra, dec, rtime)
            ts = fp * hp + fc * hc
            ts.start_time = float(ts.start_time) + dt

        # add in only the correction for the time variance in the polarization
        # due to the earth's rotation, no doppler correction applied
        else:
            raise ValueError("Unkown projection method {}".format(method))
        return ts

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

    def effective_distance(self, distance, ra, dec, pol, time, inclination):
        """ Distance scaled to account for amplitude factors

        The effective distance of the source. This scales the distance so that
        the amplitude is equal to a source which is optimally oriented with
        respect to the detector. For fixed detector-frame intrinsic parameters
        this is a measure of the expected signal strength.

        Parameters
        ----------
        distance: float
            Source luminosity distance in megaparsecs
        ra: float
            The right ascension in radians
        dec: float
            The declination in radians
        pol: float
            Polarization angle of the gravitational wave in radians
        time: float
            GPS time in seconds
        inclination:
            The inclination of the binary's orbital plane

        Returns
        -------
        eff_dist: float
            The effective distance of the source
        """
        fp, fc = self.antenna_pattern(ra, dec, pol, time)
        ic = np.cos(inclination)
        ip = 0.5 * (1. + ic * ic)
        scale = ((fp * ip) ** 2.0 + (fc * ic) ** 2.0) ** 0.5
        return distance / scale

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


"""     Space-based detector class      """

from pycbc.coordinates.space import TIME_OFFSET_20_DEGREES

class space_detector(object):
    """
    Space-based GW detector.
    """
    def __init__(self, detector='LDC', reference_time=None, orbits='EqualArmlength',
                 apply_offset=False, offset=TIME_OFFSET_20_DEGREES, use_gpu=False):
        """
        Parameters
        ----------
        detector: str (optional)
            The detector model used to generate the detector response. Default 'LDC'.
            Accepts the following arguments:

                'LDC': Uses packages written by the LISA Data Challenges working group
                to simulate a LISA-like detector. Specifically, LISA GW Response
                (10.5281/zenodo.6423435) is used to generate the detector response
                (i.e. the GW projected to the six laser links), and pyTDI
                (10.5281/zenodo.6351736) is used to generate the TDI observables.
                Orbits are generated using LISA Orbits (https://pypi.org/project/lisaorbits/).

                'FLR': Uses FastLISAResponse (arXiv:2204.06633) to simulate a LISA-like
                detector. Orbits are generated using LISA Analysis Tools
                (10.5281/zenodo.10930979). This class is GPU-compatible via CuPy.

        reference_time: float-like (optional)
            The reference time of signal in the SSB frame. Accepts any type that is
            castable to float (e.g. LIGOTimeGPS). By default, the reference time is
            set to the midpoint of the time series.

        orbits : str (optional)
            Specify which type of LISA orbit to use for projection. Accepts a filename,
            'EqualArmlength' or 'Keplerian'. If filename, orbital data is read in from
            the specified file path, assumed to be an h5 file with LISA Orbits syntax
            (https://pypi.org/project/lisaorbits/). If 'EqualArmlength', 'Keplerian',
            or 'ESA', a file is generated using LISA Orbits and saved to the working
            directory. Default 'EqualArmlength'.

        apply_offset : bool (optional)
            Specify whether to add a time offset to input times. If True,
            an offset is added to ensure LISA is oriented correctly relative to the
            Earth at t=0. If False, no offset is applied. Default False.

        offset : float (optional)
            Time offset in seconds to apply to input times if apply_offset = True.
            Default pycbc.coordinates.space.TIME_OFFSET_20_DEGREES (places LISA
            ~20 deg behind Earth).

        use_gpu : bool (optional)
            Specify whether to use GPUs in response generation. Default False. Only
            applies to GPU-compatible classes.
        """
        # initialize detector
        accept_dets = ['LDC', 'FLR']
        if detector == 'LDC':
            self.default_orbits = ['EqualArmlength', 'Keplerian']
            try:
                from lisagwresponse import ReadStrain
                ### these are probably temporary
                ### there must be a better way to import all of this
                ### on command
                self.ReadStrain = ReadStrain
            except ImportError:
                raise ImportError('lisagwresponse required to generate projections')
            try:
                from pytdi import michelson, Data
                self.michelson = michelson
                self.Data = Data
            except ImportError:
                raise ImportError('pyTDI required for TDI combinations')
            if orbits in self.default_orbits:
                try:
                    import lisaorbits
                    self.lisaorbits = lisaorbits
                except ImportError:
                    raise ImportError('LISA Orbits required if using a default orbit')
        
        elif detector == 'FLR':
            logging.warning("WARNING: This implementation of FastLISAResponse is a work " +
                            "in progress. Currently not consistent with LDC waveforms.")
            self.default_orbits = ['EqualArmlength', 'ESA']
            try:
                from fastlisaresponse import pyResponseTDI
                self.pyResponseTDI = pyResponseTDI
            except ImportError:
                raise ImportError('FastLISAResponse required for LISA projection/TDI')
            try:
                import lisatools.detector
                self.ltdet = lisatools.detector
            except ImportError:
                raise ImportError('LISAanalysistools required for loading orbit data')

        else:
            raise ValueError('Unrecognized detector argument. Currently only accepts: ' +
                             f'{accept_dets}')

        # specify whether to apply offsets to GPS times
        if apply_offset:
            self.offset = offset
        else:
            self.offset = 0.

        # allocate caches
        self.det = detector
        self.dt = None
        self.sample_times = None
        self.proj_init = None
        self.tdi_init = None
        if reference_time is not None:
            reference_time = float(reference_time)
        self.ref_time = reference_time
        self.start_time = None
        self.orbits = orbits
        self.orbit_start_time = None
        self.orbit_end_time = None
        self.use_gpu = use_gpu
        self.pad_data = False

    def orbits_init(self, orbits, **kwargs):
        """
        Initialize the orbital information for the constellation.

        Parameters
        ----------
        orbits : str
            The type of orbit to read in. If using a default orbit config,
            the corresponding orbital data is either generated or read in. 
            Else, the input is treated as a file path following LISA
            Orbits format. Default "EqualArmlength".

        Keywords
        --------
        new_file : str
            The path of the new file to be generated.

        length : int
            The number of samples to use if generating a new orbit file.
            Default 316; generates ~1 year worth of data.

        dt : float
            The time step in seconds to use if generating a new orbit file.
            Default 100000; generates ~1 year worth of data.

        t_init : float
            The start time in seconds to use if generating a new orbit file.
            Default 0; generates data starting at the LISA mission start.
        """
        if self.det == 'LDC':
            # generate a new file
            if orbits in self.default_orbits:
                if orbits == 'EqualArmlength':
                    o = self.lisaorbits.EqualArmlengthOrbits()
                if orbits == 'Keplerian':
                    o = self.lisaorbits.KeplerianOrbits()

                # read in args for orbit file generation
                gen_args = dict(new_file='orbits.h5', dt=1e5, size=316, t0=0.)
                for key in gen_args.keys():
                    try:
                        gen_args[key] = kwargs[key]
                    except (KeyError, TypeError):
                        pass

                # write to file
                o.write(gen_args['new_file'], gen_args['dt'], gen_args['size'],
                        gen_args['t0'], mode='w')
                self.orbit_start_time = gen_args['t0']
                self.orbit_end_time = gen_args['t0'] + gen_args['size']*gen_args['dt']
                self.orbits = gen_args['new_file']

            # read in from an existing file path
            else:
                import h5py
                with h5py.File(orbits, 'r') as f:
                    self.orbit_start_time = f.attrs['t0']
                    self.orbit_end_time = self.orbit_start_time + \
                                          f.attrs['dt']*f.attrs['size']
                self.orbits = orbits

        if self.det == 'FLR':
            # if orbits are already a class instance, skip this
            ### this would mean reinitializing the class if we want to change the orbits
            if type(self.orbits) is not (str or None):
                return

            # load an orbit from lisatools
            if orbits in self.default_orbits:
                if orbits == 'EqualArmlength':
                    o = self.ltdet.EqualArmlengthOrbits()
                if orbits == 'ESA':
                    o = self.ltdet.ESAOrbits()

            # create a new orbits instance for file input
            else:
                class CustomOrbits(self.ltdet.Orbits):
                    def __init__(self):
                        super().__init__(orbits)
                o = CustomOrbits()

            self.orbits = o
            self.orbit_start_time = self.orbits.t_base[0]
            self.orbit_end_time = self.orbits.t_base[-1]

    def strain_container(self, response, orbits=None):
        """
        Read in the necessary link and orbit information for generating TDI channels
        via pyTDI. Replicates the functionality of pyTDI.Data.from_gws().

        Parameters
        ----------
        response : array
            The laser link projections of the GW. Uses format of self.get_links() output.

        orbits : str, optional
            The path to the file containing orbital information for the LISA constellation.
            Default to class attribute self.orbits.

        Returns
        -------
        dict, array
            The arguments and measurements associated with the link and orbital data.
            Can be passed into pyTDI.michelson.X1.build(**args)(measurements).
        """
        compatible = ['LDC',]
        assert self.det in compatible, f"This class is only compatible with one of: {compatible}"
        
        links = ['12', '23', '31', '13', '32', '21']

        # format the measurements from link data
        measurements = {}
        for i, link in enumerate(links):
            measurements[f'isi_{link}'] = response[:, i]
            measurements[f'isi_sb_{link}'] = response[:, i]
            measurements[f'tmi_{link}'] = 0.
            measurements[f'rfi_{link}'] = 0.
            measurements[f'rfi_sb_{link}'] = 0.

        # collect attributes
        df = 1/self.dt
        t0 = self.orbit_start_time

        # call in the orbital data using pyTDI
        if orbits is None:
            orbits = self.orbits
        return self.Data.from_orbits(orbits, df, t0, 'tcb/ltt', **measurements)

    def apply_polarization(self, hp, hc, polarization):
        """
        Apply polarization rotation matrix.

        Parameters
        ----------
        hp : array
            The plus polarization of the GW.

        hc : array
            The cross polarization of the GW.

        polarization : float
            The SSB polarization angle of the GW in radians.

        Returns
        -------
        (pycbc.types.TimeSeries, pycbc.types.TimeSeries)
            The plus and cross polarizations of the GW rotated by the polarization angle.
        """
        cphi = cos(2*polarization)
        sphi = sin(2*polarization)

        hp_ssb = hp*cphi - hc*sphi
        hc_ssb = hp*sphi + hc*cphi

        return hp_ssb, hc_ssb
    
    def get_links(self, hp, hc, lamb, beta, polarization=0, reference_time=None):
        """
        Project a radiation frame waveform to the LISA constellation.

        Parameters
        ----------
        hp : pycbc.types.TimeSeries
            The plus polarization of the GW in the radiation frame.

        hc : pycbc.types.TimeSeries
            The cross polarization of the GW in the radiation frame.

        lamb : float
            The ecliptic longitude of the source in the SSB frame.

        beta : float
            The ecliptic latitude of the source in the SSB frame.

        polarization : float (optional)
            The polarization angle of the GW in radians. Default 0.

        reference_time : float (optional)
            The reference time of the GW signal in the SSB frame. Default
            behavior places start time of the signal at GPS time t=0.

        Returns
        -------
        ndarray
            The waveform projected to the LISA laser links. Shape is (6, N)
            for input waveforms with N total samples.
        """
        # rotate GW from radiation frame to SSB using polarization angle
        hp, hc = self.apply_polarization(hp, hc, polarization)
        
        # get dt and time series from waveform data
        if self.dt is None:
            self.dt = hp.delta_t

        # call the orbital data
        self.orbits_init(orbits=self.orbits)

        # calculate max/min orbital times with buffers for light travel time
        start_buffer = self.orbit_start_time
        end_buffer = self.orbit_end_time
        if self.det == 'LDC':
            lisa_arm = 2.5e9 # nominal LISA arm length (m)
            ltt_au = constants.au.value / constants.c.value
            ltt_arm = lisa_arm / constants.c.value
            start_buffer += ltt_arm + ltt_au
            end_buffer += ltt_au
        
        # make sure signal lies within orbit length
        if hp.duration + hp.start_time >= end_buffer:
            logging.warning("Time of signal end is greater than end of orbital data. " +
                            "Cutting signal at orbit end time.")
            # cut off data succeeding orbit end time
            orbit_end_idx = np.argwhere(hp.sample_times.numpy() <= end_buffer)[-1][0]
            hp = hp[:orbit_end_idx]
            hc = hc[:orbit_end_idx]
            
        if hp.start_time < start_buffer:
            logging.warning("Time of signal start is less than start of orbital data. " + 
                            "Cutting signal at orbit start time.")
            # cut off data preceding orbit start time
            orbit_start_idx = np.argwhere(hp.sample_times.numpy() >= start_buffer)[0][0]
            hp = hp[orbit_start_idx:]
            hc = hc[orbit_start_idx:]
            
            # update start time if truncating
            self.start_time = start_buffer - self.offset
            if self.pad_data:
                self.start_time += self.t0

        # configure the orbit to match signal
        self.sample_times = hp.sample_times.numpy()

        if self.det == 'LDC':
            if self.proj_init is None:
                # initialize the class
                self.proj_init = self.ReadStrain(self.sample_times, hp, hc,
                                                 gw_beta=beta, gw_lambda=lamb,
                                                 orbits=self.orbits)
            else:
                # update params in the initialized class
                self.proj_init.gw_beta = beta
                self.proj_init.gw_lambda = lamb
                self.proj_init.set_strain(self.sample_times, hp, hc)

            # project the signal
            wf_proj = self.proj_init.compute_gw_response(self.sample_times, self.proj_init.LINKS)

        if self.det == 'FLR':
            # initialize orbital data
            self.orbits.configure(t_arr = self.sample_times)
            
            # format wf to hp + i*hc
            hp = hp.numpy()
            hc = hc.numpy()
            wf = hp + 1j*hc

            ### TODO: conversions to/from CuPy arrays if use_gpu
            
            if self.proj_init is None:
                # initialize the class
                self.proj_init = self.pyResponseTDI(1/self.dt, len(wf), orbits=self.orbits,
                                                    use_gpu=self.use_gpu)
            else:
                # update params in the initialized class
                self.proj_init.sampling_frequency = 1/self.dt
                self.proj_init.num_pts = len(wf)
                self.proj_init.orbits = self.orbits
                self.proj_init.use_gpu = self.use_gpu

            # project the signal
            self.proj_init.get_projections(wf, lamb, beta, t0=self.t0)
            wf_proj = self.proj_init.y_gw
        
        return wf_proj

    def project_wave(self, hp, hc, lamb, beta, polarization, reference_time=None,
                     tdi=1, tdi_chan='AET', pad_data=False, remove_garbage=True,
                     t0=10000.):
        """
        Evaluate the TDI observables.

        In practice, TDI channels will have a "warmup time" associated with the
        light travel times between frames and satellites. By default, 10000 s
        worth of data will be excised from the start and end of the output TDI
        channels to remove glitches due to this effect. Thus, the output channels
        will be 20000 seconds shorter than the input waveforms by default. See
        pad_data, remove_garbage, and t0 to modify this behavior.

        Parameters
        ----------
        hp : pycbc.types.TimeSeries
            The plus polarization of the GW in the radiation frame.

        hc : pycbc.types.TimeSeries
            The cross polarization of the GW in the radiation frame.

        lamb : float
            The ecliptic longitude in the SSB frame.

        beta : float
            The ecliptic latitude in the SSB frame.

        polarization : float
            The polarization angle of the GW in radians.

        reference_time : float (optional)
            The reference time of the GW signal in the SSB frame. Default
            behavior places start time of the signal at GPS time t=0.

        tdi : int (optional)
            TDI channel configuration. Accepts 1 for 1st generation TDI or
            2 for 2nd generation TDI. Default 1.

        tdi_chan : str (optional)
            The TDI observables to calculate. Accepts 'XYZ', 'AET', 'AE'.
            Default 'AET'.
        
        pad_data : bool (optional)
            Flag whether to pad the start and end of the input data with zeros.
            If True, a time length t0 will be added to the start/end of the input.
            Default False.

        remove_garbage : bool or str (optional)
            Flag whether to cut data from the start and end of the output
            TDI channels. If True, a time length t0 will be excised from the start/end
            of the output. If 'zero', a time length t0 at the start/end will be zeroed;
            the total length of the waveform will not change. Default True.

        t0 : float (optional)
            Time length in seconds to pad/remove from data depending on the inputs
            pad_data and remove_garbage. Default 10000.

        Returns
        -------
        {str: pycbc.types.TimeSeries}
            The TDI observables as TimeSeries objects keyed by their
            corresponding TDI channel name.
        """
        self.pad_data = pad_data
        self.dt = hp.delta_t

        # get index corresponding to time length t0
        if t0 is not None:
            self.t0 = t0
        if pad_data or remove_garbage:
            global pad_idx
            pad_idx = int(t0/self.dt)
        
        # get reference time from class
        if reference_time is None:
            if self.ref_time is None:
                # take ref time as midpoint of signal, start time as t=0
                self.ref_time = float(hp.duration/2)
            reference_time = self.ref_time

        # get times corresponding to unpadded data
        self.ref_time = reference_time
        # set ref time at t = 0 in input waveform
        self.start_time = float(self.ref_time + hp.start_time)

        # apply times to wfs
        hp.start_time = self.start_time + self.offset
        hc.start_time = self.start_time + self.offset

        # pad the data with zeros in the SSB frame
        if pad_data:
            hp.prepend_zeros(pad_idx)
            hp.append_zeros(pad_idx)
            hc.prepend_zeros(pad_idx)
            hc.append_zeros(pad_idx)

        # generate the Doppler time series
        response = self.get_links(hp, hc, lamb, beta, polarization=polarization,
                                  reference_time=reference_time)
        
        if self.det == 'LDC':
            # set TDI generation
            if tdi == 1:
                X, Y, Z = self.michelson.X1, self.michelson.Y1, self.michelson.Z1
            elif tdi == 2:
                X, Y, Z = self.michelson.X2, self.michelson.Y2, self.michelson.Z2
            else:
                raise ValueError('Unrecognized TDI generation input. ' +
                                 'Please input either 1 or 2.')

            # load strains and orbits into a Data instance
            self.tdi_init = self.strain_container(response, self.orbits)

            # generate the XYZ TDI channels
            chanx = X.build(**self.tdi_init.args)(self.tdi_init.measurements)
            chany = Y.build(**self.tdi_init.args)(self.tdi_init.measurements)
            chanz = Z.build(**self.tdi_init.args)(self.tdi_init.measurements)

            # convert to AET if specified
            if tdi_chan == 'XYZ':
                tdi_dict = {'X': chanx, 'Y': chany, 'Z': chanz}
            elif tdi_chan == 'AET':
                chana = (chanz - chanx)/np.sqrt(2)
                chane = (chanx - 2*chany + chanz)/np.sqrt(6)
                chant = (chanx + chany + chanz)/np.sqrt(3)
                tdi_dict = {'A': chana, 'E': chane, 'T': chant}
            else:
                raise ValueError('Unrecognized TDI channel input. ' +
                                 'Please input either "XYZ" or "AET".')

        if self.det == 'FLR':
            # set TDI generation
            if tdi == 1:
                tdi_gen = '1st generation'
            elif tdi == 2:
                tdi_gen = '2nd generation'
            else:
                raise ValueError('Unrecognized TDI generation input. ' +
                                 'Please input either 1 or 2.')

            # set TDI channels
            accept_chans = ['XYZ', 'AET', 'AE']
            if tdi_chan not in accept_chans:
                raise ValueError('Unrecognized TDI channel input. ' +
                                 f'Please input one of: {accept_chans}.')

            # copy over the class initialized in get_links
            self.tdi_init = self.proj_init
            self.tdi_init.tdi = tdi_gen
            self.tdi_init.tdi_chan = tdi_chan

            # generate the TDI channels
            tdi_obs = self.tdi_init.get_tdi_delays()
            tdi_dict = {tdi_chan[i]: tdi_obs[i] for i in range(len(tdi_chan))}

        # postprocessing
        for i in range(len(tdi_chan)):
            # treat start and end gaps
            if remove_garbage:
                if remove_garbage == 'zero':
                    # zero the edge data
                    tdi_dict[tdi_chan[i]][:pad_idx] = 0
                    tdi_dict[tdi_chan[i]][-pad_idx:] = 0
                elif type(remove_garbage) == bool:
                    # cut the edge data
                    slc = slice(pad_idx, -pad_idx)
                    tdi_dict[tdi_chan[i]] = tdi_dict[tdi_chan[i]][slc]
                else:
                    raise ValueError('remove_garbage arg must be a bool ' +
                                     'or "zero"')
                
            # save as TimeSeries with SSB times
            tdi_dict[tdi_chan[i]] = TimeSeries(tdi_dict[tdi_chan[i]], delta_t=self.dt,
                                               epoch=self.start_time)
            if pad_data and (not remove_garbage or remove_garbage == 'zero'):
                # scale the start since the pads haven't been removed
                tdi_dict[tdi_chan[i]].start_time -= t0

        return tdi_dict


def ppdets(ifos, separator=', '):
    """Pretty-print a list (or set) of detectors: return a string listing
    the given detectors alphabetically and separated by the given string
    (comma by default).
    """
    if ifos:
        return separator.join(sorted(ifos))
    return 'no detectors'
