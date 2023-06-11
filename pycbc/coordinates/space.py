# Copyright (C) 2023  Shichao Wu, Alex Nitz
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
This module provides coordinate transformations related to space-borne
detectors, such as coordinate transformations between space-borne detectors
and ground-based detectors.
"""

import numpy as np
from scipy.spatial.transform import Rotation
from scipy.optimize import fsolve
from astropy import units
from astropy.constants import c, au
from astropy.time import Time
from astropy.coordinates import BarycentricMeanEcliptic, PrecessedGeocentric
from astropy.coordinates import get_body_barycentric
from astropy.coordinates import SkyCoord


# This constant makes sure LISA is behind the Earth by 19-23 degrees.
# Making this a stand-alone constant will also make it callable by
# the waveform plugin and PE config file. In the unit of 's'.
TIME_OFFSET_20_DEGREES = 7365189.431698299

# "rotation_matrix_ssb_to_lisa" and "lisa_position_ssb" should be
# more general for other detectors in the near future.


def rotation_matrix_ssb_to_lisa(alpha):
    """ The rotation matrix (of frame basis) from SSB frame to LISA frame.
    This function assumes the angle between LISA plane and the ecliptic
    is 60 degrees, and the period of LISA's self-rotation and orbital
    revolution is both one year.

    Parameters
    ----------
    alpha : float
        The angular displacement of LISA in SSB frame.
        In the unit of 'radian'.

    Returns
    -------
    r_total : numpy.array
        A 3x3 rotation matrix from SSB frame to LISA frame.
    """
    r = Rotation.from_rotvec([
        [0, 0, alpha],
        [0, -np.pi/3, 0],
        [0, 0, -alpha]
    ]).as_matrix()
    r_total = np.array(r[0]) @ np.array(r[1]) @ np.array(r[2])

    return r_total


def lisa_position_ssb(t_lisa, t0=TIME_OFFSET_20_DEGREES):
    """ Calculating the position vector and angular displacement of LISA
    in the SSB frame, at a given time. This function assumes LISA's barycenter
    is orbiting around a circular orbit within the ecliptic behind the Earth.
    The period of it is one year.

    Parameters
    ----------
    t_lisa : float
        The time when a GW signal arrives at the origin of LISA frame,
        or any other time you want.
    t0 : float
        The initial time offset of LISA, in the unit of 's',
        default is 7365189.431698299. This makes sure LISA is behind
        the Earth by 19-23 degrees.

    Returns
    -------
    (p, alpha) : tuple
    p : numpy.array
        The position vector of LISA in the SSB frame. In the unit of 'm'.
    alpha : float
        The angular displacement of LISA in the SSB frame.
        In the unit of 'radian'.
    """
    OMEGA_0 = 1.99098659277e-7
    R_ORBIT = au.value
    alpha = np.mod(OMEGA_0 * (t_lisa + t0), 2*np.pi)
    p = np.array([[R_ORBIT * np.cos(alpha)],
                  [R_ORBIT * np.sin(alpha)],
                  [0]], dtype=object)
    return (p, alpha)


def localization_to_propagation_vector(longitude, latitude):
    """ Converting the sky localization to the corresponding
    propagation unit vector of a GW signal.

    Parameters
    ----------
    longitude : float
        The longitude, in the unit of 'radian'.
    latitude : float
        The latitude, in the unit of 'radian'.

    Returns
    -------
    [[x], [y], [z]] : numpy.array
        The propagation unit vector of that GW signal.
    """
    x = -np.cos(latitude) * np.cos(longitude)
    y = -np.cos(latitude) * np.sin(longitude)
    z = -np.sin(latitude)

    return np.array([[x], [y], [z]])


def propagation_vector_to_localization(k):
    """ Converting the propagation unit vector to the corresponding
    sky localization of a GW signal.

    Parameters
    ----------
    k : numpy.array
        The propagation unit vector of a GW signal.

    Returns
    -------
    (longitude, latitude) : tuple
        The sky localization of that GW signal.
    """
    # latitude already within [-pi/2, pi/2]
    latitude = np.float64(np.arcsin(-k[2]))
    longitude = np.float64(np.arctan2(-k[1]/np.cos(latitude), -k[0]/np.cos(latitude)))
    # longitude should within [0, 2*pi)
    longitude = np.mod(longitude, 2*np.pi)

    return (longitude, latitude)


def polarization_newframe(polarization, k, rotation_matrix):
    """ Converting a polarization angle from a frame to a new frame
    by using rotation matrix method.

    Parameters
    ----------
    polarization : float
        The polarization angle in the old frame, in the unit of 'radian'.
    k : numpy.array
        The propagation unit vector of a GW signal in the old frame.
    rotation_matrix : numpy.array
        The rotation matrix (of frame basis) from the old frame to
        the new frame.

    Returns
    -------
    polarization_newframe : float
        The polarization angle in the new frame of that GW signal.
    """
    longitude, _ = propagation_vector_to_localization(k)
    u = np.array([[np.sin(longitude)], [-np.cos(longitude)], [0]])
    rotation_vector = polarization * k
    rotation_polarization = Rotation.from_rotvec(rotation_vector.T[0])
    p = rotation_polarization.apply(u.T[0]).reshape(3, 1)
    p_newframe = rotation_matrix.T @ p
    k_newframe = rotation_matrix.T @ k
    longitude_newframe, latitude_newframe = \
        propagation_vector_to_localization(k_newframe)
    u_newframe = np.array([[np.sin(longitude_newframe)],
                           [-np.cos(longitude_newframe)], [0]])
    v_newframe = np.array([[-np.sin(latitude_newframe) * np.cos(longitude_newframe)],
                           [-np.sin(latitude_newframe) * np.sin(longitude_newframe)],
                           [np.cos(latitude_newframe)]])
    p_dot_u_newframe = np.vdot(p_newframe, u_newframe)
    p_dot_v_newframe = np.vdot(p_newframe, v_newframe)
    polarization_newframe = np.arctan2(p_dot_v_newframe, p_dot_u_newframe)
    polarization_newframe = np.mod(polarization_newframe, 2*np.pi)
    # avoid the round error
    if polarization_newframe == 2*np.pi:
        polarization_newframe = 0

    return polarization_newframe


def t_lisa_from_ssb(t_ssb, longitude_ssb, latitude_ssb, t0=TIME_OFFSET_20_DEGREES):
    """ Calculating the time when a GW signal arrives at the barycenter
    of LISA, by using the time and sky localization in SSB frame.

    Parameters
    ----------
    t_ssb : float
        The time when a GW signal arrives at the origin of SSB frame.
        In the unit of 's'.
    longitude_ssb : float
        The ecliptic longitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    latitude_ssb : float
        The ecliptic latitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    t0 : float
        The initial time offset of LISA, in the unit of 's',
        default is 7365189.431698299. This makes sure LISA is behind
        the Earth by 19-23 degrees.

    Returns
    -------
    t_lisa : float
        The time when a GW signal arrives at the origin of LISA frame.
    """
    k = localization_to_propagation_vector(longitude_ssb, latitude_ssb)

    def equation(t_lisa):
        # LISA is moving, when GW arrives at LISA center,
        # time is t_lisa, not t_ssb.
        p = lisa_position_ssb(t_lisa, t0)[0]
        return t_lisa - t_ssb - np.vdot(k, p) / c.value

    return fsolve(equation, t_ssb)[0]


def t_ssb_from_t_lisa(t_lisa, longitude_ssb, latitude_ssb, t0=TIME_OFFSET_20_DEGREES):
    """ Calculating the time when a GW signal arrives at the barycenter
    of SSB, by using the time in LISA frame and sky localization in SSB frame.

    Parameters
    ----------
    t_lisa : float
        The time when a GW signal arrives at the origin of LISA frame.
        In the unit of 's'.
    longitude_ssb : float
        The ecliptic longitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    latitude_ssb : float
        The ecliptic latitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    t0 : float
        The initial time offset of LISA, in the unit of 's',
        default is 7365189.431698299. This makes sure LISA is behind
        the Earth by 19-23 degrees.

    Returns
    -------
    t_ssb : float
        The time when a GW signal arrives at the origin of SSB frame.
    """
    k = localization_to_propagation_vector(longitude_ssb, latitude_ssb)
    # LISA is moving, when GW arrives at LISA center,
    # time is t_lisa, not t_ssb.
    p = lisa_position_ssb(t_lisa, t0)[0]

    def equation(t_ssb):
        return t_lisa - t_ssb - np.vdot(k, p) / c.value

    return fsolve(equation, t_lisa)[0]


def ssb_to_lisa(t_ssb, longitude_ssb, latitude_ssb, polarization_ssb, t0=TIME_OFFSET_20_DEGREES):
    """ Converting the arrive time, the sky localization, and the polarization
    from the SSB frame to the LISA frame.

    Parameters
    ----------
    t_ssb : float or numpy.array
        The time when a GW signal arrives at the origin of SSB frame.
        In the unit of 's'.
    longitude_ssb : float or numpy.array
        The ecliptic longitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    latitude_ssb : float or numpy.array
        The ecliptic latitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    polarization_ssb : float or numpy.array
        The polarization angle of a GW signal in SSB frame.
        In the unit of 'radian'.
    t0 : float
        The initial time offset of LISA, in the unit of 's',
        default is 7365189.431698299. This makes sure LISA is behind
        the Earth by 19-23 degrees.

    Returns
    -------
    (t_lisa, longitude_lisa, latitude_lisa, polarization_lisa) : tuple
    t_lisa : float or numpy.array
        The time when a GW signal arrives at the origin of LISA frame.
        In the unit of 's'.
    longitude_lisa : float or numpy.array
        The longitude of a GW signal in LISA frame, in the unit of 'radian'.
    latitude_lisa : float or numpy.array
        The latitude of a GW signal in LISA frame, in the unit of 'radian'.
    polarization_lisa : float or numpy.array
        The polarization angle of a GW signal in LISA frame.
        In the unit of 'radian'.
    """
    if isinstance(t_ssb, np.ndarray):
        t_ssb_array = t_ssb.copy()
    else:
        t_ssb_array = np.array([t_ssb])
    if isinstance(longitude_ssb, np.ndarray):
        longitude_ssb_array = longitude_ssb.copy()
    else:
        longitude_ssb_array = np.array([longitude_ssb])
    if isinstance(latitude_ssb, np.ndarray):
        latitude_ssb_array = latitude_ssb.copy()
    else:
        latitude_ssb_array = np.array([latitude_ssb])
    if isinstance(polarization_ssb, np.ndarray):
        polarization_ssb_array = polarization_ssb.copy()
    else:
        polarization_ssb_array = np.array([polarization_ssb])
    num = len(t_ssb_array)
    t_lisa_array, longitude_lisa_array = np.zeros(num), np.zeros(num)
    latitude_lisa_array, polarization_lisa_array = np.zeros(num), np.zeros(num)

    for i in range(num):
        t_ssb = t_ssb_array[i]
        longitude_ssb = longitude_ssb_array[i]
        latitude_ssb = latitude_ssb_array[i]
        polarization_ssb = polarization_ssb_array[i]
        if longitude_ssb < 0 or longitude_ssb >= 2*np.pi:
            raise ValueError("Longitude should within [0, 2*pi).")
        if latitude_ssb < -np.pi/2 or latitude_ssb > np.pi/2:
            raise ValueError("Latitude should within [-pi/2, pi/2].")
        if polarization_ssb < 0 or polarization_ssb >= 2*np.pi:
            raise ValueError("Polarization angle should within [0, 2*pi).")
        t_lisa = t_lisa_from_ssb(t_ssb, longitude_ssb, latitude_ssb, t0)
        k_ssb = localization_to_propagation_vector(longitude_ssb, latitude_ssb)
        # Although t_lisa calculated above using the corrected LISA position vector by
        # adding t0, it corresponds to the true t_ssb, not t_ssb+t0,
        # we need to include t0 again to correct LISA position.
        alpha = lisa_position_ssb(t_lisa, t0)[1]
        rotation_matrix_lisa = rotation_matrix_ssb_to_lisa(alpha)
        k_lisa = rotation_matrix_lisa.T @ k_ssb
        longitude_lisa, latitude_lisa = propagation_vector_to_localization(k_lisa)
        polarization_lisa = polarization_newframe(polarization_ssb, k_ssb, rotation_matrix_lisa)
        t_lisa_array[i] = t_lisa
        longitude_lisa_array[i] = longitude_lisa
        latitude_lisa_array[i] = latitude_lisa
        polarization_lisa_array[i] = polarization_lisa
    if num == 1:
        params_lisa = (t_lisa_array[0], longitude_lisa_array[0],
                       latitude_lisa_array[0], polarization_lisa_array[0])
    else:
        params_lisa = (t_lisa_array, longitude_lisa_array,
                       latitude_lisa_array, polarization_lisa_array)

    return params_lisa


def lisa_to_ssb(t_lisa, longitude_lisa, latitude_lisa, polarization_lisa, t0=TIME_OFFSET_20_DEGREES):
    """ Converting the arrive time, the sky localization, and the polarization
    from the LISA frame to the SSB frame.

    Parameters
    ----------
    t_lisa : float or numpy.array
        The time when a GW signal arrives at the origin of LISA frame.
        In the unit of 's'.
    longitude_lisa : float or numpy.array
        The longitude of a GW signal in LISA frame, in the unit of 'radian'.
    latitude_lisa : float or numpy.array
        The latitude of a GW signal in LISA frame, in the unit of 'radian'.
    polarization_lisa : float or numpy.array
        The polarization angle of a GW signal in LISA frame.
        In the unit of 'radian'.
    t0 : float
        The initial time offset of LISA, in the unit of 's',
        default is 7365189.431698299. This makes sure LISA is behind
        the Earth by 19-23 degrees.

    Returns
    -------
    (t_ssb, longitude_ssb, latitude_ssb, polarization_ssb) : tuple
    t_ssb : float or numpy.array
        The time when a GW signal arrives at the origin of SSB frame.
        In the unit of 's'.
    longitude_ssb : float or numpy.array
        The ecliptic longitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    latitude_ssb : float or numpy.array
        The ecliptic latitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    polarization_ssb : float or numpy.array
        The polarization angle of a GW signal in SSB frame.
        In the unit of 'radian'.
    """
    if isinstance(t_lisa, np.ndarray):
        t_lisa_array = t_lisa.copy()
    else:
        t_lisa_array = np.array([t_lisa])
    if isinstance(longitude_lisa, np.ndarray):
        longitude_lisa_array = longitude_lisa.copy()
    else:
        longitude_lisa_array = np.array([longitude_lisa])
    if isinstance(latitude_lisa, np.ndarray):
        latitude_lisa_array = latitude_lisa.copy()
    else:
        latitude_lisa_array = np.array([latitude_lisa])
    if isinstance(polarization_lisa, np.ndarray):
        polarization_lisa_array = polarization_lisa.copy()
    else:
        polarization_lisa_array = np.array([polarization_lisa])
    num = len(t_lisa_array)
    t_ssb_array, longitude_ssb_array = np.zeros(num), np.zeros(num)
    latitude_ssb_array, polarization_ssb_array = np.zeros(num), np.zeros(num)

    for i in range(num):
        t_lisa = t_lisa_array[i]
        longitude_lisa = longitude_lisa_array[i]
        latitude_lisa = latitude_lisa_array[i]
        polarization_lisa = polarization_lisa_array[i]
        if longitude_lisa < 0 or longitude_lisa >= 2*np.pi:
            raise ValueError("Longitude should within [0, 2*pi).")
        if latitude_lisa < -np.pi/2 or latitude_lisa > np.pi/2:
            raise ValueError("Latitude should within [-pi/2, pi/2].")
        if polarization_lisa < 0 or polarization_lisa >= 2*np.pi:
            raise ValueError("Polarization angle should within [0, 2*pi).")
        k_lisa = localization_to_propagation_vector(longitude_lisa, latitude_lisa)
        alpha = lisa_position_ssb(t_lisa, t0)[1]
        rotation_matrix_lisa = rotation_matrix_ssb_to_lisa(alpha)
        k_ssb = rotation_matrix_lisa @ k_lisa
        longitude_ssb, latitude_ssb = propagation_vector_to_localization(k_ssb)
        t_ssb = t_ssb_from_t_lisa(t_lisa, longitude_ssb, latitude_ssb, t0)
        polarization_ssb = polarization_newframe(polarization_lisa, k_lisa, rotation_matrix_lisa.T)
        t_ssb_array[i] = t_ssb
        longitude_ssb_array[i] = longitude_ssb
        latitude_ssb_array[i] = latitude_ssb
        polarization_ssb_array[i] = polarization_ssb
    if num == 1:
        params_ssb = (t_ssb_array[0], longitude_ssb_array[0],
                      latitude_ssb_array[0], polarization_ssb_array[0])
    else:
        params_ssb = (t_ssb_array, longitude_ssb_array,
                      latitude_ssb_array, polarization_ssb_array)

    return params_ssb


def rotation_matrix_ssb_to_geo(epsilon=np.deg2rad(23.439281)):
    """ The rotation matrix (of frame basis) from SSB frame to
    geocentric frame.

    Parameters
    ----------
    epsilon : float
        The Earth's axial tilt (obliquity), in the unit of 'radian'.

    Returns
    -------
    r : numpy.array
        A 3x3 rotation matrix from SSB frame to geocentric frame.
    """
    r = Rotation.from_rotvec([
        [-epsilon, 0, 0]
    ]).as_matrix()

    return np.array(r[0])


def earth_position_ssb(t_geo):
    """ Calculating the position vector and angular displacement of the Earth
    in the SSB frame, at a given time. By using Astropy.

    Parameters
    ----------
    t_geo : float
        The time when a GW signal arrives at the origin of geocentric frame,
        or any other time you want.

    Returns
    -------
    (p, alpha) : tuple
    p : numpy.array
        The position vector of the Earth in the SSB frame. In the unit of 'm'.
    alpha : float
        The angular displacement of the Earth in the SSB frame.
        In the unit of 'radian'.
    """
    t = Time(t_geo, format='gps')
    pos = get_body_barycentric('earth', t)
    # BarycentricMeanEcliptic doesn't have obstime attribute,
    # it's a good inertial frame, but ICRS is not.
    icrs_coord = SkyCoord(pos, frame='icrs', obstime=t)
    bme_coord = icrs_coord.transform_to(
                    BarycentricMeanEcliptic(equinox='J2000'))
    x = bme_coord.cartesian.x.to(units.m).value
    y = bme_coord.cartesian.y.to(units.m).value
    z = bme_coord.cartesian.z.to(units.m).value
    p = np.array([[x], [y], [z]])
    alpha = bme_coord.lon.rad

    return (p, alpha)


def t_geo_from_ssb(t_ssb, longitude_ssb, latitude_ssb):
    """ Calculating the time when a GW signal arrives at the barycenter
    of the Earth, by using the time and sky localization in SSB frame.

    Parameters
    ----------
    t_ssb : float
        The time when a GW signal arrives at the origin of SSB frame.
        In the unit of 's'.
    longitude_ssb : float
        The ecliptic longitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    latitude_ssb : float
        The ecliptic latitude of a GW signal in SSB frame.
        In the unit of 'radian'.

    Returns
    -------
    t_geo : float
        The time when a GW signal arrives at the origin of geocentric frame.
    """
    k = localization_to_propagation_vector(longitude_ssb, latitude_ssb)

    def equation(t_geo):
        # Earth is moving, when GW arrives at Earth center,
        # time is t_geo, not t_ssb.
        p = earth_position_ssb(t_geo)[0]
        return t_geo - t_ssb - np.vdot(k, p) / c.value

    return fsolve(equation, t_ssb)[0]


def t_ssb_from_t_geo(t_geo, longitude_ssb, latitude_ssb):
    """ Calculating the time when a GW signal arrives at the barycenter
    of SSB, by using the time in geocentric frame and sky localization
    in SSB frame.

    Parameters
    ----------
    t_geo : float
        The time when a GW signal arrives at the origin of geocentric frame.
        In the unit of 's'.
    longitude_ssb : float
        The ecliptic longitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    latitude_ssb : float
        The ecliptic latitude of a GW signal in SSB frame.
        In the unit of 'radian'.

    Returns
    -------
    t_ssb : float
        The time when a GW signal arrives at the origin of SSB frame.
    """
    k = localization_to_propagation_vector(longitude_ssb, latitude_ssb)
    # Earth is moving, when GW arrives at Earth center,
    # time is t_geo, not t_ssb.
    p = earth_position_ssb(t_geo)[0]

    def equation(t_ssb):
        return t_geo - t_ssb - np.vdot(k, p) / c.value

    return fsolve(equation, t_geo)[0]


def ssb_to_geo(t_ssb, longitude_ssb, latitude_ssb, polarization_ssb, use_astropy=True):
    """ Converting the arrive time, the sky localization, and the polarization
    from the SSB frame to the geocentric frame.

    Parameters
    ----------
    t_ssb : float or numpy.array
        The time when a GW signal arrives at the origin of SSB frame.
        In the unit of 's'.
    longitude_ssb : float or numpy.array
        The ecliptic longitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    latitude_ssb : float or numpy.array
        The ecliptic latitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    polarization_ssb : float or numpy.array
        The polarization angle of a GW signal in SSB frame.
        In the unit of 'radian'.
    use_astropy : bool
        Using Astropy to calculate the sky localization or not.
        Default is True.

    Returns
    -------
    (t_geo, longitude_geo, latitude_geo, polarization_geo) : tuple
    t_geo : float or numpy.array
        The time when a GW signal arrives at the origin of geocentric frame.
        In the unit of 's'.
    longitude_geo : float or numpy.array
        The longitude of a GW signal in geocentric frame.
        In the unit of 'radian'.
    latitude_geo : float or numpy.array
        The latitude of a GW signal in geocentric frame.
        In the unit of 'radian'.
    polarization_geo : float or numpy.array
        The polarization angle of a GW signal in geocentric frame.
        In the unit of 'radian'.
    """
    if isinstance(t_ssb, np.ndarray):
        t_ssb_array = t_ssb.copy()
    else:
        t_ssb_array = np.array([t_ssb])
    if isinstance(longitude_ssb, np.ndarray):
        longitude_ssb_array = longitude_ssb.copy()
    else:
        longitude_ssb_array = np.array([longitude_ssb])
    if isinstance(latitude_ssb, np.ndarray):
        latitude_ssb_array = latitude_ssb.copy()
    else:
        latitude_ssb_array = np.array([latitude_ssb])
    if isinstance(polarization_ssb, np.ndarray):
        polarization_ssb_array = polarization_ssb.copy()
    else:
        polarization_ssb_array = np.array([polarization_ssb])
    num = len(t_ssb_array)
    t_geo_array, longitude_geo_array = np.zeros(num), np.zeros(num)
    latitude_geo_array, polarization_geo_array = np.zeros(num), np.zeros(num)

    for i in range(num):
        t_ssb = t_ssb_array[i]
        longitude_ssb = longitude_ssb_array[i]
        latitude_ssb = latitude_ssb_array[i]
        polarization_ssb = polarization_ssb_array[i]
        if longitude_ssb < 0 or longitude_ssb >= 2*np.pi:
            raise ValueError("Longitude should within [0, 2*pi).")
        if latitude_ssb < -np.pi/2 or latitude_ssb > np.pi/2:
            raise ValueError("Latitude should within [-pi/2, pi/2].")
        if polarization_ssb < 0 or polarization_ssb >= 2*np.pi:
            raise ValueError("Polarization angle should within [0, 2*pi).")
        t_geo = t_geo_from_ssb(t_ssb, longitude_ssb, latitude_ssb)
        k_ssb = localization_to_propagation_vector(longitude_ssb, latitude_ssb)
        rotation_matrix_geo = rotation_matrix_ssb_to_geo()
        if use_astropy:
            # BarycentricMeanEcliptic doesn't have obstime attribute,
            # it's a good inertial frame, but PrecessedGeocentric is not.
            bme_coord = BarycentricMeanEcliptic(lon=longitude_ssb*units.radian,
                                                lat=latitude_ssb*units.radian,
                                                equinox='J2000')
            geo_sky = bme_coord.transform_to(PrecessedGeocentric(
                equinox='J2000', obstime=Time(t_geo, format='gps')))
            longitude_geo = geo_sky.ra.rad
            latitude_geo = geo_sky.dec.rad
        else:
            k_geo = rotation_matrix_geo.T @ k_ssb
            longitude_geo, latitude_geo = propagation_vector_to_localization(k_geo)
        polarization_geo = polarization_newframe(polarization_ssb, k_ssb, rotation_matrix_geo)
        # As mentioned in LDC manual, the p,q vectors are opposite between 
        # LDC and LAL conventions, see Sec 4.1.5 in <LISA-LCST-SGS-MAN-001>.
        polarization_geo = np.mod(polarization_geo+np.pi, 2*np.pi)
        t_geo_array[i] = t_geo
        longitude_geo_array[i] = longitude_geo
        latitude_geo_array[i] = latitude_geo
        polarization_geo_array[i] = polarization_geo
    if num == 1:
        params_geo = (t_geo_array[0], longitude_geo_array[0],
                      latitude_geo_array[0], polarization_geo_array[0])
    else:
        params_geo = (t_geo_array, longitude_geo_array,
                      latitude_geo_array, polarization_geo_array)

    return params_geo


def geo_to_ssb(t_geo, longitude_geo, latitude_geo, polarization_geo, use_astropy=True):
    """ Converting the arrive time, the sky localization, and the polarization
    from the geocentric frame to the SSB frame.

    Parameters
    ----------
    t_geo : float or numpy.array
        The time when a GW signal arrives at the origin of geocentric frame.
        In the unit of 's'.
    longitude_geo : float or numpy.array
        The longitude of a GW signal in geocentric frame.
        In the unit of 'radian'.
    latitude_geo : float or numpy.array
        The latitude of a GW signal in geocentric frame.
        In the unit of 'radian'.
    polarization_geo : float or numpy.array
        The polarization angle of a GW signal in geocentric frame.
        In the unit of 'radian'.
    use_astropy : bool
        Using Astropy to calculate the sky localization or not.
        Default is True.

    Returns
    -------
    (t_ssb, longitude_ssb, latitude_ssb, polarization_ssb) : tuple
    t_ssb : float or numpy.array
        The time when a GW signal arrives at the origin of SSB frame.
        In the unit of 's'.
    longitude_ssb : float or numpy.array
        The ecliptic longitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    latitude_ssb : float or numpy.array
        The ecliptic latitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    polarization_ssb : float or numpy.array
        The polarization angle of a GW signal in SSB frame.
        In the unit of 'radian'.
    """
    if isinstance(t_geo, np.ndarray):
        t_geo_array = t_geo.copy()
    else:
        t_geo_array = np.array([t_geo])
    if isinstance(longitude_geo, np.ndarray):
        longitude_geo_array = longitude_geo.copy()
    else:
        longitude_geo_array = np.array([longitude_geo])
    if isinstance(latitude_geo, np.ndarray):
        latitude_geo_array = latitude_geo.copy()
    else:
        latitude_geo_array = np.array([latitude_geo])
    if isinstance(polarization_geo, np.ndarray):
        polarization_geo_array = polarization_geo.copy()
    else:
        polarization_geo_array = np.array([polarization_geo])
    num = len(t_geo_array)
    t_ssb_array, longitude_ssb_array = np.zeros(num), np.zeros(num)
    latitude_ssb_array, polarization_ssb_array = np.zeros(num), np.zeros(num)

    for i in range(num):
        t_geo = t_geo_array[i]
        longitude_geo = longitude_geo_array[i]
        latitude_geo = latitude_geo_array[i]
        polarization_geo = polarization_geo_array[i]
        if longitude_geo < 0 or longitude_geo >= 2*np.pi:
            raise ValueError("Longitude should within [0, 2*pi).")
        if latitude_geo < -np.pi/2 or latitude_geo > np.pi/2:
            raise ValueError("Latitude should within [-pi/2, pi/2].")
        if polarization_geo < 0 or polarization_geo >= 2*np.pi:
            raise ValueError("Polarization angle should within [0, 2*pi).")
        k_geo = localization_to_propagation_vector(longitude_geo, latitude_geo)
        rotation_matrix_geo = rotation_matrix_ssb_to_geo()
        k_ssb = rotation_matrix_geo @ k_geo
        longitude_ssb, latitude_ssb = propagation_vector_to_localization(k_ssb)
        polarization_ssb = polarization_newframe(polarization_geo, k_geo, rotation_matrix_geo.T)
        # As mentioned in LDC manual, the p,q vectors are opposite between 
        # LDC and LAL conventions, see Sec 4.1.5 in <LISA-LCST-SGS-MAN-001>.
        polarization_ssb = np.mod(polarization_ssb-np.pi, 2*np.pi)
        if use_astropy:
            # BarycentricMeanEcliptic doesn't have obstime attribute,
            # it's a good inertial frame, but PrecessedGeocentric is not.
            bme_coord = PrecessedGeocentric(ra=longitude_geo*units.radian,
                                            dec=latitude_geo*units.radian,
                                            equinox='J2000',
                                            obstime=Time(t_geo, format='gps'))
            ssb_sky = bme_coord.transform_to(
                        BarycentricMeanEcliptic(equinox='J2000'))
            longitude_ssb = ssb_sky.lon.rad
            latitude_ssb = ssb_sky.lat.rad
        t_ssb = t_ssb_from_t_geo(t_geo, longitude_ssb, latitude_ssb)
        t_ssb_array[i] = t_ssb
        longitude_ssb_array[i] = longitude_ssb
        latitude_ssb_array[i] = latitude_ssb
        polarization_ssb_array[i] = polarization_ssb
    if num == 1:
        params_ssb = (t_ssb_array[0], longitude_ssb_array[0],
                      latitude_ssb_array[0], polarization_ssb_array[0])
    else:
        params_ssb = (t_ssb_array, longitude_ssb_array,
                      latitude_ssb_array, polarization_ssb_array)

    return params_ssb


def lisa_to_geo(t_lisa, longitude_lisa, latitude_lisa, polarization_lisa,
                t0=TIME_OFFSET_20_DEGREES, use_astropy=True):
    """ Converting the arrive time, the sky localization, and the polarization
    from the LISA frame to the geocentric frame.

    Parameters
    ----------
    t_lisa : float or numpy.array
        The time when a GW signal arrives at the origin of LISA frame.
        In the unit of 's'.
    longitude_lisa : float or numpy.array
        The longitude of a GW signal in LISA frame, in the unit of 'radian'.
    latitude_lisa : float or numpy.array
        The latitude of a GW signal in LISA frame, in the unit of 'radian'.
    polarization_lisa : float or numpy.array
        The polarization angle of a GW signal in LISA frame.
        In the unit of 'radian'.
    t0 : float
        The initial time offset of LISA, in the unit of 's',
        default is 7365189.431698299. This makes sure LISA is behind
        the Earth by 19-23 degrees.
    use_astropy : bool
        Using Astropy to calculate the sky localization or not.
        Default is True.

    Returns
    -------
    (t_geo, longitude_geo, latitude_geo, polarization_geo) : tuple
    t_geo : float or numpy.array
        The time when a GW signal arrives at the origin of geocentric frame.
        In the unit of 's'.
    longitude_geo : float or numpy.array
        The ecliptic longitude of a GW signal in geocentric frame.
        In the unit of 'radian'.
    latitude_geo : float or numpy.array
        The ecliptic latitude of a GW signal in geocentric frame.
        In the unit of 'radian'.
    polarization_geo : float or numpy.array
        The polarization angle of a GW signal in geocentric frame.
        In the unit of 'radian'.
    """
    t_ssb, longitude_ssb, latitude_ssb, polarization_ssb = lisa_to_ssb(
        t_lisa, longitude_lisa, latitude_lisa, polarization_lisa, t0)
    t_geo, longitude_geo, latitude_geo, polarization_geo = ssb_to_geo(
        t_ssb, longitude_ssb, latitude_ssb, polarization_ssb, use_astropy)

    return (t_geo, longitude_geo, latitude_geo, polarization_geo)


def geo_to_lisa(t_geo, longitude_geo, latitude_geo, polarization_geo,
                t0=TIME_OFFSET_20_DEGREES, use_astropy=True):
    """ Converting the arrive time, the sky localization, and the polarization
    from the geocentric frame to the LISA frame.

    Parameters
    ----------
    t_geo : float or numpy.array
        The time when a GW signal arrives at the origin of geocentric frame.
        In the unit of 's'.
    longitude_geo : float or numpy.array
        The longitude of a GW signal in geocentric frame.
        In the unit of 'radian'.
    latitude_geo : float or numpy.array
        The latitude of a GW signal in geocentric frame.
        In the unit of 'radian'.
    polarization_geo : float or numpy.array
        The polarization angle of a GW signal in geocentric frame.
        In the unit of 'radian'.
    t0 : float
        The initial time offset of LISA, in the unit of 's',
        default is 7365189.431698299. This makes sure LISA is behind
        the Earth by 19-23 degrees.
    use_astropy : bool
        Using Astropy to calculate the sky localization or not.
        Default is True.

    Returns
    -------
    (t_lisa, longitude_lisa, latitude_lisa, polarization_lisa) : tuple
    t_lisa : float or numpy.array
        The time when a GW signal arrives at the origin of LISA frame.
        In the unit of 's'.
    longitude_lisa : float or numpy.array
        The longitude of a GW signal in LISA frame, in the unit of 'radian'.
    latitude_lisa : float or numpy.array
        The latitude of a GW signal in LISA frame, in the unit of 'radian'.
    polarization_geo : float or numpy.array
        The polarization angle of a GW signal in LISA frame.
        In the unit of 'radian'.
    """
    t_ssb, longitude_ssb, latitude_ssb, polarization_ssb = geo_to_ssb(
        t_geo, longitude_geo, latitude_geo, polarization_geo, use_astropy)
    t_lisa, longitude_lisa, latitude_lisa, polarization_lisa = ssb_to_lisa(
        t_ssb, longitude_ssb, latitude_ssb, polarization_ssb, t0)

    return (t_lisa, longitude_lisa, latitude_lisa, polarization_lisa)
