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
and ground-based detectors. Note that current LISA orbit used in this module
is a circular orbit, need to be replaced by a more realistic and general orbit
model in the near future.
"""

import logging
import numpy as np

from scipy.spatial.transform import Rotation
from scipy.optimize import fsolve
from astropy import units
from astropy.constants import c, au
from astropy.time import Time
from astropy.coordinates import BarycentricMeanEcliptic, PrecessedGeocentric
from astropy.coordinates import get_body_barycentric
from astropy.coordinates import SkyCoord
from astropy.coordinates.builtin_frames import ecliptic_transforms

logger = logging.getLogger('pycbc.coordinates.space')

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


def localization_to_propagation_vector(longitude, latitude,
                                       use_astropy=True, frame=None):
    """ Converting the sky localization to the corresponding
    propagation unit vector of a GW signal.

    Parameters
    ----------
    longitude : float
        The longitude, in the unit of 'radian'.
    latitude : float
        The latitude, in the unit of 'radian'.
    use_astropy : bool
        Using Astropy to calculate the sky localization or not.
        Default is True.
    frame : astropy.coordinates
        The frame from astropy.coordinates if use_astropy is True,
        the default is None.

    Returns
    -------
    [[x], [y], [z]] : numpy.array
        The propagation unit vector of that GW signal.
    """
    if use_astropy:
        x = -frame.cartesian.x.value
        y = -frame.cartesian.y.value
        z = -frame.cartesian.z.value
    else:
        x = -np.cos(latitude) * np.cos(longitude)
        y = -np.cos(latitude) * np.sin(longitude)
        z = -np.sin(latitude)
    v = np.array([[x], [y], [z]])

    return v / np.linalg.norm(v)


def propagation_vector_to_localization(k, use_astropy=True, frame=None):
    """ Converting the propagation unit vector to the corresponding
    sky localization of a GW signal.

    Parameters
    ----------
    k : numpy.array
        The propagation unit vector of a GW signal.
    use_astropy : bool
        Using Astropy to calculate the sky localization or not.
        Default is True.
    frame : astropy.coordinates
        The frame from astropy.coordinates if use_astropy is True,
        the default is None.

    Returns
    -------
    (longitude, latitude) : tuple
        The sky localization of that GW signal.
    """
    if use_astropy:
        try:
            longitude = frame.lon.rad
            latitude = frame.lat.rad
        except AttributeError:
            longitude = frame.ra.rad
            latitude = frame.dec.rad
    else:
        # latitude already within [-pi/2, pi/2]
        latitude = np.float64(np.arcsin(-k[2,0]))
        longitude = np.float64(np.arctan2(-k[1,0]/np.cos(latitude),
                               -k[0,0]/np.cos(latitude)))
        # longitude should within [0, 2*pi)
        longitude = np.mod(longitude, 2*np.pi)

    return (longitude, latitude)


def polarization_newframe(polarization, k, rotation_matrix, use_astropy=True,
                          old_frame=None, new_frame=None):
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
    use_astropy : bool
        Using Astropy to calculate the sky localization or not.
        Default is True.
    old_frame : astropy.coordinates
        The frame from astropy.coordinates if use_astropy is True,
        the default is None.
    new_frame : astropy.coordinates
        The frame from astropy.coordinates if use_astropy is True,
        the default is None. The new frame for the new polarization
        angle.

    Returns
    -------
    polarization_new_frame : float
        The polarization angle in the new frame of that GW signal.
    """
    longitude, _ = propagation_vector_to_localization(
                        k, use_astropy, old_frame)
    u = np.array([[np.sin(longitude)], [-np.cos(longitude)], [0]])
    rotation_vector = polarization * k
    rotation_polarization = Rotation.from_rotvec(rotation_vector.T[0])
    p = rotation_polarization.apply(u.T[0]).reshape(3, 1)
    p_newframe = rotation_matrix.T @ p
    k_newframe = rotation_matrix.T @ k
    longitude_newframe, latitude_newframe = \
        propagation_vector_to_localization(k_newframe, use_astropy, new_frame)
    u_newframe = np.array([[np.sin(longitude_newframe)],
                           [-np.cos(longitude_newframe)], [0]])
    v_newframe = np.array([
                    [-np.sin(latitude_newframe) * np.cos(longitude_newframe)],
                    [-np.sin(latitude_newframe) * np.sin(longitude_newframe)],
                    [np.cos(latitude_newframe)]])
    p_dot_u_newframe = np.vdot(p_newframe, u_newframe)
    p_dot_v_newframe = np.vdot(p_newframe, v_newframe)
    polarization_new_frame = np.arctan2(p_dot_v_newframe, p_dot_u_newframe)
    polarization_new_frame = np.mod(polarization_new_frame, 2*np.pi)
    # avoid the round error
    if polarization_new_frame == 2*np.pi:
        polarization_new_frame = 0

    return polarization_new_frame


def t_lisa_from_ssb(t_ssb, longitude_ssb, latitude_ssb,
                    t0=TIME_OFFSET_20_DEGREES):
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
    k = localization_to_propagation_vector(
            longitude_ssb, latitude_ssb, use_astropy=False)

    def equation(t_lisa):
        # LISA is moving, when GW arrives at LISA center,
        # time is t_lisa, not t_ssb.
        p = lisa_position_ssb(t_lisa, t0)[0]
        return t_lisa - t_ssb - np.vdot(k, p) / c.value

    return fsolve(equation, t_ssb)[0]


def t_ssb_from_t_lisa(t_lisa, longitude_ssb, latitude_ssb,
                      t0=TIME_OFFSET_20_DEGREES):
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
    k = localization_to_propagation_vector(
            longitude_ssb, latitude_ssb, use_astropy=False)
    # LISA is moving, when GW arrives at LISA center,
    # time is t_lisa, not t_ssb.
    p = lisa_position_ssb(t_lisa, t0)[0]

    def equation(t_ssb):
        return t_lisa - t_ssb - np.vdot(k, p) / c.value

    return fsolve(equation, t_lisa)[0]


def _ensure_sky_params_arrays(*arrays):
    """ Shared by the 4 X_to_ssb/ssb_to_X frame-transform functions below:
    coerce each of `arrays` (scalar or array-like) to a 1-D numpy array,
    and return them alongside `num`, the common length (from the first
    array).
    """
    coerced = tuple(
        a if isinstance(a, np.ndarray) else np.array([a]) for a in arrays)
    return coerced + (len(coerced[0]),)


def _validate_sky_params(longitude, latitude, polarization):
    """ Shared by the 4 X_to_ssb/ssb_to_X frame-transform functions below:
    range-checks the whole array at once (equivalent to the same check
    repeated per element inside their loops), raising the same
    ValueError text.
    """
    if np.any((longitude < 0) | (longitude >= 2*np.pi)):
        raise ValueError("Longitude should within [0, 2*pi).")
    if np.any((latitude < -np.pi/2) | (latitude > np.pi/2)):
        raise ValueError("Latitude should within [-pi/2, pi/2].")
    if np.any((polarization < 0) | (polarization >= 2*np.pi)):
        raise ValueError("Polarization angle should within [0, 2*pi).")


def _pack_sky_params_output(num, t, longitude, latitude, polarization):
    """ Shared by the 4 X_to_ssb/ssb_to_X frame-transform functions below:
    unwrap to scalars for a single input, otherwise keep the arrays.
    """
    if num == 1:
        return (t[0], longitude[0], latitude[0], polarization[0])
    return (t, longitude, latitude, polarization)


def _sky_ssb_space_transform(t_src, longitude_src, latitude_src,
                             polarization_src, t0, forward):
    """ Shared body of `ssb_to_lisa` and `lisa_to_ssb`: rotate the GW
    propagation vector between the SSB frame and the LISA frame at the
    appropriate LISA-frame time, and convert the arrival time to match.
    `forward=True` runs the SSB->LISA direction (`ssb_to_lisa`);
    `forward=False` runs LISA->SSB (`lisa_to_ssb`).
    """
    t_src, longitude_src, latitude_src, polarization_src, num = \
        _ensure_sky_params_arrays(
            t_src, longitude_src, latitude_src, polarization_src)
    _validate_sky_params(longitude_src, latitude_src, polarization_src)
    t_dst = np.zeros(num)
    longitude_dst, latitude_dst = np.zeros(num), np.zeros(num)
    polarization_dst = np.zeros(num)

    for i in range(num):
        k_src = localization_to_propagation_vector(
            longitude_src[i], latitude_src[i], use_astropy=False)
        if forward:
            t_dst[i] = t_lisa_from_ssb(
                t_src[i], longitude_src[i], latitude_src[i], t0)
            # t_dst was computed with the t0-corrected LISA position but
            # corresponds to the true t_ssb, not t_ssb + t0, so t0 has to
            # be applied again here to place LISA correctly.
            rot = rotation_matrix_ssb_to_lisa(
                lisa_position_ssb(t_dst[i], t0)[1])
            k_dst = rot.T @ k_src
        else:
            rot = rotation_matrix_ssb_to_lisa(
                lisa_position_ssb(t_src[i], t0)[1])
            k_dst = rot @ k_src
        longitude_dst[i], latitude_dst[i] = \
            propagation_vector_to_localization(k_dst, use_astropy=False)
        if not forward:
            t_dst[i] = t_ssb_from_t_lisa(
                t_src[i], longitude_dst[i], latitude_dst[i], t0)
        polarization_dst[i] = polarization_newframe(
            polarization_src[i], k_src, rot if forward else rot.T,
            use_astropy=False)

    return _pack_sky_params_output(
        num, t_dst, longitude_dst, latitude_dst, polarization_dst)


def ssb_to_lisa(t_ssb, longitude_ssb, latitude_ssb, polarization_ssb,
                t0=TIME_OFFSET_20_DEGREES):
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
    return _sky_ssb_space_transform(
        t_ssb, longitude_ssb, latitude_ssb, polarization_ssb, t0,
        forward=True)


def lisa_to_ssb(t_lisa, longitude_lisa, latitude_lisa, polarization_lisa,
                t0=TIME_OFFSET_20_DEGREES):
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
    return _sky_ssb_space_transform(
        t_lisa, longitude_lisa, latitude_lisa, polarization_lisa, t0,
        forward=False)




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


def t_geo_from_ssb(t_ssb, longitude_ssb, latitude_ssb,
                   use_astropy=True, frame=None):
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
    k = localization_to_propagation_vector(
            longitude_ssb, latitude_ssb, use_astropy, frame)

    def equation(t_geo):
        # Earth is moving, when GW arrives at Earth center,
        # time is t_geo, not t_ssb.
        p = earth_position_ssb(t_geo)[0]
        return t_geo - t_ssb - np.vdot(k, p) / c.value

    return fsolve(equation, t_ssb)[0]


def t_ssb_from_t_geo(t_geo, longitude_ssb, latitude_ssb,
                     use_astropy=True, frame=None):
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
    k = localization_to_propagation_vector(
            longitude_ssb, latitude_ssb, use_astropy, frame)
    # Earth is moving, when GW arrives at Earth center,
    # time is t_geo, not t_ssb.
    p = earth_position_ssb(t_geo)[0]

    def equation(t_ssb):
        return t_geo - t_ssb - np.vdot(k, p) / c.value

    return fsolve(equation, t_geo)[0]


def ssb_to_geo(t_ssb, longitude_ssb, latitude_ssb, polarization_ssb,
               use_astropy=True):
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
    t_ssb, longitude_ssb, latitude_ssb, polarization_ssb, num = \
        _ensure_sky_params_arrays(
            t_ssb, longitude_ssb, latitude_ssb, polarization_ssb)
    _validate_sky_params(longitude_ssb, latitude_ssb, polarization_ssb)
    t_geo = np.full(num, np.nan)
    longitude_geo = np.full(num, np.nan)
    latitude_geo = np.full(num, np.nan)
    polarization_geo = np.full(num, np.nan)

    for i in range(num):
        if use_astropy:
            # BarycentricMeanEcliptic doesn't have obstime attribute,
            # it's a good inertial frame, but PrecessedGeocentric is not.
            bme_coord = BarycentricMeanEcliptic(
                            lon=longitude_ssb[i]*units.radian,
                            lat=latitude_ssb[i]*units.radian,
                            equinox='J2000')
            t_geo[i] = t_geo_from_ssb(t_ssb[i], longitude_ssb[i],
                                      latitude_ssb[i], use_astropy, bme_coord)
            geo_sky = bme_coord.transform_to(PrecessedGeocentric(
                equinox='J2000', obstime=Time(t_geo[i], format='gps')))
            longitude_geo[i] = geo_sky.ra.rad
            latitude_geo[i] = geo_sky.dec.rad
            k_geo = localization_to_propagation_vector(
                        longitude_geo[i], latitude_geo[i],
                        use_astropy, geo_sky)
            k_ssb = localization_to_propagation_vector(
                        None, None, use_astropy, bme_coord)
            rotation_matrix_geo = \
                ecliptic_transforms.icrs_to_baryecliptic(
                    from_coo=None,
                    to_frame=BarycentricMeanEcliptic(equinox='J2000'))
            polarization_geo[i] = polarization_newframe(
                                    polarization_ssb[i], k_ssb,
                                    rotation_matrix_geo, use_astropy,
                                    old_frame=bme_coord,
                                    new_frame=geo_sky)
        else:
            t_geo[i] = t_geo_from_ssb(t_ssb[i], longitude_ssb[i],
                                      latitude_ssb[i], use_astropy)
            rotation_matrix_geo = rotation_matrix_ssb_to_geo()
            k_ssb = localization_to_propagation_vector(
                        longitude_ssb[i], latitude_ssb[i],
                        use_astropy)
            k_geo = rotation_matrix_geo.T @ k_ssb
            longitude_geo[i], latitude_geo[i] = \
                propagation_vector_to_localization(k_geo, use_astropy)
            polarization_geo[i] = polarization_newframe(
                                    polarization_ssb[i], k_ssb,
                                    rotation_matrix_geo, use_astropy)

        # As mentioned in LDC manual, the p,q vectors are opposite between
        # LDC and LAL conventions, see Sec 4.1.5 in <LISA-LCST-SGS-MAN-001>.
        polarization_geo[i] = np.mod(polarization_geo[i]+np.pi, 2*np.pi)

    return _pack_sky_params_output(
        num, t_geo, longitude_geo, latitude_geo, polarization_geo)


def geo_to_ssb(t_geo, longitude_geo, latitude_geo, polarization_geo,
               use_astropy=True):
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
    t_geo, longitude_geo, latitude_geo, polarization_geo, num = \
        _ensure_sky_params_arrays(
            t_geo, longitude_geo, latitude_geo, polarization_geo)
    _validate_sky_params(longitude_geo, latitude_geo, polarization_geo)
    t_ssb = np.full(num, np.nan)
    longitude_ssb = np.full(num, np.nan)
    latitude_ssb = np.full(num, np.nan)
    polarization_ssb = np.full(num, np.nan)

    for i in range(num):
        if use_astropy:
            # BarycentricMeanEcliptic doesn't have obstime attribute,
            # it's a good inertial frame, but PrecessedGeocentric is not.
            geo_coord = PrecessedGeocentric(
                            ra=longitude_geo[i]*units.radian,
                            dec=latitude_geo[i]*units.radian,
                            equinox='J2000',
                            obstime=Time(t_geo[i], format='gps'))
            ssb_sky = geo_coord.transform_to(
                        BarycentricMeanEcliptic(equinox='J2000'))
            longitude_ssb[i] = ssb_sky.lon.rad
            latitude_ssb[i] = ssb_sky.lat.rad
            k_ssb = localization_to_propagation_vector(
                        longitude_ssb[i], latitude_ssb[i],
                        use_astropy, ssb_sky)
            k_geo = localization_to_propagation_vector(
                None, None, use_astropy, geo_coord)
            rotation_matrix_geo = \
                ecliptic_transforms.icrs_to_baryecliptic(
                    from_coo=None,
                    to_frame=BarycentricMeanEcliptic(equinox='J2000'))
            t_ssb[i] = t_ssb_from_t_geo(t_geo[i], longitude_ssb[i],
                                        latitude_ssb[i], use_astropy,
                                        ssb_sky)
            polarization_ssb[i] = polarization_newframe(
                                    polarization_geo[i], k_geo,
                                    rotation_matrix_geo.T,
                                    use_astropy,
                                    old_frame=geo_coord,
                                    new_frame=ssb_sky)
        else:
            rotation_matrix_geo = rotation_matrix_ssb_to_geo()
            k_geo = localization_to_propagation_vector(
                        longitude_geo[i], latitude_geo[i], use_astropy)
            k_ssb = rotation_matrix_geo @ k_geo
            longitude_ssb[i], latitude_ssb[i] = \
                propagation_vector_to_localization(k_ssb, use_astropy)
            t_ssb[i] = t_ssb_from_t_geo(t_geo[i], longitude_ssb[i],
                                        latitude_ssb[i], use_astropy)
            polarization_ssb[i] = polarization_newframe(
                                    polarization_geo[i], k_geo,
                                    rotation_matrix_geo.T, use_astropy)

        # As mentioned in LDC manual, the p,q vectors are opposite between
        # LDC and LAL conventions, see Sec 4.1.5 in <LISA-LCST-SGS-MAN-001>.
        polarization_ssb[i] = np.mod(polarization_ssb[i]-np.pi, 2*np.pi)

    return _pack_sky_params_output(
        num, t_ssb, longitude_ssb, latitude_ssb, polarization_ssb)



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


__all__ = ['TIME_OFFSET_20_DEGREES',
           'localization_to_propagation_vector',
           'propagation_vector_to_localization', 'polarization_newframe',
           't_lisa_from_ssb', 't_ssb_from_t_lisa',
           'ssb_to_lisa', 'lisa_to_ssb',
           'rotation_matrix_ssb_to_lisa', 'rotation_matrix_ssb_to_geo',
           'lisa_position_ssb', 'earth_position_ssb',
           't_geo_from_ssb', 't_ssb_from_t_geo', 'ssb_to_geo', 'geo_to_ssb',
           'lisa_to_geo', 'geo_to_lisa',
           ]
