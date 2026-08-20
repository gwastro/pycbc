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
and ground-based detectors. The default orbit assumed throughout is LISA's
analytic circular orbit; `t_space_from_ssb`, `t_ssb_from_t_space`,
`ssb_to_space`, `space_to_ssb`, `space_to_geo` and `geo_to_space` also accept
an optional `orbit` argument (see `pycbc.coordinates.space_orbit`) to use an
arbitrary numerical or analytic constellation orbit (LISA, Taiji, TianQin,
...) instead, with no change in behavior when `orbit` is left at its default
of None.

The `*_lisa_*`/`*_to_lisa`/`lisa_to_*`-named versions of these functions
(`t_lisa_from_ssb`, `ssb_to_lisa`, ...) are deprecated aliases kept only for
backward compatibility (e.g. the BBHx waveform plugin imports and calls
several of these by name, with these exact keyword arguments) -- new code
should use the `*_space_*`/`*_to_space`/`space_to_*` names instead, which
mean exactly the same thing: "the constellation frame described by `orbit`",
not literally LISA. Both spellings always behave identically; the old names
are not scheduled for removal.
"""

import logging

import numpy as np
from astropy import units
from astropy.constants import au
from astropy.coordinates import (
    BarycentricMeanEcliptic,
    PrecessedGeocentric,
    SkyCoord,
    get_body_barycentric,
)
from astropy.coordinates.builtin_frames import ecliptic_transforms
from astropy.time import Time
from scipy.spatial.transform import Rotation

from pycbc.coordinates.space_orbit import (
    _solve_frame_arrival_time,
    constellation_frame,
    t_detector_from_ssb,
    t_ssb_from_t_detector,
)

logger = logging.getLogger('pycbc.coordinates.space')

# This constant makes sure LISA is behind the Earth by 19-23 degrees.
# Making this a stand-alone constant will also make it callable by
# the waveform plugin and PE config file. In the unit of 's'.
TIME_OFFSET_20_DEGREES = 7365189.431698299

# "rotation_matrix_ssb_to_space" and "space_position_ssb" remain specific to
# LISA's analytic circular orbit; they are the default (`orbit=None`) code
# path used by the functions below. For any other constellation orbit, pass
# an `orbit` argument to those functions instead (see
# `pycbc.coordinates.space_orbit`).


def rotation_matrix_ssb_to_space(alpha):
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


def rotation_matrix_ssb_to_lisa(alpha):
    """Deprecated alias for `rotation_matrix_ssb_to_space` (see module
    docstring)."""
    return rotation_matrix_ssb_to_space(alpha)


def space_position_ssb(t_space, t0=TIME_OFFSET_20_DEGREES):
    """ LISA's position vector and angular displacement in the SSB frame at
    a given time, assuming a circular orbit in the ecliptic, one year
    period, trailing the Earth.

    Parameters
    ----------
    t_space : float
        Any time [s]; not necessarily an actual GW arrival time.
    t0 : float, optional
        LISA's initial time offset [s], default 7365189.431698299 (keeps
        LISA 19-23 degrees behind Earth).

    Returns
    -------
    (p, alpha) : tuple
        LISA's position vector [m] and angular displacement [rad] in the
        SSB frame.
    """
    OMEGA_0 = 1.99098659277e-7  # 2*pi / sidereal year [rad/s]
    R_ORBIT = au.value
    alpha = np.mod(OMEGA_0 * (t_space + t0), 2*np.pi)
    p = np.array([[R_ORBIT * np.cos(alpha)],
                  [R_ORBIT * np.sin(alpha)],
                  [0]], dtype=object)
    return (p, alpha)


def lisa_position_ssb(t_lisa, t0=TIME_OFFSET_20_DEGREES):
    """Deprecated alias for `space_position_ssb` (see module docstring)."""
    return space_position_ssb(t_lisa, t0)


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


def _rotation_matrix_at_detector_time(t_detector, t0, orbit, sc=(1, 2, 3)):
    """ Internal helper shared by `ssb_to_space` and `space_to_ssb`: the
    rotation matrix (of frame basis) from the SSB frame to the detector
    frame at a given detector-frame time, from either the analytic circular
    LISA orbit (`orbit=None`) or a general orbit provider (see
    `pycbc.coordinates.space_orbit`).
    """
    if orbit is not None:
        _, rotation = constellation_frame([t_detector], orbit, sc=sc)
        return rotation[0]
    alpha = space_position_ssb(t_detector, t0)[1]
    return rotation_matrix_ssb_to_space(alpha)


def _sky_ssb_space_transform(t_src, longitude_src, latitude_src,
                              polarization_src, t0, orbit, sc, forward):
    """ Shared body of `ssb_to_space` and `space_to_ssb`: rotate the GW
    propagation vector between the SSB frame and the constellation frame
    at the appropriate detector-frame time, and convert the arrival time
    to match. `forward=True` runs the SSB->space direction (`ssb_to_space`);
    `forward=False` runs space->SSB (`space_to_ssb`).
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
            t_dst[i] = t_space_from_ssb(
                t_src[i], longitude_src[i], latitude_src[i], t0,
                orbit=orbit, sc=sc)
            # t0 must be reapplied here: t_dst already used it once, but
            # is LISA's arrival time for the true t_ssb, not t_ssb + t0.
            rot = _rotation_matrix_at_detector_time(t_dst[i], t0, orbit, sc)
            k_dst = rot.T @ k_src
        else:
            rot = _rotation_matrix_at_detector_time(t_src[i], t0, orbit, sc)
            k_dst = rot @ k_src
        longitude_dst[i], latitude_dst[i] = \
            propagation_vector_to_localization(k_dst, use_astropy=False)
        if not forward:
            t_dst[i] = t_ssb_from_t_space(
                t_src[i], longitude_dst[i], latitude_dst[i], t0,
                orbit=orbit, sc=sc)
        polarization_dst[i] = polarization_newframe(
            polarization_src[i], k_src, rot if forward else rot.T,
            use_astropy=False)

    return _pack_sky_params_output(
        num, t_dst, longitude_dst, latitude_dst, polarization_dst)


def t_space_from_ssb(t_ssb, longitude_ssb, latitude_ssb,
                     t0=TIME_OFFSET_20_DEGREES, orbit=None, sc=(1, 2, 3)):
    """ Arrival time at LISA's barycenter, from arrival time and sky
    localization in the SSB frame.

    Parameters
    ----------
    t_ssb : float
        Arrival time at the SSB frame origin [s].
    longitude_ssb : float
        Ecliptic longitude in the SSB frame [rad].
    latitude_ssb : float
        Ecliptic latitude in the SSB frame [rad].
    t0 : float, optional
        LISA's initial time offset [s], default 7365189.431698299 (keeps
        LISA 19-23 degrees behind Earth). Ignored if `orbit` is given.
    orbit : OrbitProvider, optional
        An object exposing `compute_position(t, sc)` (see
        `pycbc.coordinates.space_orbit`), giving the true constellation
        orbit (LISA, Taiji, TianQin, numerical, ...) to use instead of the
        analytic circular LISA orbit. Default None reproduces the
        behavior of previous versions of this function exactly.
    sc : tuple, optional
        1-indexed spacecraft labels defining the constellation. Only used
        if `orbit` is given. Default (1, 2, 3).

    Returns
    -------
    t_space : float
        Arrival time at the LISA frame origin [s].
    """
    k = localization_to_propagation_vector(
            longitude_ssb, latitude_ssb, use_astropy=False)

    if orbit is not None:
        return t_detector_from_ssb(t_ssb, k.flatten(), orbit, sc=sc)

    # LISA is moving, when GW arrives at LISA center,
    # time is t_space, not t_ssb.
    def position_fn(t):
        return space_position_ssb(t, t0)[0]
    return _solve_frame_arrival_time(t_ssb, k, position_fn, forward=True)


def t_lisa_from_ssb(t_ssb, longitude_ssb, latitude_ssb,
                    t0=TIME_OFFSET_20_DEGREES, orbit=None, sc=(1, 2, 3)):
    """Deprecated alias for `t_space_from_ssb` (see module docstring)."""
    return t_space_from_ssb(t_ssb, longitude_ssb, latitude_ssb, t0, orbit, sc)


def t_ssb_from_t_space(t_space, longitude_ssb, latitude_ssb,
                       t0=TIME_OFFSET_20_DEGREES, orbit=None, sc=(1, 2, 3)):
    """ Arrival time at the SSB frame origin, from arrival time in the LISA
    frame and sky localization in the SSB frame.

    Parameters
    ----------
    t_space : float
        Arrival time at the LISA frame origin [s].
    longitude_ssb : float
        Ecliptic longitude in the SSB frame [rad].
    latitude_ssb : float
        Ecliptic latitude in the SSB frame [rad].
    t0 : float, optional
        LISA's initial time offset [s], default 7365189.431698299 (keeps
        LISA 19-23 degrees behind Earth). Ignored if `orbit` is given.
    orbit : OrbitProvider, optional
        See `t_space_from_ssb`. Default None reproduces the behavior of
        previous versions of this function exactly.
    sc : tuple, optional
        1-indexed spacecraft labels defining the constellation. Only used
        if `orbit` is given. Default (1, 2, 3).

    Returns
    -------
    t_ssb : float
        Arrival time at the SSB frame origin [s].
    """
    k = localization_to_propagation_vector(
            longitude_ssb, latitude_ssb, use_astropy=False)

    if orbit is not None:
        return t_ssb_from_t_detector(t_space, k.flatten(), orbit, sc=sc)

    # LISA is moving, when GW arrives at LISA center,
    # time is t_space, not t_ssb.
    def position_fn(t):
        return space_position_ssb(t, t0)[0]
    return _solve_frame_arrival_time(t_space, k, position_fn, forward=False)


def t_ssb_from_t_lisa(t_lisa, longitude_ssb, latitude_ssb,
                      t0=TIME_OFFSET_20_DEGREES, orbit=None, sc=(1, 2, 3)):
    """Deprecated alias for `t_ssb_from_t_space` (see module docstring)."""
    return t_ssb_from_t_space(t_lisa, longitude_ssb, latitude_ssb, t0, orbit,
                              sc)


def ssb_to_space(t_ssb, longitude_ssb, latitude_ssb, polarization_ssb,
                 t0=TIME_OFFSET_20_DEGREES, orbit=None, sc=(1, 2, 3)):
    """ Converts arrival time, sky localization, and polarization from the
    SSB frame to the LISA frame.

    Parameters
    ----------
    t_ssb : float or numpy.array
        Arrival time at the SSB frame origin [s].
    longitude_ssb : float or numpy.array
        Ecliptic longitude in the SSB frame [rad].
    latitude_ssb : float or numpy.array
        Ecliptic latitude in the SSB frame [rad].
    polarization_ssb : float or numpy.array
        Polarization angle in the SSB frame [rad].
    t0 : float, optional
        LISA's initial time offset [s], default 7365189.431698299 (keeps
        LISA 19-23 degrees behind Earth). Ignored if `orbit` is given.
    orbit : OrbitProvider, optional
        An object exposing `compute_position(t, sc)` (see
        `pycbc.coordinates.space_orbit`), giving the true constellation
        orbit (LISA, Taiji, TianQin, numerical, ...) to use instead of the
        analytic circular LISA orbit. Default None reproduces the
        behavior of previous versions of this function exactly.
    sc : tuple, optional
        1-indexed spacecraft labels defining the constellation. Only used
        if `orbit` is given. Default (1, 2, 3).

    Returns
    -------
    (t_space, longitude_space, latitude_space, polarization_space) : tuple
        Arrival time [s], ecliptic longitude [rad], ecliptic latitude
        [rad], and polarization angle [rad] in the LISA frame.
    """
    return _sky_ssb_space_transform(
        t_ssb, longitude_ssb, latitude_ssb, polarization_ssb,
        t0, orbit, sc, forward=True)


def ssb_to_lisa(t_ssb, longitude_ssb, latitude_ssb, polarization_ssb,
                t0=TIME_OFFSET_20_DEGREES, orbit=None, sc=(1, 2, 3)):
    """Deprecated alias for `ssb_to_space` (see module docstring)."""
    return ssb_to_space(t_ssb, longitude_ssb, latitude_ssb, polarization_ssb,
                        t0, orbit, sc)


def space_to_ssb(t_space, longitude_space, latitude_space,
                 polarization_space, t0=TIME_OFFSET_20_DEGREES, orbit=None,
                 sc=(1, 2, 3)):
    """ Converts arrival time, sky localization, and polarization from the
    LISA frame to the SSB frame.

    Parameters
    ----------
    t_space : float or numpy.array
        Arrival time at the LISA frame origin [s].
    longitude_space : float or numpy.array
        Longitude in the LISA frame [rad].
    latitude_space : float or numpy.array
        Latitude in the LISA frame [rad].
    polarization_space : float or numpy.array
        Polarization angle in the LISA frame [rad].
    t0 : float, optional
        LISA's initial time offset [s], default 7365189.431698299 (keeps
        LISA 19-23 degrees behind Earth). Ignored if `orbit` is given.
    orbit : OrbitProvider, optional
        See `ssb_to_space`. Default None reproduces the behavior of
        previous versions of this function exactly.
    sc : tuple, optional
        1-indexed spacecraft labels defining the constellation. Only used
        if `orbit` is given. Default (1, 2, 3).

    Returns
    -------
    (t_ssb, longitude_ssb, latitude_ssb, polarization_ssb) : tuple
        Arrival time [s], ecliptic longitude [rad], ecliptic latitude
        [rad], and polarization angle [rad] in the SSB frame.
    """
    return _sky_ssb_space_transform(
        t_space, longitude_space, latitude_space, polarization_space,
        t0, orbit, sc, forward=False)


def lisa_to_ssb(t_lisa, longitude_lisa, latitude_lisa, polarization_lisa,
                t0=TIME_OFFSET_20_DEGREES, orbit=None, sc=(1, 2, 3)):
    """Deprecated alias for `space_to_ssb` (see module docstring) -- kept
    because the BBHx waveform plugin imports and calls this by name with
    these exact keyword arguments."""
    return space_to_ssb(t_lisa, longitude_lisa, latitude_lisa,
                        polarization_lisa, t0, orbit, sc)


def rotation_matrix_ssb_to_geo(epsilon=None):
    """ The rotation matrix (of frame basis) from SSB frame to
    geocentric frame.

    Parameters
    ----------
    epsilon : float or None, optional
        The Earth's axial tilt (obliquity), in the unit of 'radian'.
        Default None, which computes the obliquity directly from astropy
        (via the same cached ICRS <-> ecliptic rotation
        `pycbc.coordinates.space_orbit._icrs_to_ecliptic_rotation_matrix`
        uses) instead of a hard-coded constant; numerically equivalent to
        the IAU obliquity value used there to ~1e-7 rad. Pass an explicit
        value to use a simple analytic rotation by that angle instead
        (e.g. for a non-standard or historical obliquity).

    Returns
    -------
    r : numpy.array
        A 3x3 rotation matrix from SSB frame to geocentric frame.
    """
    if epsilon is None:
        # Imported lazily to avoid a circular import (space_orbit imports
        # from this module too).
        from pycbc.coordinates.space_orbit import _icrs_to_ecliptic_rotation_matrix
        return _icrs_to_ecliptic_rotation_matrix()

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

    # Earth is moving, when GW arrives at Earth center,
    # time is t_geo, not t_ssb.
    def position_fn(t):
        return earth_position_ssb(t)[0]
    return _solve_frame_arrival_time(t_ssb, k, position_fn, forward=True)


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
    def position_fn(t):
        return earth_position_ssb(t)[0]
    return _solve_frame_arrival_time(t_geo, k, position_fn, forward=False)


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


def space_to_geo(t_space, longitude_space, latitude_space,
                 polarization_space, t0=TIME_OFFSET_20_DEGREES,
                 use_astropy=True, orbit=None, sc=(1, 2, 3)):
    """ Converts arrival time, sky localization, and polarization from the
    LISA frame to the geocentric frame.

    Parameters
    ----------
    t_space : float or numpy.array
        Arrival time at the LISA frame origin [s].
    longitude_space : float or numpy.array
        Longitude in the LISA frame [rad].
    latitude_space : float or numpy.array
        Latitude in the LISA frame [rad].
    polarization_space : float or numpy.array
        Polarization angle in the LISA frame [rad].
    t0 : float, optional
        LISA's initial time offset [s], default 7365189.431698299 (keeps
        LISA 19-23 degrees behind Earth). Ignored if `orbit` is given.
    use_astropy : bool, optional
        Whether to use astropy for the sky-localization step. Default True.
    orbit : OrbitProvider, optional
        See `ssb_to_space`. Default None reproduces the behavior of
        previous versions of this function exactly. Only the LISA<->SSB
        leg of this transform depends on the constellation orbit; the
        geocentric frame itself is unaffected.
    sc : tuple, optional
        1-indexed spacecraft labels defining the constellation. Only used
        if `orbit` is given. Default (1, 2, 3).

    Returns
    -------
    (t_geo, longitude_geo, latitude_geo, polarization_geo) : tuple
        Arrival time [s], ecliptic longitude [rad], ecliptic latitude
        [rad], and polarization angle [rad] in the geocentric frame.
    """
    t_ssb, longitude_ssb, latitude_ssb, polarization_ssb = space_to_ssb(
        t_space, longitude_space, latitude_space, polarization_space, t0,
        orbit=orbit, sc=sc)
    t_geo, longitude_geo, latitude_geo, polarization_geo = ssb_to_geo(
        t_ssb, longitude_ssb, latitude_ssb, polarization_ssb, use_astropy)

    return (t_geo, longitude_geo, latitude_geo, polarization_geo)


def lisa_to_geo(t_lisa, longitude_lisa, latitude_lisa, polarization_lisa,
                t0=TIME_OFFSET_20_DEGREES, use_astropy=True, orbit=None,
                sc=(1, 2, 3)):
    """Deprecated alias for `space_to_geo` (see module docstring)."""
    return space_to_geo(t_lisa, longitude_lisa, latitude_lisa,
                        polarization_lisa, t0, use_astropy, orbit, sc)


def geo_to_space(t_geo, longitude_geo, latitude_geo, polarization_geo,
                 t0=TIME_OFFSET_20_DEGREES, use_astropy=True, orbit=None,
                 sc=(1, 2, 3)):
    """ Converts arrival time, sky localization, and polarization from the
    geocentric frame to the LISA frame.

    Parameters
    ----------
    t_geo : float or numpy.array
        Arrival time at the geocentric frame origin [s].
    longitude_geo : float or numpy.array
        Longitude in the geocentric frame [rad].
    latitude_geo : float or numpy.array
        Latitude in the geocentric frame [rad].
    polarization_geo : float or numpy.array
        Polarization angle in the geocentric frame [rad].
    t0 : float, optional
        LISA's initial time offset [s], default 7365189.431698299 (keeps
        LISA 19-23 degrees behind Earth). Ignored if `orbit` is given.
    use_astropy : bool, optional
        Whether to use astropy for the sky-localization step. Default True.
    orbit : OrbitProvider, optional
        See `ssb_to_space`. Default None reproduces the behavior of
        previous versions of this function exactly. Only the LISA<->SSB
        leg of this transform depends on the constellation orbit; the
        geocentric frame itself is unaffected.
    sc : tuple, optional
        1-indexed spacecraft labels defining the constellation. Only used
        if `orbit` is given. Default (1, 2, 3).

    Returns
    -------
    (t_space, longitude_space, latitude_space, polarization_space) : tuple
        Arrival time [s], longitude [rad], latitude [rad], and
        polarization angle [rad] in the LISA frame.
    """
    t_ssb, longitude_ssb, latitude_ssb, polarization_ssb = geo_to_ssb(
        t_geo, longitude_geo, latitude_geo, polarization_geo, use_astropy)
    t_space, longitude_space, latitude_space, polarization_space = \
        ssb_to_space(t_ssb, longitude_ssb, latitude_ssb, polarization_ssb,
                    t0, orbit=orbit, sc=sc)

    return (t_space, longitude_space, latitude_space, polarization_space)


def geo_to_lisa(t_geo, longitude_geo, latitude_geo, polarization_geo,
                t0=TIME_OFFSET_20_DEGREES, use_astropy=True, orbit=None,
                sc=(1, 2, 3)):
    """Deprecated alias for `geo_to_space` (see module docstring)."""
    return geo_to_space(t_geo, longitude_geo, latitude_geo,
                        polarization_geo, t0, use_astropy, orbit, sc)


__all__ = ['TIME_OFFSET_20_DEGREES',
           'localization_to_propagation_vector',
           'propagation_vector_to_localization', 'polarization_newframe',
           't_space_from_ssb', 't_ssb_from_t_space',
           't_lisa_from_ssb', 't_ssb_from_t_lisa',
           'ssb_to_space', 'space_to_ssb',
           'ssb_to_lisa', 'lisa_to_ssb',
           'rotation_matrix_ssb_to_space', 'rotation_matrix_ssb_to_lisa',
           'rotation_matrix_ssb_to_geo',
           'space_position_ssb', 'lisa_position_ssb', 'earth_position_ssb',
           't_geo_from_ssb', 't_ssb_from_t_geo', 'ssb_to_geo', 'geo_to_ssb',
           'space_to_geo', 'geo_to_space',
           'lisa_to_geo', 'geo_to_lisa',
           ]
