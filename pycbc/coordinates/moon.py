# Copyright (C) 2026  Shichao Wu, Alex Nitz
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
This module extends `pycbc.coordinates.space`'s SSB-hub coordinate
transforms to a lunar detector (e.g. LGWA), for arrival time, sky
localization, and polarization angle -- the pieces needed for coherent
multiband parameter estimation combining a lunar detector with LISA-family
and ground-based detectors. It does not implement a lunar antenna-pattern
or single-link response function (the analogue of that work for LISA/Taiji/
TianQin is also still pending); this is purely the coordinate layer.

Design, following Tissino et al. 2026 ("The geometry of lunar gravitational
wave detection", arXiv:2606.04918), which shows (their Eq. 7-9) that a
source's sky position and polarization angle do not need to be re-expressed
between solar-system-barycentric-comoving frames -- "these vectors are
always constant for a given source, as motion within the Solar System is
negligible compared to the source distance". Accordingly, the SSB -> Moon
rotation used here is the identity: `longitude`/`latitude`/`polarization`
come out numerically unchanged between the SSB and Moon frames, and only
the arrival time changes. This is analogous to how `ssb_to_lisa` does not
apply any parallax correction to sky position either -- only to arrival
time.

What *does* need the Moon's precise position (and, for a specific surface
site, its libration) is the arrival time itself (Tissino et al. Eq. 12-13):
`moon_site_position_ssb` supports both the Moon's center (real ephemeris
only, no extra dependency) and a specific surface site (via the optional
`lunarsky` package, which models real lunar libration). The lunar body's
*orientation* (needed for an actual antenna-pattern/response function, as
in Tissino et al. Eq. 13's `n`, `b` vectors) is not modeled here at all --
that is response-function work, not coordinate-layer work.
"""

import numpy as np
from astropy.constants import c
from scipy.optimize import fsolve

from pycbc.coordinates.space import (
    geo_to_ssb,
    lisa_to_ssb,
    localization_to_propagation_vector,
    polarization_newframe,
    propagation_vector_to_localization,
    ssb_to_geo,
    ssb_to_lisa,
)

__all__ = [
    'rotation_matrix_ssb_to_moon', 'moon_site_position_ssb',
    't_moon_from_ssb', 't_ssb_from_t_moon', 'ssb_to_moon', 'moon_to_ssb',
    'moon_to_geo', 'geo_to_moon', 'moon_to_lisa', 'lisa_to_moon',
]


def rotation_matrix_ssb_to_moon():
    """The rotation matrix (of frame basis) from the SSB frame to the
    lunar frame used throughout this module.

    This is the identity: the lunar frame here is ecliptic-axis-aligned
    with the SSB frame (translated to the Moon's position, not rotated) --
    the Moon's own body orientation/libration is not modeled at this
    coordinate-unification layer (see module docstring). A future lunar
    antenna-pattern/response function would need a real, libration-aware,
    body-fixed rotation instead; this placeholder is not it.

    Returns
    -------
    r : numpy.array
        A 3x3 identity matrix.
    """
    return np.eye(3)


def moon_site_position_ssb(t_moon, longitude=None, latitude=None,
                           height=0.0):
    """Position vector of a point on (or above) the Moon in the SSB frame,
    at a given time.

    Default `longitude`/`latitude` (None) gives the Moon's barycenter, via
    astropy's real solar-system ephemeris -- no extra dependency. An
    explicit selenodetic `longitude`/`latitude` instead gives the precise
    position of that surface site, including real lunar libration, via
    the optional `lunarsky` package (only imported if a site is given).

    Parameters
    ----------
    t_moon : float
        GPS time [s].
    longitude, latitude : float or None, optional
        Selenodetic longitude/latitude of a surface site [rad]. Default
        None (both must be None, or both given), which gives the Moon's
        barycenter.
    height : float, optional
        Height above the reference selenoid [m]. Only used if
        `longitude`/`latitude` are given. Default 0.0.

    Returns
    -------
    p : numpy.array
        Position vector in the SSB frame, shape (3, 1), in meters.
    """
    if longitude is None and latitude is None:
        from pycbc.coordinates.space_orbit import _real_body_position_velocity
        pos, _ = _real_body_position_velocity(t_moon, 'moon')
        return pos.reshape(3, 1)

    if longitude is None or latitude is None:
        raise ValueError(
            "longitude and latitude must either both be given (a specific "
            "surface site), or both left as None (the Moon's barycenter).")

    try:
        from lunarsky import MoonLocation
    except ImportError as exc:
        raise ImportError(
            "lunarsky is required for moon_site_position_ssb with an "
            "explicit site location (it models real lunar libration); "
            "pip install lunarsky, or omit longitude/latitude to use the "
            "Moon's barycenter instead (no lunarsky needed).") from exc
    from astropy import units as apy_units
    from astropy.coordinates import ICRS
    from astropy.time import Time

    from pycbc.coordinates.space_orbit import _icrs_to_ecliptic_rotation_matrix

    site = MoonLocation.from_selenodetic(
        lon=longitude * apy_units.rad, lat=latitude * apy_units.rad,
        height=height * apy_units.m)
    icrs = site.get_mcmf(Time(t_moon, format='gps')).transform_to(ICRS())
    pos_icrs = np.array([
        icrs.cartesian.x.to(apy_units.m).value,
        icrs.cartesian.y.to(apy_units.m).value,
        icrs.cartesian.z.to(apy_units.m).value])
    rotation = _icrs_to_ecliptic_rotation_matrix()
    return (rotation @ pos_icrs).reshape(3, 1)


def t_moon_from_ssb(t_ssb, longitude_ssb, latitude_ssb,
                    longitude_site=None, latitude_site=None):
    """ Arrival time at a point on the Moon, from arrival time and sky
    localization in the SSB frame.

    Parameters
    ----------
    t_ssb : float
        Arrival time at the SSB frame origin [s].
    longitude_ssb : float
        Ecliptic longitude in the SSB frame [rad].
    latitude_ssb : float
        Ecliptic latitude in the SSB frame [rad].
    longitude_site, latitude_site : float or None, optional
        See `moon_site_position_ssb`. Default None (the Moon's
        barycenter).

    Returns
    -------
    t_moon : float
        Arrival time at the given point on the Moon [s].
    """
    k = localization_to_propagation_vector(
            longitude_ssb, latitude_ssb, use_astropy=False)

    def equation(t_moon):
        # the Moon is moving, when GW arrives at the given point,
        # time is t_moon, not t_ssb.
        p = moon_site_position_ssb(t_moon, longitude_site, latitude_site)
        return t_moon - t_ssb - np.vdot(k, p) / c.value

    return fsolve(equation, t_ssb)[0]


def t_ssb_from_t_moon(t_moon, longitude_ssb, latitude_ssb,
                      longitude_site=None, latitude_site=None):
    """ Arrival time at the SSB frame origin, from arrival time at a point
    on the Moon and sky localization in the SSB frame.

    Parameters
    ----------
    t_moon : float
        Arrival time at the given point on the Moon [s].
    longitude_ssb : float
        Ecliptic longitude in the SSB frame [rad].
    latitude_ssb : float
        Ecliptic latitude in the SSB frame [rad].
    longitude_site, latitude_site : float or None, optional
        See `moon_site_position_ssb`. Default None (the Moon's
        barycenter).

    Returns
    -------
    t_ssb : float
        Arrival time at the SSB frame origin [s].
    """
    k = localization_to_propagation_vector(
            longitude_ssb, latitude_ssb, use_astropy=False)

    # the Moon is moving, when GW arrives at the given point,
    # time is t_moon, not t_ssb.
    p = moon_site_position_ssb(t_moon, longitude_site, latitude_site)

    def equation(t_ssb):
        return t_moon - t_ssb - np.vdot(k, p) / c.value

    return fsolve(equation, t_moon)[0]


def ssb_to_moon(t_ssb, longitude_ssb, latitude_ssb, polarization_ssb,
                longitude_site=None, latitude_site=None):
    """ Converts arrival time, sky localization, and polarization from the
    SSB frame to a lunar frame.

    Because `rotation_matrix_ssb_to_moon` is the identity (see module
    docstring), `longitude`/`latitude`/`polarization` come out numerically
    identical to the SSB-frame input; only the arrival time changes.

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
    longitude_site, latitude_site : float or None, optional
        See `moon_site_position_ssb`. Default None (the Moon's
        barycenter).

    Returns
    -------
    (t_moon, longitude_moon, latitude_moon, polarization_moon) : tuple
        Arrival time [s], longitude [rad], latitude [rad], and
        polarization angle [rad] in the lunar frame.
    """
    if not isinstance(t_ssb, np.ndarray):
        t_ssb = np.array([t_ssb])
    if not isinstance(longitude_ssb, np.ndarray):
        longitude_ssb = np.array([longitude_ssb])
    if not isinstance(latitude_ssb, np.ndarray):
        latitude_ssb = np.array([latitude_ssb])
    if not isinstance(polarization_ssb, np.ndarray):
        polarization_ssb = np.array([polarization_ssb])
    num = len(t_ssb)
    t_moon = np.full(num, np.nan)
    longitude_moon = np.full(num, np.nan)
    latitude_moon = np.full(num, np.nan)
    polarization_moon = np.full(num, np.nan)

    rotation_matrix = rotation_matrix_ssb_to_moon()

    for i in range(num):
        if longitude_ssb[i] < 0 or longitude_ssb[i] >= 2*np.pi:
            raise ValueError("Longitude should within [0, 2*pi).")
        if latitude_ssb[i] < -np.pi/2 or latitude_ssb[i] > np.pi/2:
            raise ValueError("Latitude should within [-pi/2, pi/2].")
        if polarization_ssb[i] < 0 or polarization_ssb[i] >= 2*np.pi:
            raise ValueError("Polarization angle should within [0, 2*pi).")
        t_moon[i] = t_moon_from_ssb(
            t_ssb[i], longitude_ssb[i], latitude_ssb[i],
            longitude_site, latitude_site)
        k_ssb = localization_to_propagation_vector(
                    longitude_ssb[i], latitude_ssb[i], use_astropy=False)
        k_moon = rotation_matrix.T @ k_ssb
        longitude_moon[i], latitude_moon[i] = \
            propagation_vector_to_localization(k_moon, use_astropy=False)
        polarization_moon[i] = polarization_newframe(
            polarization_ssb[i], k_ssb, rotation_matrix, use_astropy=False)

    if num == 1:
        params_moon = (t_moon[0], longitude_moon[0],
                      latitude_moon[0], polarization_moon[0])
    else:
        params_moon = (t_moon, longitude_moon,
                      latitude_moon, polarization_moon)

    return params_moon


def moon_to_ssb(t_moon, longitude_moon, latitude_moon, polarization_moon,
                longitude_site=None, latitude_site=None):
    """ Converts arrival time, sky localization, and polarization from a
    lunar frame to the SSB frame. Inverse of `ssb_to_moon`.

    Parameters
    ----------
    t_moon : float or numpy.array
        Arrival time at the given point on the Moon [s].
    longitude_moon : float or numpy.array
        Longitude in the lunar frame [rad].
    latitude_moon : float or numpy.array
        Latitude in the lunar frame [rad].
    polarization_moon : float or numpy.array
        Polarization angle in the lunar frame [rad].
    longitude_site, latitude_site : float or None, optional
        See `moon_site_position_ssb`. Default None (the Moon's
        barycenter).

    Returns
    -------
    (t_ssb, longitude_ssb, latitude_ssb, polarization_ssb) : tuple
        Arrival time [s], ecliptic longitude [rad], ecliptic latitude
        [rad], and polarization angle [rad] in the SSB frame.
    """
    if not isinstance(t_moon, np.ndarray):
        t_moon = np.array([t_moon])
    if not isinstance(longitude_moon, np.ndarray):
        longitude_moon = np.array([longitude_moon])
    if not isinstance(latitude_moon, np.ndarray):
        latitude_moon = np.array([latitude_moon])
    if not isinstance(polarization_moon, np.ndarray):
        polarization_moon = np.array([polarization_moon])
    num = len(t_moon)
    t_ssb = np.full(num, np.nan)
    longitude_ssb = np.full(num, np.nan)
    latitude_ssb = np.full(num, np.nan)
    polarization_ssb = np.full(num, np.nan)

    rotation_matrix = rotation_matrix_ssb_to_moon()

    for i in range(num):
        if longitude_moon[i] < 0 or longitude_moon[i] >= 2*np.pi:
            raise ValueError("Longitude should within [0, 2*pi).")
        if latitude_moon[i] < -np.pi/2 or latitude_moon[i] > np.pi/2:
            raise ValueError("Latitude should within [-pi/2, pi/2].")
        if polarization_moon[i] < 0 or polarization_moon[i] >= 2*np.pi:
            raise ValueError("Polarization angle should within [0, 2*pi).")
        k_moon = localization_to_propagation_vector(
                    longitude_moon[i], latitude_moon[i], use_astropy=False)
        k_ssb = rotation_matrix @ k_moon
        longitude_ssb[i], latitude_ssb[i] = \
            propagation_vector_to_localization(k_ssb, use_astropy=False)
        polarization_ssb[i] = polarization_newframe(
            polarization_moon[i], k_moon, rotation_matrix.T,
            use_astropy=False)
        t_ssb[i] = t_ssb_from_t_moon(
            t_moon[i], longitude_ssb[i], latitude_ssb[i],
            longitude_site, latitude_site)

    if num == 1:
        params_ssb = (t_ssb[0], longitude_ssb[0],
                     latitude_ssb[0], polarization_ssb[0])
    else:
        params_ssb = (t_ssb, longitude_ssb,
                     latitude_ssb, polarization_ssb)

    return params_ssb


def moon_to_geo(t_moon, longitude_moon, latitude_moon, polarization_moon,
                longitude_site=None, latitude_site=None, use_astropy=True):
    """ Converting the arrive time, the sky localization, and the
    polarization from a lunar frame to the geocentric frame.

    Parameters
    ----------
    t_moon, longitude_moon, latitude_moon, polarization_moon :
        See `moon_to_ssb`.
    longitude_site, latitude_site : float or None, optional
        See `moon_site_position_ssb`. Default None (the Moon's
        barycenter).
    use_astropy : bool, optional
        Using Astropy to calculate the sky localization or not.
        Default is True.

    Returns
    -------
    (t_geo, longitude_geo, latitude_geo, polarization_geo) : tuple
        See `pycbc.coordinates.space.ssb_to_geo`.
    """
    t_ssb, longitude_ssb, latitude_ssb, polarization_ssb = moon_to_ssb(
        t_moon, longitude_moon, latitude_moon, polarization_moon,
        longitude_site, latitude_site)
    return ssb_to_geo(
        t_ssb, longitude_ssb, latitude_ssb, polarization_ssb, use_astropy)


def geo_to_moon(t_geo, longitude_geo, latitude_geo, polarization_geo,
                longitude_site=None, latitude_site=None, use_astropy=True):
    """ Converting the arrive time, the sky localization, and the
    polarization from the geocentric frame to a lunar frame.

    Parameters
    ----------
    t_geo, longitude_geo, latitude_geo, polarization_geo :
        See `pycbc.coordinates.space.geo_to_ssb`.
    longitude_site, latitude_site : float or None, optional
        See `moon_site_position_ssb`. Default None (the Moon's
        barycenter).
    use_astropy : bool, optional
        Using Astropy to calculate the sky localization or not.
        Default is True.

    Returns
    -------
    (t_moon, longitude_moon, latitude_moon, polarization_moon) : tuple
        See `ssb_to_moon`.
    """
    t_ssb, longitude_ssb, latitude_ssb, polarization_ssb = geo_to_ssb(
        t_geo, longitude_geo, latitude_geo, polarization_geo, use_astropy)
    return ssb_to_moon(
        t_ssb, longitude_ssb, latitude_ssb, polarization_ssb,
        longitude_site, latitude_site)


def moon_to_lisa(t_moon, longitude_moon, latitude_moon, polarization_moon,
                 longitude_site=None, latitude_site=None,
                 t0=None, orbit=None, sc=(1, 2, 3)):
    """ Converting the arrive time, the sky localization, and the
    polarization from a lunar frame to the LISA (or, via `orbit`, any
    constellation) frame.

    Parameters
    ----------
    t_moon, longitude_moon, latitude_moon, polarization_moon :
        See `moon_to_ssb`.
    longitude_site, latitude_site : float or None, optional
        See `moon_site_position_ssb`. Default None (the Moon's
        barycenter).
    t0 : float or None, optional
        See `pycbc.coordinates.space.ssb_to_lisa`. Default None, which
        uses that function's own default.
    orbit : OrbitProvider, optional
        See `pycbc.coordinates.space.ssb_to_lisa`. Default None (analytic
        circular LISA orbit); passing a Taiji/TianQin/numeric orbit here
        makes this a Moon<->Taiji or Moon<->TianQin transform.
    sc : tuple, optional
        1-indexed spacecraft labels defining the constellation. Only used
        if `orbit` is given. Default (1, 2, 3).

    Returns
    -------
    (t_lisa, longitude_lisa, latitude_lisa, polarization_lisa) : tuple
        See `pycbc.coordinates.space.ssb_to_lisa`.
    """
    t_ssb, longitude_ssb, latitude_ssb, polarization_ssb = moon_to_ssb(
        t_moon, longitude_moon, latitude_moon, polarization_moon,
        longitude_site, latitude_site)
    kwargs = {'orbit': orbit, 'sc': sc}
    if t0 is not None:
        kwargs['t0'] = t0
    return ssb_to_lisa(
        t_ssb, longitude_ssb, latitude_ssb, polarization_ssb, **kwargs)


def lisa_to_moon(t_lisa, longitude_lisa, latitude_lisa, polarization_lisa,
                 longitude_site=None, latitude_site=None,
                 t0=None, orbit=None, sc=(1, 2, 3)):
    """ Converting the arrive time, the sky localization, and the
    polarization from the LISA (or, via `orbit`, any constellation) frame
    to a lunar frame.

    Parameters
    ----------
    t_lisa, longitude_lisa, latitude_lisa, polarization_lisa :
        See `pycbc.coordinates.space.lisa_to_ssb`.
    longitude_site, latitude_site : float or None, optional
        See `moon_site_position_ssb`. Default None (the Moon's
        barycenter).
    t0 : float or None, optional
        See `pycbc.coordinates.space.lisa_to_ssb`. Default None, which
        uses that function's own default.
    orbit : OrbitProvider, optional
        See `pycbc.coordinates.space.lisa_to_ssb`. Default None (analytic
        circular LISA orbit); passing a Taiji/TianQin/numeric orbit here
        makes this a Taiji<->Moon or TianQin<->Moon transform.
    sc : tuple, optional
        1-indexed spacecraft labels defining the constellation. Only used
        if `orbit` is given. Default (1, 2, 3).

    Returns
    -------
    (t_moon, longitude_moon, latitude_moon, polarization_moon) : tuple
        See `ssb_to_moon`.
    """
    kwargs = {'orbit': orbit, 'sc': sc}
    if t0 is not None:
        kwargs['t0'] = t0
    t_ssb, longitude_ssb, latitude_ssb, polarization_ssb = lisa_to_ssb(
        t_lisa, longitude_lisa, latitude_lisa, polarization_lisa, **kwargs)
    return ssb_to_moon(
        t_ssb, longitude_ssb, latitude_ssb, polarization_ssb,
        longitude_site, latitude_site)
