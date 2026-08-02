# Copyright (C) 2026  Shichao Wu, Alex Nitz, Jacopo Tissino, Jan Harms
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
transforms to an arbitrary user-specified *fixed* reference point ("Ref"
frame) -- for arrival time, sky localization, and polarization angle, the
same three quantities `pycbc.coordinates.moon`/`space` unify across
SSB/GEO/LISA-family/Moon frames.

Motivation: for coherent multiband parameter estimation combining detectors
that are very far apart (e.g. a lunar detector, LISA, and a ground-based
detector), the choice of timing reference point matters a lot for sampling
efficiency, even though every such choice is physically equivalent. Tissino
et al. 2026 ("The geometry of lunar gravitational wave detection",
arXiv:2606.04918), Table 1 and Sec. 3.2, show the arrival-time uncertainty
for the same lunar-detector source varies from ~11.66s (SSB-referenced) to
~0.25s (Moon-barycenter-referenced) to ~0.1s (a numerically-optimized
reference point) -- and stress that "the Earth and Moon are not comoving
with the SSB, and are therefore not valid origins" for that kind of
optimization: only a *fixed* (SSB-comoving) point is physically equivalent
to any other fixed point, which is exactly the property this module
requires and exploits.

Because the "Ref" frame is defined to be a fixed point (not, e.g., a
detector's own instantaneous, moving position), this is *simpler* than the
Moon/GEO/LISA frames already implemented elsewhere in this package:

- Sky localization and polarization angle are unchanged between SSB and Ref
  (the same Eq. 7-9 argument as `pycbc.coordinates.moon`'s module
  docstring): the rotation used here is the identity.
- Arrival time conversion is a closed-form linear light-travel-time
  correction (`t_ref = t_ssb + dot(k, ref_position) / c`) -- unlike
  `pycbc.coordinates.moon`/`space`'s `t_X_from_ssb` functions, which need
  `scipy.optimize.fsolve` because the Moon/Earth/LISA constellation are all
  *moving* targets. A fixed point needs no root-finding at all.

`ref_position` must be a fixed (time-independent) 3-vector position, in
meters, in the SSB ecliptic frame -- e.g. a specific numerically-optimized
point (Tissino et al. Sec. 3.2), or a named body's position frozen at one
fixed reference time (e.g. `pycbc.coordinates.space.earth_position_ssb`
evaluated once). Resolving a *name* (like "the optimal point" or
"geocenter") into a concrete `ref_position` vector is the caller's
responsibility -- this module only implements the (frame-generic) transform
machinery once that vector is known.
"""

import numpy as np
from astropy.constants import c

from pycbc.coordinates.moon import moon_to_ssb, ssb_to_moon
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
    "rotation_matrix_ssb_to_ref",
    "t_ref_from_ssb",
    "t_ssb_from_t_ref",
    "ssb_to_ref",
    "ref_to_ssb",
    "ref_to_geo",
    "geo_to_ref",
    "ref_to_moon",
    "moon_to_ref",
    "ref_to_lisa",
    "lisa_to_ref",
]


def rotation_matrix_ssb_to_ref():
    """The rotation matrix (of frame basis) from the SSB frame to the
    fixed "Ref" reference-point frame used throughout this module.

    This is the identity, for the same reason as
    `pycbc.coordinates.moon.rotation_matrix_ssb_to_moon`: a fixed point's
    sky localization/polarization representation does not depend on its
    (fixed) position within the Solar System (Tissino et al. 2026 Eq.
    7-9). See the module docstring.

    Returns
    -------
    r : numpy.array
        A 3x3 identity matrix.
    """
    return np.eye(3)


def t_ref_from_ssb(t_ssb, longitude_ssb, latitude_ssb, ref_position):
    """Calculating the time when a GW signal arrives at a fixed
    reference point, by using the time and sky localization in the SSB
    frame.

    Closed form: unlike `pycbc.coordinates.moon.t_moon_from_ssb` (which
    needs `scipy.optimize.fsolve` because the Moon is moving),
    `ref_position` is fixed, so this is a single light-travel-time
    correction with no root-finding needed.

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
    ref_position : numpy.array
        The fixed (time-independent) position of the reference point in
        the SSB frame, shape (3,), in the unit of 'm'. Must be a point
        comoving with the SSB (i.e. genuinely fixed, not e.g. a moving
        detector's instantaneous position) -- see the module docstring.

    Returns
    -------
    t_ref : float or numpy.array
        The time when a GW signal arrives at the reference point.
    """
    if not isinstance(t_ssb, np.ndarray):
        t_ssb = np.array([t_ssb])
        scalar = True
    else:
        scalar = False
    longitude_ssb = np.atleast_1d(longitude_ssb)
    latitude_ssb = np.atleast_1d(latitude_ssb)
    ref_position = np.asarray(ref_position, dtype=float)

    t_ref = np.empty_like(t_ssb, dtype=float)
    for i in range(len(t_ssb)):
        k = localization_to_propagation_vector(
            longitude_ssb[i], latitude_ssb[i], use_astropy=False
        )
        t_ref[i] = t_ssb[i] + np.vdot(k, ref_position) / c.value

    if scalar or len(t_ref) == 1:
        return t_ref[0]
    return t_ref


def t_ssb_from_t_ref(t_ref, longitude_ssb, latitude_ssb, ref_position):
    """Calculating the time when a GW signal arrives at the barycenter
    of SSB, by using the time at a fixed reference point and sky
    localization in the SSB frame. Inverse of `t_ref_from_ssb`, also
    closed form.

    Parameters
    ----------
    t_ref : float or numpy.array
        The time when a GW signal arrives at the reference point.
        In the unit of 's'.
    longitude_ssb : float or numpy.array
        The ecliptic longitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    latitude_ssb : float or numpy.array
        The ecliptic latitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    ref_position : numpy.array
        See `t_ref_from_ssb`.

    Returns
    -------
    t_ssb : float or numpy.array
        The time when a GW signal arrives at the origin of SSB frame.
    """
    if not isinstance(t_ref, np.ndarray):
        t_ref = np.array([t_ref])
        scalar = True
    else:
        scalar = False
    longitude_ssb = np.atleast_1d(longitude_ssb)
    latitude_ssb = np.atleast_1d(latitude_ssb)
    ref_position = np.asarray(ref_position, dtype=float)

    t_ssb = np.empty_like(t_ref, dtype=float)
    for i in range(len(t_ref)):
        k = localization_to_propagation_vector(
            longitude_ssb[i], latitude_ssb[i], use_astropy=False
        )
        t_ssb[i] = t_ref[i] - np.vdot(k, ref_position) / c.value

    if scalar or len(t_ssb) == 1:
        return t_ssb[0]
    return t_ssb


def ssb_to_ref(t_ssb, longitude_ssb, latitude_ssb, polarization_ssb, ref_position):
    """Converting the arrive time, the sky localization, and the
    polarization from the SSB frame to a fixed reference-point frame.

    Because `rotation_matrix_ssb_to_ref` is the identity (see module
    docstring), `longitude`/`latitude`/`polarization` come out
    numerically identical to the SSB-frame input; only the arrival time
    changes.

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
    ref_position : numpy.array
        See `t_ref_from_ssb`.

    Returns
    -------
    (t_ref, longitude_ref, latitude_ref, polarization_ref) : tuple
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
    t_ref = np.full(num, np.nan)
    longitude_ref = np.full(num, np.nan)
    latitude_ref = np.full(num, np.nan)
    polarization_ref = np.full(num, np.nan)

    rotation_matrix = rotation_matrix_ssb_to_ref()

    for i in range(num):
        if longitude_ssb[i] < 0 or longitude_ssb[i] >= 2 * np.pi:
            raise ValueError("Longitude should within [0, 2*pi).")
        if latitude_ssb[i] < -np.pi / 2 or latitude_ssb[i] > np.pi / 2:
            raise ValueError("Latitude should within [-pi/2, pi/2].")
        if polarization_ssb[i] < 0 or polarization_ssb[i] >= 2 * np.pi:
            raise ValueError("Polarization angle should within [0, 2*pi).")
        t_ref[i] = t_ref_from_ssb(
            t_ssb[i], longitude_ssb[i], latitude_ssb[i], ref_position
        )
        k_ssb = localization_to_propagation_vector(
            longitude_ssb[i], latitude_ssb[i], use_astropy=False
        )
        k_ref = rotation_matrix.T @ k_ssb
        longitude_ref[i], latitude_ref[i] = propagation_vector_to_localization(
            k_ref, use_astropy=False
        )
        polarization_ref[i] = polarization_newframe(
            polarization_ssb[i], k_ssb, rotation_matrix, use_astropy=False
        )

    if num == 1:
        return (t_ref[0], longitude_ref[0], latitude_ref[0], polarization_ref[0])
    return (t_ref, longitude_ref, latitude_ref, polarization_ref)


def ref_to_ssb(t_ref, longitude_ref, latitude_ref, polarization_ref, ref_position):
    """Converting the arrive time, the sky localization, and the
    polarization from a fixed reference-point frame to the SSB frame.
    Inverse of `ssb_to_ref`.

    Parameters
    ----------
    t_ref, longitude_ref, latitude_ref, polarization_ref :
        See `ssb_to_ref`'s return values.
    ref_position : numpy.array
        See `t_ref_from_ssb`.

    Returns
    -------
    (t_ssb, longitude_ssb, latitude_ssb, polarization_ssb) : tuple
    """
    if not isinstance(t_ref, np.ndarray):
        t_ref = np.array([t_ref])
    if not isinstance(longitude_ref, np.ndarray):
        longitude_ref = np.array([longitude_ref])
    if not isinstance(latitude_ref, np.ndarray):
        latitude_ref = np.array([latitude_ref])
    if not isinstance(polarization_ref, np.ndarray):
        polarization_ref = np.array([polarization_ref])
    num = len(t_ref)
    t_ssb = np.full(num, np.nan)
    longitude_ssb = np.full(num, np.nan)
    latitude_ssb = np.full(num, np.nan)
    polarization_ssb = np.full(num, np.nan)

    rotation_matrix = rotation_matrix_ssb_to_ref()

    for i in range(num):
        if longitude_ref[i] < 0 or longitude_ref[i] >= 2 * np.pi:
            raise ValueError("Longitude should within [0, 2*pi).")
        if latitude_ref[i] < -np.pi / 2 or latitude_ref[i] > np.pi / 2:
            raise ValueError("Latitude should within [-pi/2, pi/2].")
        if polarization_ref[i] < 0 or polarization_ref[i] >= 2 * np.pi:
            raise ValueError("Polarization angle should within [0, 2*pi).")
        k_ref = localization_to_propagation_vector(
            longitude_ref[i], latitude_ref[i], use_astropy=False
        )
        k_ssb = rotation_matrix @ k_ref
        longitude_ssb[i], latitude_ssb[i] = propagation_vector_to_localization(
            k_ssb, use_astropy=False
        )
        polarization_ssb[i] = polarization_newframe(
            polarization_ref[i], k_ref, rotation_matrix.T, use_astropy=False
        )
        t_ssb[i] = t_ssb_from_t_ref(
            t_ref[i], longitude_ssb[i], latitude_ssb[i], ref_position
        )

    if num == 1:
        return (t_ssb[0], longitude_ssb[0], latitude_ssb[0], polarization_ssb[0])
    return (t_ssb, longitude_ssb, latitude_ssb, polarization_ssb)


def ref_to_geo(
    t_ref, longitude_ref, latitude_ref, polarization_ref, ref_position, use_astropy=True
):
    """Converting the arrive time, the sky localization, and the
    polarization from a fixed reference-point frame to the geocentric
    frame.

    Parameters
    ----------
    t_ref, longitude_ref, latitude_ref, polarization_ref :
        See `ssb_to_ref`'s return values.
    ref_position : numpy.array
        See `t_ref_from_ssb`.
    use_astropy : bool, optional
        Using Astropy to calculate the sky localization or not.
        Default is True.

    Returns
    -------
    (t_geo, longitude_geo, latitude_geo, polarization_geo) : tuple
        See `pycbc.coordinates.space.ssb_to_geo`.
    """
    t_ssb, longitude_ssb, latitude_ssb, polarization_ssb = ref_to_ssb(
        t_ref, longitude_ref, latitude_ref, polarization_ref, ref_position
    )
    return ssb_to_geo(t_ssb, longitude_ssb, latitude_ssb, polarization_ssb, use_astropy)


def geo_to_ref(
    t_geo, longitude_geo, latitude_geo, polarization_geo, ref_position, use_astropy=True
):
    """Converting the arrive time, the sky localization, and the
    polarization from the geocentric frame to a fixed reference-point
    frame.

    Parameters
    ----------
    t_geo, longitude_geo, latitude_geo, polarization_geo :
        See `pycbc.coordinates.space.geo_to_ssb`.
    ref_position : numpy.array
        See `t_ref_from_ssb`.
    use_astropy : bool, optional
        Using Astropy to calculate the sky localization or not.
        Default is True.

    Returns
    -------
    (t_ref, longitude_ref, latitude_ref, polarization_ref) : tuple
        See `ssb_to_ref`.
    """
    t_ssb, longitude_ssb, latitude_ssb, polarization_ssb = geo_to_ssb(
        t_geo, longitude_geo, latitude_geo, polarization_geo, use_astropy
    )
    return ssb_to_ref(
        t_ssb, longitude_ssb, latitude_ssb, polarization_ssb, ref_position
    )


def ref_to_moon(
    t_ref,
    longitude_ref,
    latitude_ref,
    polarization_ref,
    ref_position,
    longitude_site=None,
    latitude_site=None,
):
    """Converting the arrive time, the sky localization, and the
    polarization from a fixed reference-point frame to a lunar frame.

    Parameters
    ----------
    t_ref, longitude_ref, latitude_ref, polarization_ref :
        See `ssb_to_ref`'s return values.
    ref_position : numpy.array
        See `t_ref_from_ssb`.
    longitude_site, latitude_site : float or None, optional
        See `pycbc.coordinates.moon.moon_site_position_ssb`. Default
        None (the Moon's barycenter).

    Returns
    -------
    (t_moon, longitude_moon, latitude_moon, polarization_moon) : tuple
        See `pycbc.coordinates.moon.ssb_to_moon`.
    """
    t_ssb, longitude_ssb, latitude_ssb, polarization_ssb = ref_to_ssb(
        t_ref, longitude_ref, latitude_ref, polarization_ref, ref_position
    )
    return ssb_to_moon(
        t_ssb,
        longitude_ssb,
        latitude_ssb,
        polarization_ssb,
        longitude_site,
        latitude_site,
    )


def moon_to_ref(
    t_moon,
    longitude_moon,
    latitude_moon,
    polarization_moon,
    ref_position,
    longitude_site=None,
    latitude_site=None,
):
    """Converting the arrive time, the sky localization, and the
    polarization from a lunar frame to a fixed reference-point frame.

    Parameters
    ----------
    t_moon, longitude_moon, latitude_moon, polarization_moon :
        See `pycbc.coordinates.moon.moon_to_ssb`.
    ref_position : numpy.array
        See `t_ref_from_ssb`.
    longitude_site, latitude_site : float or None, optional
        See `pycbc.coordinates.moon.moon_site_position_ssb`. Default
        None (the Moon's barycenter).

    Returns
    -------
    (t_ref, longitude_ref, latitude_ref, polarization_ref) : tuple
        See `ssb_to_ref`.
    """
    t_ssb, longitude_ssb, latitude_ssb, polarization_ssb = moon_to_ssb(
        t_moon,
        longitude_moon,
        latitude_moon,
        polarization_moon,
        longitude_site,
        latitude_site,
    )
    return ssb_to_ref(
        t_ssb, longitude_ssb, latitude_ssb, polarization_ssb, ref_position
    )


def ref_to_lisa(
    t_ref,
    longitude_ref,
    latitude_ref,
    polarization_ref,
    ref_position,
    t0=None,
    orbit=None,
    sc=(1, 2, 3),
):
    """Converting the arrive time, the sky localization, and the
    polarization from a fixed reference-point frame to the LISA (or, via
    `orbit`, any constellation) frame.

    Parameters
    ----------
    t_ref, longitude_ref, latitude_ref, polarization_ref :
        See `ssb_to_ref`'s return values.
    ref_position : numpy.array
        See `t_ref_from_ssb`.
    t0 : float or None, optional
        See `pycbc.coordinates.space.ssb_to_lisa`. Default None, which
        uses that function's own default.
    orbit : OrbitProvider, optional
        See `pycbc.coordinates.space.ssb_to_lisa`. Default None (analytic
        circular LISA orbit); passing a Taiji/TianQin/numeric orbit here
        makes this a Ref<->Taiji or Ref<->TianQin transform.
    sc : tuple, optional
        1-indexed spacecraft labels defining the constellation. Only used
        if `orbit` is given. Default (1, 2, 3).

    Returns
    -------
    (t_lisa, longitude_lisa, latitude_lisa, polarization_lisa) : tuple
        See `pycbc.coordinates.space.ssb_to_lisa`.
    """
    t_ssb, longitude_ssb, latitude_ssb, polarization_ssb = ref_to_ssb(
        t_ref, longitude_ref, latitude_ref, polarization_ref, ref_position
    )
    kwargs = {"orbit": orbit, "sc": sc}
    if t0 is not None:
        kwargs["t0"] = t0
    return ssb_to_lisa(t_ssb, longitude_ssb, latitude_ssb, polarization_ssb, **kwargs)


def lisa_to_ref(
    t_lisa,
    longitude_lisa,
    latitude_lisa,
    polarization_lisa,
    ref_position,
    t0=None,
    orbit=None,
    sc=(1, 2, 3),
):
    """Converting the arrive time, the sky localization, and the
    polarization from the LISA (or, via `orbit`, any constellation) frame
    to a fixed reference-point frame.

    Parameters
    ----------
    t_lisa, longitude_lisa, latitude_lisa, polarization_lisa :
        See `pycbc.coordinates.space.lisa_to_ssb`.
    ref_position : numpy.array
        See `t_ref_from_ssb`.
    t0 : float or None, optional
        See `pycbc.coordinates.space.lisa_to_ssb`. Default None, which
        uses that function's own default.
    orbit : OrbitProvider, optional
        See `pycbc.coordinates.space.lisa_to_ssb`. Default None (analytic
        circular LISA orbit); passing a Taiji/TianQin/numeric orbit here
        makes this a Taiji<->Ref or TianQin<->Ref transform.
    sc : tuple, optional
        1-indexed spacecraft labels defining the constellation. Only used
        if `orbit` is given. Default (1, 2, 3).

    Returns
    -------
    (t_ref, longitude_ref, latitude_ref, polarization_ref) : tuple
        See `ssb_to_ref`.
    """
    kwargs = {"orbit": orbit, "sc": sc}
    if t0 is not None:
        kwargs["t0"] = t0
    t_ssb, longitude_ssb, latitude_ssb, polarization_ssb = lisa_to_ssb(
        t_lisa, longitude_lisa, latitude_lisa, polarization_lisa, **kwargs
    )
    return ssb_to_ref(
        t_ssb, longitude_ssb, latitude_ssb, polarization_ssb, ref_position
    )
