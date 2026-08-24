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
This module generalizes the constellation-frame machinery in
`pycbc.coordinates.space` (which assumes a fixed, analytic, circular LISA
orbit) to work with an arbitrary triangular constellation orbit, analytic or
numerical, for any mission (LISA, Taiji, TianQin, ...).

Any object exposing `compute_position(t, sc)` with the same call signature
and return shape as `lisaorbits.Orbits` (https://pypi.org/project/lisaorbits)
can be used as an orbit provider here, including a real `lisaorbits.Orbits`
instance passed in directly -- this module does not import or depend on
`lisaorbits` itself.

This module does not change or replace anything in `pycbc.coordinates.space`;
it is purely additive.
"""
import numpy as np
from astropy.constants import c as SPEED_OF_LIGHT
from scipy.optimize import fsolve

logger = __import__('logging').getLogger(__name__)


def constellation_frame(t, orbit, sc=(1, 2, 3)):
    """Compute the instantaneous constellation centroid and the rotation
    matrix from the SSB frame to a frame comoving with the constellation,
    directly from the spacecraft positions supplied by `orbit`.

    Generalizes `rotation_matrix_ssb_to_lisa`/`lisa_position_ssb` in
    `pycbc.coordinates.space`, which assume a fixed circular LISA orbit, to
    any 3-spacecraft triangular constellation.

    The frame is built at each time sample from the three positions as:

    * the centroid of the three spacecraft;
    * the plane normal (right-handed with respect to spacecraft ordering
      1 -> 2 -> 3), taken as the z-axis;
    * the projection onto that plane of the direction from the centroid to
      spacecraft 1, taken as the x-axis;
    * completing a right-handed orthonormal triad for the y-axis.

    Parameters
    ----------
    t : (N,) array-like
        SSB time(s) [s] at which to evaluate the frame.
    orbit : OrbitProvider
        Any object exposing `compute_position(t, sc)` with the calling
        convention of `lisaorbits.Orbits` (e.g. a `lisaorbits.Orbits`
        instance, or a `NumericOrbits` instance).
    sc : (3,) array-like, optional
        1-indexed spacecraft labels identifying the 3 constellation
        vertices, in the order used to define the plane normal and the
        reference (x-axis) spacecraft. Default (1, 2, 3).

    Returns
    -------
    centroid : (N, 3) ndarray
        Constellation centroid position in the SSB frame [m].
    rotation : (N, 3, 3) ndarray
        Rotation matrix from the SSB frame to the constellation frame at
        each time sample, such that ``v_constellation = rotation @ v_ssb``.
    """
    t = np.atleast_1d(np.asarray(t, dtype=float))
    pos = orbit.compute_position(t, sc)  # (N, 3, 3)
    r1, r2, r3 = pos[:, 0], pos[:, 1], pos[:, 2]

    centroid = (r1 + r2 + r3) / 3.0

    normal = np.cross(r2 - r1, r3 - r1)
    normal = normal / np.linalg.norm(normal, axis=-1, keepdims=True)

    x_raw = r1 - centroid
    # remove any out-of-plane component so the x-axis is exactly in-plane
    x_axis = x_raw - np.sum(x_raw * normal, axis=-1, keepdims=True) * normal
    x_axis = x_axis / np.linalg.norm(x_axis, axis=-1, keepdims=True)

    y_axis = np.cross(normal, x_axis)

    # rows are the new (constellation-frame) basis vectors, expressed in the
    # SSB basis, so that `rotation @ v_ssb` gives the components of `v_ssb`
    # in the constellation frame.
    rotation = np.stack([x_axis, y_axis, normal], axis=1)  # (N, 3, 3)

    # The constellation geometry alone only pins the frame down to a
    # rotation about its own normal and the sign of the in-plane axes; fix
    # that by matching rotation_matrix_ssb_to_lisa's convention (a 180
    # degree rotation about the normal). Pure axis-labeling, no physics
    # effect -- verified to reproduce rotation_matrix_ssb_to_lisa exactly
    # in the circular-orbit limit, see test_coordinates_space_orbit.py.
    axis_convention = np.diag([-1.0, -1.0, 1.0])
    rotation = rotation @ axis_convention

    return centroid, rotation


def _solve_frame_arrival_time(t_known, k_ssb, position_fn, forward):
    """Solve `t - t_other - dot(k, position(t)) / c = 0` for whichever side
    is unknown, i.e. the light-travel-time relationship between the SSB
    origin and a moving frame's origin.

    Parameters
    ----------
    t_known : float
        The time at the known side of the relationship [s].
    k_ssb : array-like
        The unit propagation vector of the GW signal in the SSB frame.
    position_fn : callable
        `position_fn(t)` returns the moving frame's origin position [m] in
        the SSB frame at SSB time `t`.
    forward : bool
        True solves for the frame's own time given `t_known = t_ssb`; the
        frame position then depends on the unknown, so `position_fn` is
        re-evaluated inside the root solve. False solves for `t_ssb`, where
        the frame position is known up front and evaluated once.

    Returns
    -------
    float
        The time at the unknown side of the relationship [s].
    """
    k_ssb = np.ravel(np.asarray(k_ssb, dtype=float))
    if forward:
        def equation(t):
            p = np.ravel(np.asarray(position_fn(t[0]), dtype=float))
            return t[0] - t_known - np.dot(k_ssb, p) / SPEED_OF_LIGHT.value
    else:
        p = np.ravel(np.asarray(position_fn(t_known), dtype=float))

        def equation(t):
            return t_known - t[0] - np.dot(k_ssb, p) / SPEED_OF_LIGHT.value
    return fsolve(equation, t_known)[0]


def t_detector_from_ssb(t_ssb, k_ssb, orbit, sc=(1, 2, 3)):
    """Compute the time at which a GW signal arrives at the constellation
    centroid, given the time and propagation direction in the SSB frame.

    Generalizes `pycbc.coordinates.space.t_lisa_from_ssb` to any orbit
    provider. The constellation moves, so the arrival time solves the
    implicit light-travel-time equation, as in the original.

    Parameters
    ----------
    t_ssb : float
        Arrival time at the SSB frame origin [s].
    k_ssb : (3,) array-like
        The unit propagation vector of the GW signal in the SSB frame.
    orbit : OrbitProvider
        Any object exposing `compute_position(t, sc)`, see
        `constellation_frame`.
    sc : (3,) array-like, optional
        1-indexed spacecraft labels defining the constellation. Default
        (1, 2, 3).

    Returns
    -------
    float
        Arrival time at the constellation centroid [s].
    """
    def position_fn(t):
        return constellation_frame([t], orbit, sc=sc)[0][0]
    return _solve_frame_arrival_time(t_ssb, k_ssb, position_fn, forward=True)


def t_ssb_from_t_detector(t_detector, k_ssb, orbit, sc=(1, 2, 3)):
    """Compute the time at which a GW signal arrives at the origin of the
    SSB frame, given the arrival time at the constellation centroid and the
    propagation direction in the SSB frame.

    Inverse of `t_detector_from_ssb`; see that function for parameter and
    return conventions.
    """
    def position_fn(t):
        return constellation_frame([t], orbit, sc=sc)[0][0]
    return _solve_frame_arrival_time(
        t_detector, k_ssb, position_fn, forward=False)

__all__ = [
    'constellation_frame',
    't_detector_from_ssb',
    't_ssb_from_t_detector',
]
