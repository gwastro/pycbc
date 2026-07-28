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
instance (`EqualArmlengthOrbits`, `KeplerianOrbits`, `InterpolatedOrbits`,
`OEMOrbits`) passed in directly -- this module does not import or depend on
`lisaorbits` itself. Users who do not have `lisaorbits` installed can instead
use `NumericOrbits` below, which reads the same (times, positions[,
velocities]) input as `lisaorbits.InterpolatedOrbits` and interpolates using
only `scipy`.

This module does not change or replace anything in `pycbc.coordinates.space`;
it is purely additive.
"""
import numpy as np
from scipy.interpolate import make_interp_spline
from scipy.optimize import fsolve
from astropy.constants import c as SPEED_OF_LIGHT

logger = __import__('logging').getLogger(__name__)


class NumericOrbits:
    """Interpolate a numerical constellation orbit given as arrays of
    spacecraft position (and, optionally, velocity) samples in the SSB frame.

    This mirrors the input contract of `lisaorbits.InterpolatedOrbits`, using
    only `scipy` for the spline interpolation, so that it can be used as a
    drop-in `OrbitProvider` without requiring `lisaorbits` to be installed.

    Parameters
    ----------
    t_interp : (N,) array-like
        Interpolating SSB times [s] (must be strictly increasing).
    positions : (N, M, 3) array-like
        Spacecraft positions in the SSB frame [m], with dimensions
        (time, spacecraft, xyz).
    velocities : (N, M, 3) array-like or None, optional
        Spacecraft velocities in the SSB frame [m/s], same shape as
        `positions`. If None (the default), velocities are computed as the
        analytic derivative of the position spline.
    interp_order : int, optional
        Spline interpolation order. Default 5 (quintic), matching the
        default used by `lisaorbits.InterpolatedOrbits`.
    """

    def __init__(self, t_interp, positions, velocities=None, interp_order=5):
        self.t_interp = np.asarray(t_interp, dtype=float)
        self.positions = np.asarray(positions, dtype=float)

        if self.t_interp.ndim != 1:
            raise ValueError(
                f't_interp must be 1-dimensional, got shape '
                f'{self.t_interp.shape}')
        if self.positions.ndim != 3 or self.positions.shape[2] != 3:
            raise ValueError(
                'positions must have shape (N, M, 3), got shape '
                f'{self.positions.shape}')
        if self.positions.shape[0] != self.t_interp.shape[0]:
            raise ValueError(
                'positions and t_interp must agree along the time axis: '
                f'got {self.positions.shape[0]} and {self.t_interp.shape[0]}')

        self.num_sc = self.positions.shape[1]
        self.interp_order = int(interp_order)
        if self.interp_order < 3:
            raise ValueError(
                f'interp_order must be at least 3, got {interp_order}')

        def interpolate(y):
            return make_interp_spline(self.t_interp, y, k=self.interp_order)

        # one spline per (spacecraft, coordinate) pair
        self._pos_splines = [
            [interpolate(self.positions[:, i, j]) for j in range(3)]
            for i in range(self.num_sc)
        ]

        if velocities is not None:
            velocities = np.asarray(velocities, dtype=float)
            if velocities.shape != self.positions.shape:
                raise ValueError(
                    'velocities must have the same shape as positions: '
                    f'got {velocities.shape} and {self.positions.shape}')
            self._vel_splines = [
                [interpolate(velocities[:, i, j]) for j in range(3)]
                for i in range(self.num_sc)
            ]
        else:
            self._vel_splines = [
                [spline.derivative() for spline in row]
                for row in self._pos_splines
            ]

    def _sc_indices(self, sc):
        """Map 1-indexed spacecraft labels (matching the `lisaorbits`
        convention) to 0-indexed array positions. `sc=None` selects all
        spacecraft.
        """
        if sc is None:
            return np.arange(self.num_sc)
        sc = np.atleast_1d(np.asarray(sc))
        if np.any(sc < 1) or np.any(sc > self.num_sc):
            raise ValueError(
                f'spacecraft labels must be in [1, {self.num_sc}], got {sc}')
        return sc - 1

    def _evaluate(self, splines, t, sc):
        t = np.atleast_1d(np.asarray(t, dtype=float))
        idx = self._sc_indices(sc)
        out = np.empty((len(t), len(idx), 3))
        for k, i in enumerate(idx):
            for j in range(3):
                out[:, k, j] = splines[i][j](t)
        return out

    def compute_position(self, t, sc=None):
        """Spacecraft position(s) at time(s) `t`.

        Parameters
        ----------
        t : (N,) array-like
            SSB time [s].
        sc : (M,) array-like or None, optional
            1-indexed spacecraft label(s). Default all spacecraft.

        Returns
        -------
        (N, M, 3) ndarray
            Spacecraft position(s) in the SSB frame [m].
        """
        return self._evaluate(self._pos_splines, t, sc)

    def compute_velocity(self, t, sc=None):
        """Spacecraft velocity/ies at time(s) `t`. Same conventions as
        `compute_position`.

        Returns
        -------
        (N, M, 3) ndarray
            Spacecraft velocity/ies in the SSB frame [m/s].
        """
        return self._evaluate(self._vel_splines, t, sc)


def constellation_frame(t, orbit, sc=(1, 2, 3)):
    """Compute the instantaneous constellation centroid and the rotation
    matrix from the SSB frame to a frame comoving with the constellation,
    directly from the spacecraft positions supplied by `orbit`.

    This generalizes the fixed circular-orbit
    `rotation_matrix_ssb_to_lisa`/`lisa_position_ssb` functions in
    `pycbc.coordinates.space` (which assume an idealized, rigid, circular
    LISA orbit) to work with any 3-spacecraft triangular constellation, given
    either analytic or numerical spacecraft positions -- e.g. LISA, Taiji, or
    TianQin.

    For each time sample, the instantaneous frame is built directly from the
    three spacecraft positions as:

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

    # Fix the residual axis-labeling freedom (the constellation geometry
    # alone only pins down the frame up to a rotation about its own normal
    # and the sign of the in-plane axes) by matching the convention already
    # used by `rotation_matrix_ssb_to_lisa` in the circular-orbit limit: a
    # fixed 180 degree rotation of the in-plane (x, y) axes about the plane
    # normal. This has no effect on the physics (light travel times, sky
    # localization magnitude, ...), only on the labeling of the frame axes,
    # and is verified to reproduce `rotation_matrix_ssb_to_lisa` exactly for
    # the analytic circular orbit -- see test_coordinates_space_orbit.py.
    axis_convention = np.diag([-1.0, -1.0, 1.0])
    rotation = rotation @ axis_convention

    return centroid, rotation


def t_detector_from_ssb(t_ssb, k_ssb, orbit, sc=(1, 2, 3)):
    """Compute the time at which a GW signal arrives at the constellation
    centroid, given the time and propagation direction in the SSB frame.

    This generalizes `pycbc.coordinates.space.t_lisa_from_ssb` to accept an
    arbitrary orbit provider instead of the fixed circular LISA orbit. The
    constellation is moving, so the arrival time is found by solving the
    implicit light-travel-time equation self-consistently, exactly as in the
    original function.

    Parameters
    ----------
    t_ssb : float
        The time at which the GW signal arrives at the origin of the SSB
        frame [s].
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
        The time at which the GW signal arrives at the constellation
        centroid [s].
    """
    k_ssb = np.asarray(k_ssb, dtype=float)

    def equation(t_detector):
        centroid, _ = constellation_frame([t_detector[0]], orbit, sc=sc)
        return t_detector[0] - t_ssb - np.dot(k_ssb, centroid[0]) \
            / SPEED_OF_LIGHT.value

    return fsolve(equation, t_ssb)[0]


def t_ssb_from_t_detector(t_detector, k_ssb, orbit, sc=(1, 2, 3)):
    """Compute the time at which a GW signal arrives at the origin of the
    SSB frame, given the arrival time at the constellation centroid and the
    propagation direction in the SSB frame.

    Inverse of `t_detector_from_ssb`; see that function for parameter and
    return conventions.
    """
    k_ssb = np.asarray(k_ssb, dtype=float)
    centroid, _ = constellation_frame([t_detector], orbit, sc=sc)

    def equation(t_ssb):
        return t_detector - t_ssb[0] - np.dot(k_ssb, centroid[0]) \
            / SPEED_OF_LIGHT.value

    return fsolve(equation, t_detector)[0]


__all__ = [
    'NumericOrbits',
    'constellation_frame',
    't_detector_from_ssb',
    't_ssb_from_t_detector',
]
