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
can be used as an orbit provider, including a real `lisaorbits.Orbits`
instance passed in directly -- this module does not import or depend on
`lisaorbits` itself. Users who do not have `lisaorbits` installed can instead
use `NumericOrbits` below, which reads the same (times, positions[,
velocities]) input as `lisaorbits.InterpolatedOrbits` and interpolates using
only `scipy`.

`compute_velocity(t, sc)` and `compute_acceleration(t, sc)` (same calling
convention, [m/s] and [m/s^2]) are analytic derivatives of the interpolating
spline, mirroring `lisaorbits.InterpolatedOrbits`.
"""
import h5py
import numpy as np
from scipy.interpolate import make_interp_spline

logger = __import__('logging').getLogger(__name__)


class NumericOrbits:
    """Interpolate a numerical constellation orbit given as arrays of
    spacecraft position (and, optionally, velocity) samples in the SSB frame.

    Mirrors the input contract of `lisaorbits.InterpolatedOrbits`, using
    only `scipy` for the spline interpolation, so it works as a drop-in
    `OrbitProvider` without requiring `lisaorbits` to be installed.

    Parameters
    ----------
    t_interp : (N,) array-like
        Interpolating SSB times [s] (must be strictly increasing).
    positions : (N, M, 3) array-like
        Spacecraft positions in the SSB frame [m]: N time samples, M
        spacecraft (any number, not just 3 -- see `num_sc` below), 3
        Cartesian (x, y, z) coordinates each.
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
            # One vector-valued spline over all (spacecraft, xyz) columns
            # flattened together, not num_sc*3 separate scalar splines:
            # scipy shares the basis-function computation across columns
            # in a single BSpline call, ~5-8x faster with identical results.
            return make_interp_spline(
                self.t_interp, y.reshape(len(self.t_interp), -1),
                k=self.interp_order)

        self._pos_spline = interpolate(self.positions)

        if velocities is not None:
            velocities = np.asarray(velocities, dtype=float)
            if velocities.shape != self.positions.shape:
                raise ValueError(
                    'velocities must have the same shape as positions: '
                    f'got {velocities.shape} and {self.positions.shape}')
            self._vel_spline = interpolate(velocities)
        else:
            self._vel_spline = self._pos_spline.derivative()
        # Kept only so `to_file` can round-trip whether velocities were
        # explicitly supplied (written out) or derived (omitted, so a
        # reloaded instance re-derives them the same way this one did).
        self.velocities = velocities

        # Acceleration is always the derivative of the velocity spline
        # (provided or derived), matching lisaorbits.InterpolatedOrbits.
        self._acc_spline = self._vel_spline.derivative()

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

    def _evaluate(self, spline, t, sc):
        t = np.atleast_1d(np.asarray(t, dtype=float))
        idx = self._sc_indices(sc)
        values = spline(t).reshape(len(t), self.num_sc, 3)
        return values[:, idx, :]

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
        return self._evaluate(self._pos_spline, t, sc)

    def compute_velocity(self, t, sc=None):
        """Spacecraft velocity/ies at time(s) `t`. Same conventions as
        `compute_position`.

        Returns
        -------
        (N, M, 3) ndarray
            Spacecraft velocity/ies in the SSB frame [m/s].
        """
        return self._evaluate(self._vel_spline, t, sc)

    def compute_acceleration(self, t, sc=None):
        """Spacecraft acceleration(s) at time(s) `t`, as the analytic
        derivative of `compute_velocity`'s spline. Same conventions as
        `compute_position`.

        Returns
        -------
        (N, M, 3) ndarray
            Spacecraft acceleration(s) in the SSB frame [m/s^2].
        """
        return self._evaluate(self._acc_spline, t, sc)

    @classmethod
    def from_file(cls, path, group=None, interp_order=5):
        """Load a numerical constellation orbit from an HDF5 file.

        The file (or the given group within it) must contain datasets `t`
        (shape (N,), SSB time [s]) and `positions` (shape (N, M, 3), SSB-
        frame spacecraft positions [m]). An optional `velocities` dataset
        (same shape as `positions`) is used if present; otherwise
        velocities are the analytic derivative of the position spline, as
        in `__init__`. This is the file format LISA, Taiji, or TianQin (or
        a numerically-propagated orbit for any of them) can be supplied
        in, e.g. from a PE config's `orbit-file` option (see
        `pycbc.transforms`).

        Parameters
        ----------
        path : str
            Path to the HDF5 orbit file.
        group : str, optional
            Name of an HDF5 group within the file to read the datasets
            from. Default None (read from the file's root).
        interp_order : int, optional
            See `__init__`. Default 5.

        Returns
        -------
        NumericOrbits
            A `NumericOrbits` instance built from the file's contents.
        """
        with h5py.File(path, 'r') as f:
            node = f[group] if group is not None else f
            t_interp = node['t'][:]
            positions = node['positions'][:]
            velocities = node['velocities'][:] \
                if 'velocities' in node else None
        return cls(t_interp, positions, velocities=velocities,
                   interp_order=interp_order)

    def to_file(self, path, group=None, mode='w'):
        """Write this orbit to an HDF5 file in the format read by
        `from_file`, e.g. for use with a PE config's `orbit-file` option
        (see `pycbc.transforms`).

        Writes datasets `t` (SSB time [s]) and `positions` (SSB-frame
        spacecraft positions [m]). If this instance was constructed with
        explicit `velocities`, those are written too; otherwise the
        `velocities` dataset is omitted, so a round trip through
        `from_file` re-derives velocities from the position spline the
        same way this instance does.
        never stored -- as in `__init__`, it's always the analytic
        derivative of the velocity spline.

        Parameters
        ----------
        path : str
            Path to the HDF5 orbit file to write.
        group : str, optional
            Name of an HDF5 group within the file to write the datasets
            to. Default None (write to the file's root).
        mode : str, optional
            `h5py.File` mode. Default 'w' (create the file, truncating any
            existing file at `path`); use 'a' to add a group to an
            existing file.
        """
        with h5py.File(path, mode) as f:
            node = f.create_group(group) if group is not None else f
            node['t'] = self.t_interp
            node['positions'] = self.positions
            if self.velocities is not None:
                node['velocities'] = self.velocities


__all__ = [
    'NumericOrbits',
]
