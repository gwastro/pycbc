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

`compute_velocity(t, sc)` and `compute_acceleration(t, sc)` (same calling
convention, [m/s] and [m/s^2] respectively) are also provided everywhere in
this module -- `NumericOrbits` and the analytic mission classes below --
needed by single-link response models that depend on spacecraft velocity
(e.g. the sinc-type finite-armlength/light-travel-time correction) or
acceleration. For the analytic classes these are exact closed-form
derivatives of `compute_position` (the same formula
`lisaorbits.EqualArmlengthOrbits` differentiates, verified to agree with it
at essentially machine precision in the test suite), not finite-difference
approximations; for `NumericOrbits` they are analytic derivatives of the
interpolating spline, exactly mirroring `lisaorbits.InterpolatedOrbits`.

This module does not change or replace anything in `pycbc.coordinates.space`;
it is purely additive.
"""
import h5py
import numpy as np
from astropy import units
from astropy.constants import au as ASTRONOMICAL_UNIT
from astropy.constants import c as SPEED_OF_LIGHT
from astropy.coordinates import (
    ICRS,
    BarycentricMeanEcliptic,
    get_body_barycentric_posvel,
)
from astropy.time import Time, TimeDelta
from scipy.interpolate import make_interp_spline
from scipy.optimize import fsolve

logger = __import__('logging').getLogger(__name__)


def _parse_oem_file(path):
    """Minimal, dependency-free parser for a CCSDS OEM (Orbit Ephemeris
    Message) v2.0 file, e.g. as published in
    https://github.com/esa/lisa-orbit-files and read by
    `lisaorbits.OEMOrbits` (via the external `oem` package). Only the
    first ephemeris segment is distinguished; a file with more than one
    `META_START`/`META_STOP` block would have its data rows silently
    concatenated (not a concern for the single-segment files tested
    against here).

    Parameters
    ----------
    path : str
        Path to the OEM file.

    Returns
    -------
    meta : dict
        Metadata key/value pairs (e.g. `REF_FRAME`, `CENTER_NAME`,
        `TIME_SYSTEM`), collected from both the file header and the
        META_START/META_STOP block.
    epochs_iso : list of str
        ISO 8601 timestamp string for each data row.
    rows : (N, 6) or (N, 9) ndarray
        Position [km], velocity [km/s], and (if present) acceleration
        [km/s^2] for each data row, in that column order.
    """
    meta = {}
    epochs_iso = []
    values = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line in ('META_START', 'META_STOP'):
                continue
            if line[0].isdigit():
                fields = line.split()
                epochs_iso.append(fields[0])
                values.append([float(x) for x in fields[1:]])
            elif '=' in line:
                key, _, value = line.partition('=')
                meta[key.strip()] = value.strip()
    if not epochs_iso:
        raise ValueError(f'{path}: no data rows found')
    version = meta.get('CCSDS_OEM_VERS')
    if version is not None and version != '2.0':
        raise ValueError(
            f'{path}: unsupported CCSDS_OEM_VERS {version!r} (expected '
            "'2.0')")
    return meta, epochs_iso, np.array(values)


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
        `from_file` re-derives velocities from the position spline exactly
        as this instance does, rather than through an extra layer of
        spline interpolation of re-sampled velocities. Acceleration is
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

    @classmethod
    def from_oem_files(cls, oem_1, oem_2, oem_3, interp_order=5):
        """Load a numerical constellation orbit from three CCSDS OEM
        (Orbit Ephemeris Message) files, one per spacecraft -- the format
        used by ESA's official LISA orbit files
        (https://github.com/esa/lisa-orbit-files), readable by
        `lisaorbits.OEMOrbits`, parsed here with a small self-contained
        parser (no dependency on `lisaorbits` or its `oem` package).

        OEM files give position/velocity in the ICRS/EME2000 (equatorial)
        frame, km and km/s, ISO 8601 timestamps in TCB or TDB; this
        converts to pycbc's convention (SSB BarycentricMeanEcliptic,
        meters, GPS seconds) via `_icrs_to_ecliptic_rotation_matrix`. Any
        acceleration triplet in the file is ignored -- acceleration is
        always derived from the velocity spline, as in `__init__`.

        Only the first ephemeris segment of each file is read (matching
        `lisaorbits.OEMOrbits`); raises `ValueError` if `REF_FRAME` isn't
        `EME2000`, or if the three files don't share identical epochs.

        Parameters
        ----------
        oem_1, oem_2, oem_3 : str
            Paths to the OEM files for spacecraft 1, 2, and 3.
        interp_order : int, optional
            See `__init__`. Default 5.

        Returns
        -------
        NumericOrbits
            A `NumericOrbits` instance, in pycbc's SSB-ecliptic convention,
            built from the three files' contents.
        """
        t_gps = None
        positions_icrs_km = []
        for path in (oem_1, oem_2, oem_3):
            meta, epochs_iso, rows = _parse_oem_file(path)
            if meta.get('REF_FRAME') != 'EME2000':
                raise ValueError(
                    f"{path}: expected REF_FRAME = EME2000, got "
                    f"{meta.get('REF_FRAME')!r}")
            scale = 'tdb' if meta.get('TIME_SYSTEM') == 'TDB' else 'tcb'
            this_t_gps = Time(
                epochs_iso, format='isot', scale=scale).gps
            if t_gps is None:
                t_gps = this_t_gps
            elif not np.array_equal(t_gps, this_t_gps):
                raise ValueError(
                    'input OEM files do not share identical epochs')
            positions_icrs_km.append(rows[:, 0:3])

        positions_icrs_m = np.stack(positions_icrs_km, axis=1) * 1e3  # (N,3,3)
        rotation = _icrs_to_ecliptic_rotation_matrix()
        positions_ecliptic_m = (
            positions_icrs_m.reshape(-1, 3) @ rotation.T
        ).reshape(positions_icrs_m.shape)
        return cls(t_gps, positions_ecliptic_m, interp_order=interp_order)

    @classmethod
    def from_lisaorbits_file(cls, path, interp_order=5):
        """Load a numerical constellation orbit from an orbit file
        generated by `lisaorbits.Orbits.write` (the format used by the
        wider LISA community's simulation tools, e.g. LISA Instrument),
        without requiring `lisaorbits` itself to be installed.

        Positions come from the `tcb/x` dataset, ICRS frame, sampled
        uniformly in TCB time as `t0 + arange(size)*dt` (from the file's
        `t0`/`dt`/`size` attrs, TCB seconds since
        `lisaconstants.LISA_EPOCH_TCB` = 2035-01-01T00:00:00 TCB);
        converted to pycbc's convention the same way as `from_oem_files`.
        As there, `tcb/v`/`tcb/a` are ignored -- velocity/acceleration
        come from the position spline instead.

        Parameters
        ----------
        path : str
            Path to a `lisaorbits`-format orbit HDF5 file.
        interp_order : int, optional
            See `__init__`. Default 5.

        Returns
        -------
        NumericOrbits
            A `NumericOrbits` instance, in pycbc's SSB-ecliptic convention,
            built from the file's contents.
        """
        with h5py.File(path, 'r') as f:
            t0 = float(f.attrs['t0'])
            dt = float(f.attrs['dt'])
            size = int(f.attrs['size'])
            positions_icrs_m = f['tcb/x'][:]

        t_tcb = t0 + np.arange(size) * dt
        epoch = Time('2035-01-01T00:00:00.000', format='isot', scale='tcb')
        t_gps = (epoch + TimeDelta(t_tcb, format='sec')).gps

        rotation = _icrs_to_ecliptic_rotation_matrix()
        positions_ecliptic_m = (
            positions_icrs_m.reshape(-1, 3) @ rotation.T
        ).reshape(positions_icrs_m.shape)
        return cls(t_gps, positions_ecliptic_m, interp_order=interp_order)

    @classmethod
    def from_triangle_dat_files(cls, orbit_dir, sc_labels=('1', '2', '3'),
                                t0=0.0, dt=86400.0, max_rows=None,
                                interp_order=5):
        """Load a numerical constellation orbit from Taiji Data Challenge
        (TDC) Triangle-Simulator-format orbit files: a directory
        containing `SCP{label}.dat`/`SCV{label}.dat` per spacecraft
        (position/velocity, one row per sample, ASCII, AU / AU-per-day,
        ecliptic SSB frame -- Triangle-Simulator's own `Orbit` class
        format, shipped with its `OrbitData/` datasets).

        Unlike OEM, these files carry no absolute timestamps, just a row
        index -- the caller supplies the SSB time of the first row
        (`t0`) and the row spacing (`dt`, default 1 day, matching every
        `OrbitData/` dataset Triangle-Simulator ships).

        Parameters
        ----------
        orbit_dir : str
            Directory containing `SCP{label}.dat`/`SCV{label}.dat` files.
        sc_labels : sequence of str, optional
            Spacecraft file-name labels, in the order they should be
            assigned to 1-indexed spacecraft. Default `('1', '2', '3')`,
            matching Triangle-Simulator's own convention.
        t0 : float, optional
            SSB time [s] (pycbc convention: GPS seconds) of the first row.
            Default 0.0 -- the caller is responsible for supplying the
            correct absolute epoch if one is needed; these files carry no
            timestamp of their own.
        dt : float, optional
            Time spacing between rows [s]. Default 86400.0 (1 day),
            matching every dataset shipped with Triangle-Simulator.
        max_rows : int or None, optional
            Passed to `numpy.loadtxt` to read only the first `max_rows`
            samples (for quick tests on large files). Default None (read
            all rows).
        interp_order : int, optional
            See `__init__`. Default 5.

        Returns
        -------
        NumericOrbits
            A `NumericOrbits` instance built from the three files' contents.
        """
        positions = []
        velocities = []
        n_rows = None
        for label in sc_labels:
            pos = np.loadtxt(
                f'{orbit_dir}/SCP{label}.dat',
                max_rows=max_rows) * ASTRONOMICAL_UNIT.value
            vel = np.loadtxt(
                f'{orbit_dir}/SCV{label}.dat',
                max_rows=max_rows) * ASTRONOMICAL_UNIT.value / 86400.0
            if n_rows is None:
                n_rows = pos.shape[0]
            elif pos.shape[0] != n_rows or vel.shape[0] != n_rows:
                raise ValueError(
                    f'{orbit_dir}: SCP/SCV*.dat files have mismatched '
                    'row counts across spacecraft')
            positions.append(pos)
            velocities.append(vel)

        t_interp = t0 + np.arange(n_rows) * dt
        positions = np.stack(positions, axis=1)  # (N, M, 3)
        velocities = np.stack(velocities, axis=1)
        return cls(t_interp, positions, velocities=velocities,
                   interp_order=interp_order)


_ICRS_TO_ECLIPTIC_ROTATION = None


def _icrs_to_ecliptic_rotation_matrix():
    """The fixed 3x3 rotation matrix from ICRS (equatorial, e.g. EME2000)
    to the BarycentricMeanEcliptic(equinox='J2000') convention this module
    and `pycbc.coordinates.space` use. Both frames are barycentric (no
    origin shift), and the equinox is pinned to a fixed epoch, so this is a
    pure, time-independent axis rotation; it is computed once (by rotating
    the ICRS unit basis vectors) and cached, rather than invoking the full
    astropy frame-transform machinery on every call.
    """
    global _ICRS_TO_ECLIPTIC_ROTATION
    if _ICRS_TO_ECLIPTIC_ROTATION is None:
        basis = np.eye(3)
        icrs = ICRS(x=basis[0] * units.m, y=basis[1] * units.m,
                    z=basis[2] * units.m, representation_type='cartesian')
        ecl = icrs.transform_to(
            BarycentricMeanEcliptic(equinox='J2000')).cartesian
        # transform(e_k) gives the k-th column of the ICRS->ecliptic
        # rotation matrix directly (no transpose): stacking ecl.x/y/z as
        # rows already places ecl.{x,y,z}[k] (the i-th ecliptic coordinate
        # of ICRS basis vector e_k) at matrix position [i, k].
        _ICRS_TO_ECLIPTIC_ROTATION = np.array([
            ecl.x.to(units.m).value,
            ecl.y.to(units.m).value,
            ecl.z.to(units.m).value,
        ])
    return _ICRS_TO_ECLIPTIC_ROTATION


def _real_earth_ecliptic_longitude(t=0.0):
    """Real Earth's ecliptic longitude at SSB time `t` [s], from
    `pycbc.coordinates.space.earth_position_ssb` (astropy-based real
    ephemeris). Imported lazily (not at module level) to avoid a circular
    import: `pycbc.coordinates.space` imports from this module.
    """
    from pycbc.coordinates.space import earth_position_ssb
    return earth_position_ssb(t)[1]


def _real_body_position_velocity(t, body='earth'):
    """Real position and velocity of any astropy-recognized solar-system
    body (e.g. 'earth', 'moon') in the SSB ecliptic frame, at SSB time(s)
    `t` [s] (GPS seconds), from astropy's JPL-based ephemeris (includes
    real perturbations, unlike any closed-form orbit approximation). Uses
    the cached `_icrs_to_ecliptic_rotation_matrix` rather than
    `pycbc.coordinates.space.earth_position_ssb`'s own per-call astropy
    frame transform, so repeated queries stay fast. Imported lazily to
    avoid a circular import.

    Parameters
    ----------
    t : array-like
        SSB time(s) [s] (GPS seconds).
    body : str, optional
        Name of the body, as recognized by
        `astropy.coordinates.get_body_barycentric_posvel`. Default
        `'earth'`.

    Returns
    -------
    position : (N, 3) ndarray
        The body's position in the SSB ecliptic frame [m].
    velocity : (N, 3) ndarray
        The body's velocity in the SSB ecliptic frame [m/s].
    """
    time = Time(np.atleast_1d(t), format='gps')
    pos, vel = get_body_barycentric_posvel(body, time)
    rotation = _icrs_to_ecliptic_rotation_matrix()
    pos_icrs = np.stack([pos.x.to(units.m).value,
                        pos.y.to(units.m).value,
                        pos.z.to(units.m).value], axis=-1)
    vel_icrs = np.stack([vel.x.to(units.m / units.s).value,
                        vel.y.to(units.m / units.s).value,
                        vel.z.to(units.m / units.s).value], axis=-1)
    return pos_icrs @ rotation.T, vel_icrs @ rotation.T


def _real_earth_position_velocity(t):
    """Real Earth position/velocity in the SSB ecliptic frame. See
    `_real_body_position_velocity` (this is a thin wrapper for `'earth'`,
    kept for the existing call sites in this module).
    """
    return _real_body_position_velocity(t, 'earth')


EARTH_ORBIT_ANGULAR_FREQUENCY = 1.99098659277e-7  # [rad/s]
# = 2*pi / (sidereal year in seconds, 31558149.7635456); pinned as a
# literal (not derived via `import lal`) per project convention of
# avoiding a lal dependency in new code.

# Named constants instead of np.deg2rad() calls directly in the default
# constructor arguments below, to avoid a ruff B008 warning.
_TAIJI_LEAD_ANGLE = np.deg2rad(20.0)
_TIANQIN_LAMBDA_S = np.deg2rad(120.5)
_TIANQIN_BETA_S = np.deg2rad(-4.7)
# The theoretical value (not the 1.7e5 km engineering rounding also seen in
# the literature), for consistency with pycbc.psd.analytical_space's own
# TianQin PSD functions, which all default to this same value.
_TIANQIN_ARMLENGTH = np.sqrt(3) * 1e8


def _equal_arm_orbit_position(alpha, armlength, sc):
    """Shared first-order-in-eccentricity Keplerian expansion (Rubbo,
    Cornish & Poujade 2004, Phys. Rev. D 69, 082003) underlying
    `LisaEqualArmOrbit` and `TaijiEqualArmOrbit`: a rigid, circular
    triangular constellation at guiding-center phase `alpha` [rad], with
    the given `armlength` [m].

    `sin(alpha)`/`cos(alpha)` (the only per-element transcendental calls
    needed -- everything else reduces to `beta_n`-dependent scalars via
    the angle-subtraction identity for `cos(alpha - beta_n)`) are
    evaluated once and reused across spacecraft, not recomputed inside the
    loop: for large `alpha` arrays this is the dominant cost, so this
    matters for performance parity with `lisaorbits`, not just style.
    """
    a = ASTRONOMICAL_UNIT.value
    e = armlength / (2 * a * np.sqrt(3))
    sin_a, cos_a = np.sin(alpha), np.cos(alpha)
    sin_a_cos_a = sin_a * cos_a
    sin_a2, cos_a2 = sin_a ** 2, cos_a ** 2
    out = np.empty((len(alpha), len(sc), 3))
    for k, n in enumerate(sc):
        beta_n = (n - 1) * 2 * np.pi / 3.0
        sin_b, cos_b = np.sin(beta_n), np.cos(beta_n)
        out[:, k, 0] = a * cos_a + a * e * (
            sin_a_cos_a * sin_b - (1 + sin_a2) * cos_b)
        out[:, k, 1] = a * sin_a + a * e * (
            sin_a_cos_a * cos_b - (1 + cos_a2) * sin_b)
        out[:, k, 2] = -np.sqrt(3) * a * e * (
            cos_a * cos_b + sin_a * sin_b)  # cos(alpha - beta_n)
    return out


def _equal_arm_orbit_velocity(alpha, omega, armlength, sc):
    """d/dt of `_equal_arm_orbit_position`, given the (constant) guiding-
    center angular frequency `omega` = d(alpha)/dt [rad/s]. Exact analytic
    derivative of the same closed-form position expansion (chain rule
    through `alpha(t)`), not a finite-difference approximation -- the same
    precision-in-principle as `lisaorbits.EqualArmlengthOrbits.
    compute_velocity`, which differentiates the identical formula. See
    `_equal_arm_orbit_position` for the per-spacecraft caching rationale.
    """
    a = ASTRONOMICAL_UNIT.value
    e = armlength / (2 * a * np.sqrt(3))
    sin_a, cos_a = np.sin(alpha), np.cos(alpha)
    sin_a_cos_a = sin_a * cos_a
    cos2_minus_sin2 = cos_a ** 2 - sin_a ** 2
    out = np.empty((len(alpha), len(sc), 3))
    for k, n in enumerate(sc):
        beta_n = (n - 1) * 2 * np.pi / 3.0
        sin_b, cos_b = np.sin(beta_n), np.cos(beta_n)
        out[:, k, 0] = omega * (
            -a * sin_a + a * e * (
                cos2_minus_sin2 * sin_b - 2 * sin_a_cos_a * cos_b))
        out[:, k, 1] = omega * (
            a * cos_a + a * e * (
                cos2_minus_sin2 * cos_b + 2 * sin_a_cos_a * sin_b))
        out[:, k, 2] = omega * (
            np.sqrt(3) * a * e * (
                sin_a * cos_b - cos_a * sin_b))  # sin(alpha - beta_n)
    return out


def _equal_arm_orbit_acceleration(alpha, omega, armlength, sc):
    """d^2/dt^2 of `_equal_arm_orbit_position`, i.e. d/dt of
    `_equal_arm_orbit_velocity`. See `_equal_arm_orbit_position` for the
    precision/performance rationale.
    """
    a = ASTRONOMICAL_UNIT.value
    e = armlength / (2 * a * np.sqrt(3))
    sin_a, cos_a = np.sin(alpha), np.cos(alpha)
    sin_a_cos_a = sin_a * cos_a
    sin_a2, cos_a2 = sin_a ** 2, cos_a ** 2
    omega2 = omega ** 2
    out = np.empty((len(alpha), len(sc), 3))
    for k, n in enumerate(sc):
        beta_n = (n - 1) * 2 * np.pi / 3.0
        sin_b, cos_b = np.sin(beta_n), np.cos(beta_n)
        out[:, k, 0] = omega2 * (
            -a * cos_a - 4 * a * e * (
                sin_a_cos_a * sin_b + (0.5 - sin_a2) * cos_b))
        out[:, k, 1] = omega2 * (
            -a * sin_a - 4 * a * e * (
                sin_a_cos_a * cos_b + (0.5 - cos_a2) * sin_b))
        out[:, k, 2] = omega2 * (
            np.sqrt(3) * a * e * (
                cos_a * cos_b + sin_a * sin_b))  # cos(alpha - beta_n)
    return out


def _kepler_orbit_elements(armlength, semi_major_axis):
    """Eccentricity and inclination for the second-order-in-eccentricity
    equal-arm Kepler constellation: three spacecraft on independent
    two-body Kepler ellipses (same semi-major axis and eccentricity,
    tilted out of the ecliptic by a common inclination, each rotated 120
    degrees in mean anomaly and ascending node from the others), chosen
    so mutual distances stay close to `armlength` to second order in
    eccentricity -- one order better than `_equal_arm_orbit_position`'s
    flat construction.

    Standard "tilted formation" result (same construction as
    `lisaorbits.KeplerianOrbits`); implemented here independently and
    validated against its output to near machine precision.
    """
    alpha = armlength / (2.0 * semi_major_axis)
    delta = 5.0 / 8.0
    nu = np.pi / 3.0 + delta * alpha
    e = np.sqrt(
        1 + 4 * alpha * np.cos(nu) / np.sqrt(3) + 4 * alpha ** 2 / 3) - 1
    tan_i = (alpha * np.sin(nu)
             / (np.sqrt(3) / 2.0 + alpha * np.cos(nu)))
    cos_i = 1.0 / np.sqrt(1 + tan_i ** 2)
    sin_i = tan_i * cos_i
    return e, cos_i, sin_i


def _kepler_eccentric_anomaly(mean_anomaly, e, kepler_order=2):
    """Solve Kepler's equation `psi - e*sin(psi) = mean_anomaly` for the
    eccentric anomaly `psi`: start from the low-eccentricity "equation
    of the center" series (accurate to O(e^4)) and refine with
    `kepler_order` Newton-Raphson iterations. For LISA/Taiji's small
    eccentricities (~1e-2), one iteration already reaches machine
    precision.
    """
    m = mean_anomaly
    psi = (m + (e - e ** 3 / 8.0) * np.sin(m)
           + 0.5 * e ** 2 * np.sin(2 * m)
           + 3.0 / 8.0 * e ** 3 * np.sin(3 * m))
    for _ in range(kepler_order):
        error = psi - e * np.sin(psi) - m
        psi = psi - error / (1 - e * np.cos(psi))
    return psi


def _kepler_orbit_position(alpha, armlength, semi_major_axis, sc,
                           kepler_order=2):
    """Spacecraft position(s) for the second-order-in-eccentricity
    tilted-Kepler-ellipse equal-arm constellation (see
    `_kepler_orbit_elements`), at guiding-center phase `alpha` [rad].
    """
    e, cos_i, sin_i = _kepler_orbit_elements(armlength, semi_major_axis)
    a = semi_major_axis
    out = np.empty((len(alpha), len(sc), 3))
    for k, n in enumerate(sc):
        theta_n = (n - 1) * 2 * np.pi / 3.0
        psi = _kepler_eccentric_anomaly(
            alpha - theta_n, e, kepler_order=kepler_order)
        cos_psi, sin_psi = np.cos(psi), np.sin(psi)
        ref_x = a * cos_i * (cos_psi - e)
        ref_y = a * np.sqrt(1 - e ** 2) * sin_psi
        ref_z = -a * sin_i * (cos_psi - e)
        cos_lam, sin_lam = np.cos(theta_n), np.sin(theta_n)
        out[:, k, 0] = cos_lam * ref_x - sin_lam * ref_y
        out[:, k, 1] = sin_lam * ref_x + cos_lam * ref_y
        out[:, k, 2] = ref_z
    return out


def _kepler_orbit_velocity(alpha, omega, armlength, semi_major_axis, sc,
                           kepler_order=2):
    """d/dt of `_kepler_orbit_position`, given the (constant) guiding-
    center angular frequency `omega` = d(alpha)/dt [rad/s].
    """
    e, cos_i, sin_i = _kepler_orbit_elements(armlength, semi_major_axis)
    a = semi_major_axis
    out = np.empty((len(alpha), len(sc), 3))
    for k, n in enumerate(sc):
        theta_n = (n - 1) * 2 * np.pi / 3.0
        psi = _kepler_eccentric_anomaly(
            alpha - theta_n, e, kepler_order=kepler_order)
        cos_psi, sin_psi = np.cos(psi), np.sin(psi)
        dpsi_dt = omega / (1 - e * cos_psi)
        ref_vx = -a * dpsi_dt * cos_i * sin_psi
        ref_vy = a * dpsi_dt * np.sqrt(1 - e ** 2) * cos_psi
        ref_vz = a * dpsi_dt * sin_i * sin_psi
        cos_lam, sin_lam = np.cos(theta_n), np.sin(theta_n)
        out[:, k, 0] = cos_lam * ref_vx - sin_lam * ref_vy
        out[:, k, 1] = sin_lam * ref_vx + cos_lam * ref_vy
        out[:, k, 2] = ref_vz
    return out


def _kepler_orbit_acceleration(alpha, omega, armlength, semi_major_axis,
                               sc, kepler_order=2):
    """d^2/dt^2 of `_kepler_orbit_position`. For an unperturbed two-body
    Kepler ellipse, the acceleration is simply the Newtonian gravitational
    acceleration towards the central body, `-n^2 a^3 r / |r|^3` (using
    `n = omega`, `a = semi_major_axis`, matching the standard
    :math:`\\ddot r = -GM r/|r|^3` with `GM = n^2 a^3` for this orbit's own
    period) -- this holds regardless of eccentricity or orbit orientation,
    so it is evaluated directly from the position, not by an independent
    per-component derivative.
    """
    pos = _kepler_orbit_position(
        alpha, armlength, semi_major_axis, sc, kepler_order=kepler_order)
    r = np.linalg.norm(pos, axis=-1)
    factor = -(omega ** 2) * (semi_major_axis ** 3) / r ** 3
    return pos * factor[:, :, np.newaxis]


class _LisaGuidingCenter:
    """Mixin providing `_phase(t)` for LISA's guiding-center convention:
    `t0` is added to *time* (`omega*(t+t0)`), a legacy convention kept for
    BBHx compatibility. Combine with `_EqualArmConstellation`/
    `_KeplerConstellation`; the concrete class's own `__init__` must set
    `self.t0`.
    """

    def _phase(self, t):
        return EARTH_ORBIT_ANGULAR_FREQUENCY * (t + self.t0)


class _TaijiGuidingCenter:
    """Mixin providing `_phase(t)` for Taiji's guiding-center convention:
    `lead_angle` is added directly to *phase* (`omega*t + kappa0 +
    lead_angle`) -- unlike LISA's time-offset convention above, the two
    are not directly comparable term-for-term. Combine with
    `_EqualArmConstellation`/`_KeplerConstellation`; the concrete class's
    own `__init__` must set `self.kappa0`/`self.lead_angle`.
    """

    def _phase(self, t):
        return EARTH_ORBIT_ANGULAR_FREQUENCY * t + self.kappa0 \
            + self.lead_angle


class _EqualArmConstellation:
    """Mixin implementing `compute_position`/`compute_velocity`/
    `compute_acceleration` for the rigid, circular, first-order-in-
    eccentricity equal-arm constellation (Rubbo, Cornish & Poujade 2004,
    Phys. Rev. D 69, 082003) -- shared by `LisaEqualArmOrbit` and
    `TaijiEqualArmOrbit`, which differ only in their guiding-center phase
    convention (see `_LisaGuidingCenter`/`_TaijiGuidingCenter`, combined
    in via multiple inheritance) and their `armlength` default. The
    concrete class's own `__init__` must set `self.armlength`.
    """

    def compute_position(self, t, sc=(1, 2, 3)):
        """Spacecraft position(s) at time(s) `t`. See
        `NumericOrbits.compute_position` for the calling convention.

        Returns
        -------
        (N, M, 3) ndarray
            Spacecraft position(s) in the SSB frame [m].
        """
        t = np.atleast_1d(np.asarray(t, dtype=float))
        sc = np.atleast_1d(sc)
        alpha = self._phase(t)
        return _equal_arm_orbit_position(alpha, self.armlength, sc)

    def compute_velocity(self, t, sc=(1, 2, 3)):
        """Spacecraft velocity/ies at time(s) `t`, as the exact analytic
        derivative of `compute_position` (not a finite-difference
        approximation). See `NumericOrbits.compute_position` for the
        calling convention.

        Returns
        -------
        (N, M, 3) ndarray
            Spacecraft velocity/ies in the SSB frame [m/s].
        """
        t = np.atleast_1d(np.asarray(t, dtype=float))
        sc = np.atleast_1d(sc)
        alpha = self._phase(t)
        return _equal_arm_orbit_velocity(
            alpha, EARTH_ORBIT_ANGULAR_FREQUENCY, self.armlength, sc)

    def compute_acceleration(self, t, sc=(1, 2, 3)):
        """Spacecraft acceleration(s) at time(s) `t`, as the exact
        analytic second derivative of `compute_position`. See
        `NumericOrbits.compute_position` for the calling convention.

        Returns
        -------
        (N, M, 3) ndarray
            Spacecraft acceleration(s) in the SSB frame [m/s^2].
        """
        t = np.atleast_1d(np.asarray(t, dtype=float))
        sc = np.atleast_1d(sc)
        alpha = self._phase(t)
        return _equal_arm_orbit_acceleration(
            alpha, EARTH_ORBIT_ANGULAR_FREQUENCY, self.armlength, sc)


class _KeplerConstellation:
    """Mixin implementing `compute_position`/`compute_velocity`/
    `compute_acceleration` for the second-order-in-eccentricity tilted-
    Kepler-ellipse equal-arm constellation (see `_kepler_orbit_elements`)
    -- shared by `LisaKeplerianOrbit` and `TaijiKeplerianOrbit`, same
    combination pattern as `_EqualArmConstellation`. The concrete class's
    own `__init__` must set `self.armlength`/`self.semi_major_axis`/
    `self.kepler_order`.
    """

    def compute_position(self, t, sc=(1, 2, 3)):
        """Spacecraft position(s) at time(s) `t`. See
        `NumericOrbits.compute_position` for the calling convention.

        Returns
        -------
        (N, M, 3) ndarray
            Spacecraft position(s) in the SSB frame [m].
        """
        t = np.atleast_1d(np.asarray(t, dtype=float))
        sc = np.atleast_1d(sc)
        alpha = self._phase(t)
        return _kepler_orbit_position(
            alpha, self.armlength, self.semi_major_axis, sc,
            kepler_order=self.kepler_order)

    def compute_velocity(self, t, sc=(1, 2, 3)):
        """Spacecraft velocity/ies at time(s) `t`. See
        `NumericOrbits.compute_position` for the calling convention.

        Returns
        -------
        (N, M, 3) ndarray
            Spacecraft velocity/ies in the SSB frame [m/s].
        """
        t = np.atleast_1d(np.asarray(t, dtype=float))
        sc = np.atleast_1d(sc)
        alpha = self._phase(t)
        return _kepler_orbit_velocity(
            alpha, EARTH_ORBIT_ANGULAR_FREQUENCY, self.armlength,
            self.semi_major_axis, sc, kepler_order=self.kepler_order)

    def compute_acceleration(self, t, sc=(1, 2, 3)):
        """Spacecraft acceleration(s) at time(s) `t`. See
        `NumericOrbits.compute_position` for the calling convention.

        Returns
        -------
        (N, M, 3) ndarray
            Spacecraft acceleration(s) in the SSB frame [m/s^2].
        """
        t = np.atleast_1d(np.asarray(t, dtype=float))
        sc = np.atleast_1d(sc)
        alpha = self._phase(t)
        return _kepler_orbit_acceleration(
            alpha, EARTH_ORBIT_ANGULAR_FREQUENCY, self.armlength,
            self.semi_major_axis, sc, kepler_order=self.kepler_order)


class LisaEqualArmOrbit(_LisaGuidingCenter, _EqualArmConstellation):
    """Idealized LISA heliocentric orbit: the same rigid, circular
    (first-order-in-eccentricity) equal-arm triangle already used by
    `pycbc.coordinates.space.lisa_position_ssb`/
    `rotation_matrix_ssb_to_lisa` (the `orbit=None` default of
    `ssb_to_lisa`/`lisa_to_ssb`/etc.), exposed as an explicit
    `OrbitProvider` -- e.g. to pass to `constellation_frame`, or as a
    baseline against a numeric or other-mission orbit.
    `LisaEqualArmOrbit()` with no arguments reproduces the
    existing `orbit=None` behavior exactly.

    Unlike `TaijiEqualArmOrbit`/`TianQinAnalyticOrbit`, `t0` does *not*
    default to a real-Earth-anchored epoch -- it's the pre-existing
    constant tuned for BBHx compatibility, and changing that default
    here would silently disagree with `pycbc.coordinates.space`'s own
    default. Pass a different `t0` for a different reference epoch.

    Parameters
    ----------
    armlength : float, optional
        Constellation arm length [m]. Default 2.5e9 (design value).
    t0 : float or None, optional
        Reference time offset [s], with the same meaning as in
        `pycbc.coordinates.space.lisa_position_ssb`. Default None, which
        uses `pycbc.coordinates.space.TIME_OFFSET_20_DEGREES`.
    """

    def __init__(self, armlength=2.5e9, t0=None):
        self.armlength = float(armlength)
        if t0 is None:
            from pycbc.coordinates.space import TIME_OFFSET_20_DEGREES
            t0 = TIME_OFFSET_20_DEGREES
        self.t0 = float(t0)


class LisaKeplerianOrbit(_LisaGuidingCenter, _KeplerConstellation):
    """LISA heliocentric orbit built from a genuine two-body Kepler
    ellipse per spacecraft, rather than `LisaEqualArmOrbit`'s flat,
    first-order-in-eccentricity expansion: three spacecraft on
    independent, common-inclination, common-eccentricity ellipses, each
    rotated 120 degrees from the others in mean anomaly and ascending
    node, with eccentricity/inclination from a standard "tilted
    formation" equal-arm construction (see `_kepler_orbit_elements`),
    accurate to second order in eccentricity. Independently implemented
    and validated to near machine precision against
    `lisaorbits.KeplerianOrbits`, not ported from it.

    Not strictly "more accurate" than `LisaEqualArmOrbit`: that first-
    order construction happens to keep LISA's arm length essentially
    exactly constant over a year, while this true-Kepler construction has
    a real ~0.2% arm-length variation (an inherent property of the
    construction, not a numerical error). Use this to match reference
    data built with a true-Kepler-ellipse convention (e.g. some LDC/TDC
    products); use `LisaEqualArmOrbit` otherwise (it's also this
    module's pre-existing default).

    Parameters
    ----------
    armlength : float, optional
        Constellation arm length [m]. Default 2.5e9 (design value).
    semi_major_axis : float, optional
        Guiding-center semi-major axis [m]. Default 1 AU.
    t0 : float or None, optional
        Reference time offset [s], with the same meaning as in
        `LisaEqualArmOrbit`. Default None, which uses
        `pycbc.coordinates.space.TIME_OFFSET_20_DEGREES`.
    kepler_order : int, optional
        Number of Newton-Raphson iterations used to solve Kepler's
        equation for the eccentric anomaly. Default 2 (converges to
        machine precision for LISA's eccentricity well within this).
    """

    def __init__(self, armlength=2.5e9, semi_major_axis=None, t0=None,
                kepler_order=2):
        self.armlength = float(armlength)
        self.semi_major_axis = (
            ASTRONOMICAL_UNIT.value if semi_major_axis is None
            else float(semi_major_axis))
        if t0 is None:
            from pycbc.coordinates.space import TIME_OFFSET_20_DEGREES
            t0 = TIME_OFFSET_20_DEGREES
        self.t0 = float(t0)
        self.kepler_order = int(kepler_order)


class TaijiEqualArmOrbit(_TaijiGuidingCenter, _EqualArmConstellation):
    """Idealized Taiji heliocentric orbit: a rigid, circular (first-order-
    in-eccentricity) equal-arm triangle, per Rubbo, Cornish & Poujade
    2004 (Phys. Rev. D 69, 082003) -- the same functional form as the
    LISA orbit in `pycbc.coordinates.space.lisa_position_ssb`/
    `rotation_matrix_ssb_to_lisa`, but leading the Earth-like guiding
    center by `lead_angle` (design value 20 degrees) instead of trailing
    it, with Taiji's own arm length.

    For prototyping (e.g. single-link response and TDI work) ahead of an
    official numerical orbit product -- not a substitute for real
    mission ephemeris; use `NumericOrbits.from_file` once one exists.

    Parameters
    ----------
    armlength : float, optional
        Constellation arm length [m]. Default 3.0e9 (design value).
    lead_angle : float, optional
        Angle by which the constellation leads the Earth-like guiding
        center [rad] -- added directly to orbital phase, so positive means
        ahead of the guiding center (unlike LISA's `t0`, which is added to
        *time* to produce a trailing offset; the two are not directly
        comparable term-for-term). Default ``deg2rad(20)`` (design value).
    kappa0 : float or None, optional
        Reference ecliptic longitude of the Earth-like guiding center at
        `t=0` [rad], before `lead_angle` is added. Default None, which
        anchors it to the real Earth's ecliptic longitude at SSB time 0
        (via `pycbc.coordinates.space.earth_position_ssb`), so that
        `TaijiEqualArmOrbit()` with no arguments is roughly realistic
        "today". Pass an explicit value for an arbitrary or
        scenario-specific reference epoch instead.
    """

    def __init__(self, armlength=3.0e9, lead_angle=_TAIJI_LEAD_ANGLE,
                kappa0=None):
        self.armlength = float(armlength)
        self.lead_angle = float(lead_angle)
        if kappa0 is None:
            kappa0 = _real_earth_ecliptic_longitude(0.0)
        self.kappa0 = float(kappa0)


class TaijiKeplerianOrbit(_TaijiGuidingCenter, _KeplerConstellation):
    """Taiji heliocentric orbit, accurate to second order in eccentricity
    (rather than `TaijiEqualArmOrbit`'s first order) -- see
    `LisaKeplerianOrbit` for the underlying construction. Leads the
    Earth-like guiding center by `lead_angle` (design value 20 degrees),
    with Taiji's own arm length.

    Parameters
    ----------
    armlength : float, optional
        Constellation arm length [m]. Default 3.0e9 (design value).
    semi_major_axis : float, optional
        Guiding-center semi-major axis [m]. Default 1 AU.
    lead_angle : float, optional
        Angle by which the constellation leads the Earth-like guiding
        center [rad] -- added directly to orbital phase, so positive means
        ahead of the guiding center (unlike LISA's `t0`, which is added to
        *time* to produce a trailing offset; the two are not directly
        comparable term-for-term). Default ``deg2rad(20)`` (design value).
    kappa0 : float or None, optional
        Reference ecliptic longitude of the Earth-like guiding center at
        `t=0` [rad], before `lead_angle` is added. Default None, which
        anchors it to the real Earth's ecliptic longitude at SSB time 0,
        as in `TaijiEqualArmOrbit`.
    kepler_order : int, optional
        See `LisaKeplerianOrbit`. Default 2.
    """

    def __init__(self, armlength=3.0e9, semi_major_axis=None,
                lead_angle=_TAIJI_LEAD_ANGLE, kappa0=None, kepler_order=2):
        self.armlength = float(armlength)
        self.semi_major_axis = (
            ASTRONOMICAL_UNIT.value if semi_major_axis is None
            else float(semi_major_axis))
        self.lead_angle = float(lead_angle)
        if kappa0 is None:
            kappa0 = _real_earth_ecliptic_longitude(0.0)
        self.kappa0 = float(kappa0)
        self.kepler_order = int(kepler_order)


class TianQinAnalyticOrbit:
    """Idealized analytic TianQin geocentric orbit: a fast, rigidly-
    rotating triangular constellation whose plane is fixed in inertial
    space (pointing towards the calibration source RX J0806.3+1527),
    around a guiding center coincident with the Earth, per Hu et al 2018
    (Class. Quantum Grav. 35, 095008). Unlike LISA/Taiji, the triangle's
    shape and orientation are held fixed rather than letting each
    spacecraft flex on its own ellipse.

    `guiding_center` picks how the Earth-like guiding center moves:
    `'circular'` (default) is a pure circular heliocentric orbit,
    neglecting Earth's ~1.7% eccentricity; `'real_earth'` uses astropy's
    JPL ephemeris instead, more accurate but needing a lookup (and, for
    acceleration, a finite difference, since astropy has no direct
    acceleration) per call rather than a closed form.

    Not a substitute for real mission ephemeris, same caveat as
    `TaijiEqualArmOrbit`.

    Parameters
    ----------
    armlength : float, optional
        Constellation arm length [m]. Default `sqrt(3) * 1e8` (the
        theoretical design value; matches `pycbc.psd.analytical_space`'s
        TianQin PSD functions).
    lambda_s, beta_s : float, optional
        Ecliptic longitude/latitude of the calibration source RX
        J0806.3+1527 [rad], which fixes the constellation plane's
        orientation. Defaults are the design values (120.5 deg, -4.7 deg).
    rotation_period : float, optional
        Constellation self-rotation period [s]. Default 3.65 days (design
        value).
    initial_orbit_phase : float, optional
        Initial fast-rotation phase of the first spacecraft [rad].
        Default 0 -- unlike `kappa0` below, there is no "real" reference
        for this (it is an arbitrary, mission-dependent choice).
    kappa0 : float or None, optional
        Reference ecliptic longitude of the Earth-like guiding center at
        `t=0` [rad], used only when `guiding_center='circular'`. Default
        None, which anchors it to the real Earth's ecliptic longitude at
        SSB time 0 (via `pycbc.coordinates.space.earth_position_ssb`), so
        that `TianQinAnalyticOrbit()` with no arguments is roughly
        realistic "today". Pass an explicit value for an arbitrary or
        scenario-specific reference epoch instead.
    guiding_center : {'circular', 'real_earth'}, optional
        See above. Default `'circular'` (unchanged default behavior).
    """

    def __init__(self, armlength=_TIANQIN_ARMLENGTH, lambda_s=_TIANQIN_LAMBDA_S,
                beta_s=_TIANQIN_BETA_S, rotation_period=3.65 * 86400.0,
                initial_orbit_phase=0.0, kappa0=None,
                guiding_center='circular'):
        self.armlength = float(armlength)
        self.lambda_s = float(lambda_s)
        self.beta_s = float(beta_s)
        self.rotation_period = float(rotation_period)
        self.initial_orbit_phase = float(initial_orbit_phase)
        if guiding_center not in ('circular', 'real_earth'):
            raise ValueError(
                "guiding_center must be 'circular' or 'real_earth', got "
                f"{guiding_center!r}")
        self.guiding_center = guiding_center
        if kappa0 is None:
            kappa0 = _real_earth_ecliptic_longitude(0.0)
        self.kappa0 = float(kappa0)

    def _guiding_center_position(self, t):
        if self.guiding_center == 'real_earth':
            pos, _ = _real_earth_position_velocity(t)
            return pos
        a = ASTRONOMICAL_UNIT.value
        alpha_earth = EARTH_ORBIT_ANGULAR_FREQUENCY * t + self.kappa0
        return a * np.stack([
            np.cos(alpha_earth), np.sin(alpha_earth),
            np.zeros_like(t)], axis=-1)

    def _guiding_center_position_velocity(self, t):
        if self.guiding_center == 'real_earth':
            return _real_earth_position_velocity(t)
        a = ASTRONOMICAL_UNIT.value
        omega_e = EARTH_ORBIT_ANGULAR_FREQUENCY
        alpha_earth = omega_e * t + self.kappa0
        pos = a * np.stack([
            np.cos(alpha_earth), np.sin(alpha_earth),
            np.zeros_like(t)], axis=-1)
        vel = a * omega_e * np.stack([
            -np.sin(alpha_earth), np.cos(alpha_earth),
            np.zeros_like(t)], axis=-1)
        return pos, vel

    def _guiding_center_acceleration(self, t):
        if self.guiding_center == 'real_earth':
            # astropy does not expose Earth's acceleration directly;
            # a small central finite difference of its (real, ephemeris-
            # based) velocity is far more accurate here than falling back
            # to an idealized two-body -GM_sun*r/|r|^3 approximation,
            # which would reintroduce the very approximation this mode
            # exists to avoid.
            dt = 1.0
            _, vel_plus = _real_earth_position_velocity(t + dt)
            _, vel_minus = _real_earth_position_velocity(t - dt)
            return (vel_plus - vel_minus) / (2 * dt)
        a = ASTRONOMICAL_UNIT.value
        omega_e = EARTH_ORBIT_ANGULAR_FREQUENCY
        alpha_earth = omega_e * t + self.kappa0
        return -a * omega_e ** 2 * np.stack([
            np.cos(alpha_earth), np.sin(alpha_earth),
            np.zeros_like(t)], axis=-1)

    def compute_position(self, t, sc=(1, 2, 3)):
        """Spacecraft position(s) at time(s) `t`. See
        `NumericOrbits.compute_position` for the calling convention.

        Returns
        -------
        (N, M, 3) ndarray
            Spacecraft position(s) in the SSB frame [m].
        """
        t = np.atleast_1d(np.asarray(t, dtype=float))
        sc = np.atleast_1d(sc)
        earth = self._guiding_center_position(t)
        f_sc = 1.0 / self.rotation_period
        out = np.empty((len(t), len(sc), 3))
        for k, n in enumerate(sc):
            kappa_n = 2 * np.pi / 3.0 * (n - 1) + self.initial_orbit_phase
            phase = 2 * np.pi * f_sc * t + kappa_n
            xn = (self.armlength / np.sqrt(3)) * (
                np.sin(self.beta_s) * np.cos(self.lambda_s)
                * np.sin(phase)
                + np.sin(self.lambda_s) * np.cos(phase))
            yn = (self.armlength / np.sqrt(3)) * (
                np.sin(self.beta_s) * np.sin(self.lambda_s)
                * np.sin(phase)
                - np.cos(self.lambda_s) * np.cos(phase))
            zn = -(self.armlength / np.sqrt(3)) \
                * np.cos(self.beta_s) * np.sin(phase)
            out[:, k, 0] = earth[:, 0] + xn
            out[:, k, 1] = earth[:, 1] + yn
            out[:, k, 2] = earth[:, 2] + zn
        return out

    def compute_velocity(self, t, sc=(1, 2, 3)):
        """Spacecraft velocity/ies at time(s) `t`, exact (not a finite-
        difference approximation): with `guiding_center='circular'`, both
        the guiding center and the fast-rotation term have constant
        angular frequency, so each is a simple chain-rule derivative of
        `compute_position`; with `'real_earth'`, the guiding-center
        velocity comes directly from the ephemeris instead. See
        `NumericOrbits.compute_position` for the calling convention.

        Returns
        -------
        (N, M, 3) ndarray
            Spacecraft velocity/ies in the SSB frame [m/s].
        """
        t = np.atleast_1d(np.asarray(t, dtype=float))
        sc = np.atleast_1d(sc)
        _, earth_vel = self._guiding_center_position_velocity(t)
        omega_sc = 2 * np.pi / self.rotation_period
        out = np.empty((len(t), len(sc), 3))
        for k, n in enumerate(sc):
            kappa_n = 2 * np.pi / 3.0 * (n - 1) + self.initial_orbit_phase
            phase = omega_sc * t + kappa_n
            vxn = omega_sc * (self.armlength / np.sqrt(3)) * (
                np.sin(self.beta_s) * np.cos(self.lambda_s)
                * np.cos(phase)
                - np.sin(self.lambda_s) * np.sin(phase))
            vyn = omega_sc * (self.armlength / np.sqrt(3)) * (
                np.sin(self.beta_s) * np.sin(self.lambda_s)
                * np.cos(phase)
                + np.cos(self.lambda_s) * np.sin(phase))
            vzn = -omega_sc * (self.armlength / np.sqrt(3)) \
                * np.cos(self.beta_s) * np.cos(phase)
            out[:, k, 0] = earth_vel[:, 0] + vxn
            out[:, k, 1] = earth_vel[:, 1] + vyn
            out[:, k, 2] = earth_vel[:, 2] + vzn
        return out

    def compute_acceleration(self, t, sc=(1, 2, 3)):
        """Spacecraft acceleration(s) at time(s) `t`, exact: with
        `guiding_center='circular'`, both terms are uniform circular
        motion, so each term's acceleration is -omega^2 times that
        term's position; with `'real_earth'`, the guiding-center
        contribution is a finite difference of the real ephemeris
        velocity instead (see `_guiding_center_acceleration`). See
        `NumericOrbits.compute_position` for the calling convention.

        Returns
        -------
        (N, M, 3) ndarray
            Spacecraft acceleration(s) in the SSB frame [m/s^2].
        """
        t = np.atleast_1d(np.asarray(t, dtype=float))
        sc = np.atleast_1d(sc)
        earth_acc = self._guiding_center_acceleration(t)
        omega_sc = 2 * np.pi / self.rotation_period
        out = np.empty((len(t), len(sc), 3))
        for k, n in enumerate(sc):
            kappa_n = 2 * np.pi / 3.0 * (n - 1) + self.initial_orbit_phase
            phase = omega_sc * t + kappa_n
            xn = (self.armlength / np.sqrt(3)) * (
                np.sin(self.beta_s) * np.cos(self.lambda_s)
                * np.sin(phase)
                + np.sin(self.lambda_s) * np.cos(phase))
            yn = (self.armlength / np.sqrt(3)) * (
                np.sin(self.beta_s) * np.sin(self.lambda_s)
                * np.sin(phase)
                - np.cos(self.lambda_s) * np.cos(phase))
            zn = -(self.armlength / np.sqrt(3)) \
                * np.cos(self.beta_s) * np.sin(phase)
            out[:, k, 0] = earth_acc[:, 0] - omega_sc ** 2 * xn
            out[:, k, 1] = earth_acc[:, 1] - omega_sc ** 2 * yn
            out[:, k, 2] = earth_acc[:, 2] - omega_sc ** 2 * zn
        return out


def constellation_frame(t, orbit, sc=(1, 2, 3)):
    """Compute the instantaneous constellation centroid and the rotation
    matrix from the SSB frame to a frame comoving with the constellation,
    directly from the spacecraft positions supplied by `orbit`.

    Generalizes `rotation_matrix_ssb_to_lisa`/`lisa_position_ssb` in
    `pycbc.coordinates.space` (fixed circular LISA orbit) to any
    3-spacecraft triangular constellation given analytic or numerical
    positions -- e.g. LISA, Taiji, or TianQin.

    For each time sample, the instantaneous frame is built from the three
    spacecraft positions as:

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
    """Shared `fsolve` pattern behind `t_lisa_from_ssb`/`t_ssb_from_t_lisa`,
    `t_geo_from_ssb`/`t_ssb_from_t_geo` (in `pycbc.coordinates.space`), and
    `t_detector_from_ssb`/`t_ssb_from_t_detector` below: solve
    `t - t_other - dot(k, position(t)) / c = 0` for whichever side is
    unknown, given a light-travel-time relationship between the SSB origin
    and a moving frame's origin.

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
        True solves for the frame's own time given `t_known = t_ssb` (the
        frame's position depends on the unknown itself, so `position_fn`
        is re-evaluated inside the root-solve). False solves for `t_ssb`
        given `t_known` = the frame's own time (the frame's position is
        evaluated once, at the known time, since it doesn't depend on the
        unknown `t_ssb`).

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
    'NumericOrbits',
    'LisaEqualArmOrbit',
    'LisaKeplerianOrbit',
    'TaijiEqualArmOrbit',
    'TaijiKeplerianOrbit',
    'TianQinAnalyticOrbit',
    'constellation_frame',
    't_detector_from_ssb',
    't_ssb_from_t_detector',
]
