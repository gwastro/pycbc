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
import h5py
import numpy as np
from scipy.interpolate import make_interp_spline
from scipy.optimize import fsolve
from astropy.constants import c as SPEED_OF_LIGHT
from astropy.constants import au as ASTRONOMICAL_UNIT

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

    @classmethod
    def from_file(cls, path, group=None, interp_order=5):
        """Load a numerical constellation orbit from an HDF5 file.

        The file (or the given group within it) must contain datasets `t`
        (shape (N,), SSB time [s]) and `positions` (shape (N, M, 3), SSB-
        frame spacecraft positions [m]). An optional `velocities` dataset
        (same shape as `positions`) is used if present; otherwise
        velocities are computed as the analytic derivative of the position
        spline, as in `__init__`. This is the file format any of LISA,
        Taiji, TianQin (or a numerically-propagated orbit for any of them)
        can be supplied in, e.g. from a PE config's `orbit-file` option
        (see `pycbc.transforms`).

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
        from astropy import units
        from astropy.coordinates import ICRS, BarycentricMeanEcliptic
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


class ICRSOrbitAdapter:
    """Wrap an `OrbitProvider` given in the ICRS (equatorial) frame so it
    can be used as a drop-in `OrbitProvider` in pycbc's own
    BarycentricMeanEcliptic convention.

    Real spacecraft ephemerides are most often published in an equatorial
    frame (e.g. CCSDS OEM files use EME2000, which `lisaorbits.OEMOrbits`
    reads without further rotation, so any `lisaorbits.Orbits` instance
    built from such files reports positions in ICRS). This module and
    `pycbc.coordinates.space` assume the SSB ecliptic frame throughout, so
    such an orbit cannot be passed to `constellation_frame`/`link_vector`/
    etc. directly without producing silently-wrong frame-dependent
    quantities (e.g. `constellation_frame`'s rotation matrix and centroid
    components) -- even though frame-independent quantities like arm
    length happen to come out correct either way.

    This class only wraps `compute_position`/`compute_velocity`; both are
    rotated with the same fixed matrix, since the ICRS -> ecliptic
    transform (fixed equinox, no origin shift) commutes with time
    differentiation.

    Parameters
    ----------
    orbit : OrbitProvider
        Any object exposing `compute_position(t, sc)` (and, optionally,
        `compute_velocity(t, sc)`) with positions/velocities in the ICRS
        frame -- for example a real `lisaorbits.OEMOrbits` instance built
        from official ESA lisa-orbit-files data.
    """
    def __init__(self, orbit):
        self._orbit = orbit
        self._rotation = _icrs_to_ecliptic_rotation_matrix()

    def compute_position(self, t, sc=(1, 2, 3)):
        return self._orbit.compute_position(t, sc) @ self._rotation.T

    def compute_velocity(self, t, sc=(1, 2, 3)):
        return self._orbit.compute_velocity(t, sc) @ self._rotation.T


def _real_earth_ecliptic_longitude(t=0.0):
    """Real Earth's ecliptic longitude at SSB time `t` [s], from
    `pycbc.coordinates.space.earth_position_ssb` (astropy-based real
    ephemeris). Imported lazily (not at module level) to avoid a circular
    import: `pycbc.coordinates.space` imports from this module.
    """
    from pycbc.coordinates.space import earth_position_ssb
    return earth_position_ssb(t)[1]


EARTH_ORBIT_ANGULAR_FREQUENCY = 1.99098659277e-7  # [rad/s], ~1 sidereal year


def _equal_arm_orbit_position(alpha, armlength, sc):
    """Shared first-order-in-eccentricity Keplerian expansion (Rubbo,
    Cornish & Poujade 2004, Phys. Rev. D 69, 082003) underlying
    `LisaAnalyticOrbit` and `TaijiAnalyticOrbit`: a rigid, circular
    triangular constellation at guiding-center phase `alpha` [rad], with
    the given `armlength` [m].
    """
    a = ASTRONOMICAL_UNIT.value
    e = armlength / (2 * a * np.sqrt(3))
    out = np.empty((len(alpha), len(sc), 3))
    for k, n in enumerate(sc):
        beta_n = (n - 1) * 2 * np.pi / 3.0
        out[:, k, 0] = a * np.cos(alpha) + a * e * (
            np.sin(alpha) * np.cos(alpha) * np.sin(beta_n)
            - (1 + np.sin(alpha) ** 2) * np.cos(beta_n))
        out[:, k, 1] = a * np.sin(alpha) + a * e * (
            np.sin(alpha) * np.cos(alpha) * np.cos(beta_n)
            - (1 + np.cos(alpha) ** 2) * np.sin(beta_n))
        out[:, k, 2] = -np.sqrt(3) * a * e * np.cos(alpha - beta_n)
    return out


class LisaAnalyticOrbit:
    """Idealized analytic LISA heliocentric orbit: the same rigid,
    circular (first order in eccentricity) triangular constellation
    already used by `pycbc.coordinates.space.lisa_position_ssb`/
    `rotation_matrix_ssb_to_lisa` (the `orbit=None` default of
    `ssb_to_lisa`/`lisa_to_ssb`/etc.), exposed here as an explicit
    `OrbitProvider` object -- e.g. to pass directly to
    `constellation_frame`/`link_vector`, or as a drop-in comparison
    baseline against a numeric or other-mission orbit.

    `LisaAnalyticOrbit()` with no arguments reproduces the existing
    `orbit=None` default behavior exactly (same `t0`, same formula).
    Unlike `TaijiAnalyticOrbit`/`TianQinAnalyticOrbit`, it does *not*
    default to a real-Earth-anchored reference epoch: `t0` is already a
    real, pre-existing constant tuned for compatibility with the BBHx
    waveform plugin, and changing that default here would silently
    disagree with `pycbc.coordinates.space`'s own default. Pass a
    different `t0` for a different reference epoch.

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
        alpha = EARTH_ORBIT_ANGULAR_FREQUENCY * (t + self.t0)
        return _equal_arm_orbit_position(alpha, self.armlength, sc)


class TaijiAnalyticOrbit:
    """Idealized analytic Taiji heliocentric orbit: a rigid, circular
    (first order in eccentricity) triangular constellation, per Rubbo,
    Cornish & Poujade 2004 (Phys. Rev. D 69, 082003) -- the same
    functional form underlying the LISA orbit in
    `pycbc.coordinates.space.lisa_position_ssb`/`rotation_matrix_ssb_to_lisa`
    -- leading the Earth-like guiding center by `lead_angle` (design value
    20 degrees) instead of trailing it, with Taiji's own arm length.

    This is an idealized reference orbit intended for prototyping and
    methods development (e.g. single-link response and TDI work) ahead of
    an official numerical orbit product. It is not a substitute for real
    mission ephemeris in science-quality analysis -- use
    `NumericOrbits.from_file` with an official orbit file for that once
    one exists.

    Parameters
    ----------
    armlength : float, optional
        Constellation arm length [m]. Default 3.0e9 (design value).
    lead_angle : float, optional
        Angle by which the constellation leads the Earth-like guiding
        center [rad]. Default ``deg2rad(20)`` (design value).
    kappa0 : float or None, optional
        Reference ecliptic longitude of the Earth-like guiding center at
        `t=0` [rad], before `lead_angle` is added. Default None, which
        anchors it to the real Earth's ecliptic longitude at SSB time 0
        (via `pycbc.coordinates.space.earth_position_ssb`), so that
        `TaijiAnalyticOrbit()` with no arguments is roughly realistic
        "today". Pass an explicit value for an arbitrary or
        scenario-specific reference epoch instead.
    """

    def __init__(self, armlength=3.0e9, lead_angle=np.deg2rad(20.0),
                kappa0=None):
        self.armlength = float(armlength)
        self.lead_angle = float(lead_angle)
        if kappa0 is None:
            kappa0 = _real_earth_ecliptic_longitude(0.0)
        self.kappa0 = float(kappa0)

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
        alpha = EARTH_ORBIT_ANGULAR_FREQUENCY * t + self.kappa0 \
            + self.lead_angle
        return _equal_arm_orbit_position(alpha, self.armlength, sc)


class TianQinAnalyticOrbit:
    """Idealized analytic TianQin geocentric orbit: a fast, rigidly-
    rotating triangular constellation whose plane is fixed in inertial
    space (pointing towards the calibration source RX J0806.3+1527),
    around a guiding center coincident with the Earth, per Hu et al 2018
    (Class. Quantum Grav. 35, 095008). The guiding center itself is
    approximated here as a pure circular heliocentric orbit (Earth's real
    ~1.7% orbital eccentricity is neglected).

    This is an idealized reference orbit intended for prototyping and
    methods development ahead of an official numerical orbit product. It
    is not a substitute for real mission ephemeris in science-quality
    analysis -- use `NumericOrbits.from_file` with an official orbit file
    for that once one exists.

    Parameters
    ----------
    armlength : float, optional
        Constellation arm length [m]. Default 1.7e8 (design value).
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
        `t=0` [rad]. Default None, which anchors it to the real Earth's
        ecliptic longitude at SSB time 0 (via
        `pycbc.coordinates.space.earth_position_ssb`), so that
        `TianQinAnalyticOrbit()` with no arguments is roughly realistic
        "today". Pass an explicit value for an arbitrary or
        scenario-specific reference epoch instead.
    """

    def __init__(self, armlength=1.7e8, lambda_s=np.deg2rad(120.5),
                beta_s=np.deg2rad(-4.7), rotation_period=3.65 * 86400.0,
                initial_orbit_phase=0.0, kappa0=None):
        self.armlength = float(armlength)
        self.lambda_s = float(lambda_s)
        self.beta_s = float(beta_s)
        self.rotation_period = float(rotation_period)
        self.initial_orbit_phase = float(initial_orbit_phase)
        if kappa0 is None:
            kappa0 = _real_earth_ecliptic_longitude(0.0)
        self.kappa0 = float(kappa0)

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
        a = ASTRONOMICAL_UNIT.value
        alpha_earth = EARTH_ORBIT_ANGULAR_FREQUENCY * t + self.kappa0
        earth = a * np.stack([
            np.cos(alpha_earth), np.sin(alpha_earth),
            np.zeros_like(t)], axis=-1)
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


def link_vector(t, orbit, sc_emitter, sc_receiver):
    """Compute the instantaneous unit vector and arm length for a single
    laser link between two spacecraft -- the geometric input needed by
    single-link response formulas, for any constellation orbit (LISA,
    Taiji, TianQin, analytic or numerical).

    The unit vector points from the emitter spacecraft to the receiver
    spacecraft, both evaluated at the same time `t` (i.e. this neglects the
    light travel time between emission and reception). This matches the
    convention of `lisaorbits.Orbits.compute_unit_vector`, except that
    function retards the emitter position by the link's actual light travel
    time; the same-time approximation here is what is needed for the
    spatial projection factors (e.g. n_hat . k_hat, or (n_hat x n_hat):P)
    that enter single-link response formulas, and is accurate to the same
    order those formulas already are. Applications that need the true
    retarded emission time (e.g. exact TDI time-delay operators) can solve
    for it explicitly with the same self-consistent approach
    `t_detector_from_ssb` uses for the constellation centroid.

    Parameters
    ----------
    t : (N,) array-like
        SSB time(s) [s].
    orbit : OrbitProvider
        Any object exposing `compute_position(t, sc)`, see
        `constellation_frame`.
    sc_emitter : int or (M,) array-like
        1-indexed spacecraft label(s) of the emitter.
    sc_receiver : int or (M,) array-like
        1-indexed spacecraft label(s) of the receiver, same shape as
        `sc_emitter`. `link_vector(t, orbit, i, j)` and
        `link_vector(t, orbit, j, i)` give antiparallel unit vectors and
        the same arm length.

    Returns
    -------
    unit_vector : (N, M, 3) ndarray
        Unit vector(s) pointing from emitter to receiver.
    arm_length : (N, M) ndarray
        Instantaneous distance(s) between emitter and receiver [m].
    """
    t = np.atleast_1d(np.asarray(t, dtype=float))
    sc_emitter = np.atleast_1d(np.asarray(sc_emitter))
    sc_receiver = np.atleast_1d(np.asarray(sc_receiver))
    if sc_emitter.shape != sc_receiver.shape:
        raise ValueError(
            'sc_emitter and sc_receiver must have the same shape, got '
            f'{sc_emitter.shape} and {sc_receiver.shape}')

    pos_emitter = orbit.compute_position(t, sc_emitter)  # (N, M, 3)
    pos_receiver = orbit.compute_position(t, sc_receiver)  # (N, M, 3)
    delta = pos_receiver - pos_emitter
    arm_length = np.linalg.norm(delta, axis=-1)  # (N, M)
    unit_vector = delta / arm_length[..., np.newaxis]
    return unit_vector, arm_length


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
    'LisaAnalyticOrbit',
    'TaijiAnalyticOrbit',
    'TianQinAnalyticOrbit',
    'ICRSOrbitAdapter',
    'constellation_frame',
    'link_vector',
    't_detector_from_ssb',
    't_ssb_from_t_detector',
]
