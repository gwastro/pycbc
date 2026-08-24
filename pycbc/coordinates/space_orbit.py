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
Analytic, closed-form constellation orbits for space-borne detectors.

`LisaEqualArmOrbit` and `TaijiEqualArmOrbit` are rigid, circular (first
order in eccentricity) equal-arm triangular constellations following the
expansion of Rubbo, Cornish & Poujade 2004 (Phys. Rev. D 69, 082003) -- the
same functional form as the LISA orbit hard-coded in
`pycbc.coordinates.space.lisa_position_ssb`/`rotation_matrix_ssb_to_lisa`,
here exposed as explicit, parameterised orbit providers.

Each class implements `compute_position(t, sc)`, `compute_velocity(t, sc)`
and `compute_acceleration(t, sc)` with the same calling convention and
return shapes as `lisaorbits.Orbits`, so they can be passed anywhere an
orbit provider is expected. Velocity and acceleration are exact analytic
derivatives of `compute_position`, not finite differences.

These are idealised reference orbits for prototyping ahead of an official
numeric orbit product; they are not a substitute for real mission
ephemerides.

This module does not change or replace anything in `pycbc.coordinates.space`;
it is purely additive.
"""
import numpy as np
from astropy.constants import au as ASTRONOMICAL_UNIT

logger = __import__('logging').getLogger(__name__)


EARTH_ORBIT_ANGULAR_FREQUENCY = 1.99098659277e-7  # [rad/s]
# = 2*pi / (sidereal year in seconds, 31558149.7635456); pinned as a
# literal (not derived via `import lal`) per project convention of
# avoiding a lal dependency in new code.

# Named constant instead of an np.deg2rad() call directly in the default
# constructor argument below, to avoid a ruff B008 warning.
_TAIJI_LEAD_ANGLE = np.deg2rad(20.0)


def _real_earth_ecliptic_longitude(t=0.0):
    """Real Earth's ecliptic longitude at SSB time `t` [s], from
    `pycbc.coordinates.space.earth_position_ssb` (astropy-based real
    ephemeris). Imported lazily (not at module level) to avoid a circular
    import: `pycbc.coordinates.space` imports from this module.
    """
    from pycbc.coordinates.space import earth_position_ssb
    return earth_position_ssb(t)[1]


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


__all__ = [
    'LisaEqualArmOrbit',
    'TaijiEqualArmOrbit',
]
