# Copyright (C) 2026  Shichao Wu, Alex Nitz
#
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
"""
Regression tests for the analytic equal-arm constellation orbits in
`pycbc.coordinates.space_orbit`.

`LisaEqualArmOrbit` and `TaijiEqualArmOrbit` are checked against a
hand-written reimplementation of the same closed-form orbit (the LISA Data
Challenge manual's "Equal arm analytic orbit", LISA-LCST-SGS-MAN-001 Sec.
8.1.1 Eq. 48-52; Rubbo, Cornish & Poujade 2004, Phys. Rev. D 69, 082003),
written independently in this file so the tests compare two implementations
rather than comparing the code to itself. `LisaEqualArmOrbit` is
additionally required to reproduce `pycbc.coordinates.space`'s pre-existing
hard-coded LISA orbit with its default `t0`, since that default is meant to
stand in for the existing code path exactly.
"""
import numpy
import unittest
from astropy.constants import au

from pycbc.coordinates import space
from pycbc.coordinates import space_orbit
from utils import simple_exit


seed = 8202
numpy.random.seed(seed)

OMEGA_0 = 1.99098659277e-7  # 2*pi / sidereal year [rad/s]
ARMLENGTH = 2.5e9  # matches pycbc.detector.space._space_detectors['LISA']
SEMI_MAJOR_AXIS = au.value
ECCENTRICITY = ARMLENGTH / (2 * SEMI_MAJOR_AXIS * numpy.sqrt(3))
T0 = space.TIME_OFFSET_20_DEGREES

# Real Earth ecliptic longitude at t=0, the epoch Taiji's guiding centre is
# anchored to by default.
PHI0_REAL = space.earth_position_ssb(0.0)[1]

TAIJI_ARMLENGTH = 3.0e9  # matches pycbc.detector.space._space_detectors
TAIJI_ECCENTRICITY = TAIJI_ARMLENGTH / (2 * SEMI_MAJOR_AXIS * numpy.sqrt(3))
TAIJI_LEAD_ANGLE = numpy.deg2rad(20.0)


class _LisaGuidingCenterFixture:
    """Mixin providing `_phase(t)` for the LISA fixture's guiding-center
    convention (`OMEGA_0 * (t + T0)`) -- mirrors the shape of
    `space_orbit._LisaGuidingCenter`, but independently written (not
    imported from production), so a bug in the production formula
    wouldn't be invisible here. Combine with
    `_EqualArmConstellationFixture`.
    """
    def _phase(self, t):
        return OMEGA_0 * (t + T0)


class _TaijiGuidingCenterFixture:
    """Mixin providing `_phase(t)` for the Taiji fixture's guiding-center
    convention (`OMEGA_0 * t + PHI0_REAL + TAIJI_LEAD_ANGLE`, anchored to
    real Earth's longitude at t=0 so "leads the Earth by 20 degrees" can
    be checked against the real Earth position, not an arbitrarily-phased
    circular proxy) -- mirrors `space_orbit._TaijiGuidingCenter`,
    independently written. Combine with `_EqualArmConstellationFixture`.
    """
    def _phase(self, t):
        return OMEGA_0 * t + PHI0_REAL + TAIJI_LEAD_ANGLE


class _EqualArmConstellationFixture:
    """Mixin implementing `compute_position` for the LDC manual's 'Equal
    arm analytic orbit' (LISA-LCST-SGS-MAN-001, Sec. 8.1.1, Eq. 48-52) --
    mirrors `space_orbit._EqualArmConstellation`, independently written.
    Shared by `AnalyticEqualArmOrbit`/`AnalyticTaijiOrbit` below (same
    functional form, per Rubbo, Cornish & Poujade 2004, Phys. Rev. D 69,
    082003); the concrete class's own `__init__` must set
    `self.eccentricity`.
    """
    def compute_position(self, t, sc=(1, 2, 3)):
        t = numpy.atleast_1d(numpy.asarray(t, dtype=float))
        sc = numpy.atleast_1d(sc)
        alpha = self._phase(t)
        out = numpy.empty((len(t), len(sc), 3))
        for k, n in enumerate(sc):
            beta_n = (n - 1) * 2 * numpy.pi / 3.0
            out[:, k, 0] = SEMI_MAJOR_AXIS * numpy.cos(alpha) \
                + SEMI_MAJOR_AXIS * self.eccentricity * (
                    numpy.sin(alpha) * numpy.cos(alpha) * numpy.sin(beta_n)
                    - (1 + numpy.sin(alpha) ** 2) * numpy.cos(beta_n))
            out[:, k, 1] = SEMI_MAJOR_AXIS * numpy.sin(alpha) \
                + SEMI_MAJOR_AXIS * self.eccentricity * (
                    numpy.sin(alpha) * numpy.cos(alpha) * numpy.cos(beta_n)
                    - (1 + numpy.cos(alpha) ** 2) * numpy.sin(beta_n))
            out[:, k, 2] = -numpy.sqrt(3) * SEMI_MAJOR_AXIS \
                * self.eccentricity * numpy.cos(alpha - beta_n)
        return out


class AnalyticEqualArmOrbit(_LisaGuidingCenterFixture,
                             _EqualArmConstellationFixture):
    """Reference implementation of the LDC manual's 'Equal arm analytic
    orbit', LISA flavor. Used only to validate `space_orbit.
    constellation_frame` and related functions against the existing
    analytic functions in `pycbc.coordinates.space`, which assume this
    same orbit but only track the (eccentricity-independent) guiding
    center.

    Parameters
    ----------
    eccentricity : float, optional
        Constellation eccentricity. Default `ECCENTRICITY` (LISA's own
        arm length).
    """
    def __init__(self, eccentricity=None):
        self.eccentricity = (
            ECCENTRICITY if eccentricity is None else eccentricity)


class AnalyticTaijiOrbit(_TaijiGuidingCenterFixture,
                          _EqualArmConstellationFixture):
    """Reference implementation of the Taiji heliocentric orbit --
    `AnalyticEqualArmOrbit`'s formula with Taiji's own guiding-center
    convention and eccentricity. Used only to check that
    `constellation_frame`/`NumericOrbits` handle a constellation genuinely
    different from the LISA fixture above (different arm length,
    eccentricity, and heliocentric lead angle).
    """
    def __init__(self):
        self.eccentricity = TAIJI_ECCENTRICITY


def _arm_lengths(orbit, t):
    pos = orbit.compute_position(t)
    r1, r2, r3 = pos[:, 0], pos[:, 1], pos[:, 2]
    return (numpy.linalg.norm(r1 - r2, axis=-1),
            numpy.linalg.norm(r2 - r3, axis=-1),
            numpy.linalg.norm(r1 - r3, axis=-1))


class TestEqualArmAnalyticOrbits(unittest.TestCase):
    """`LisaEqualArmOrbit`/`TaijiEqualArmOrbit` are the production (not
    test-only) analytic reference orbits, intended for prototyping ahead of
    an official numeric orbit product. The key check is that each
    separately implemented production class agrees with the hand-written
    fixture above when given the same reference phase, i.e. that promoting
    the fixtures' math to a real, documented, parameterised class did not
    introduce a transcription error.
    """
    def test_lisa_matches_hardcoded_space_functions(self):
        """With its default `t0`, `LisaEqualArmOrbit` must reproduce the
        constellation centroid of `pycbc.coordinates.space`'s pre-existing
        hard-coded LISA orbit, since that default is meant to be a drop-in
        stand-in for the existing `orbit=None` code path.
        """
        orbit = space_orbit.LisaEqualArmOrbit()
        times = numpy.random.uniform(0.0, 3.15e7, size=20)
        centroid = orbit.compute_position(times).mean(axis=1)
        expected = numpy.array([
            space.lisa_position_ssb(t, orbit.t0)[0].flatten().astype(float)
            for t in times
        ])
        self.assertLess(numpy.max(numpy.abs(centroid - expected)),
                        1e-6 * SEMI_MAJOR_AXIS)

    def test_lisa_matches_independent_fixture_implementation(self):
        production = space_orbit.LisaEqualArmOrbit()
        fixture = AnalyticEqualArmOrbit()
        times = numpy.random.uniform(0.0, 3.15e7, size=20)
        self.assertLess(
            numpy.max(numpy.abs(production.compute_position(times)
                                - fixture.compute_position(times))), 1e-3,
            'LisaEqualArmOrbit does not match the independent '
            'AnalyticEqualArmOrbit test fixture')

    def test_lisa_default_t0_matches_time_offset_20_degrees(self):
        orbit = space_orbit.LisaEqualArmOrbit()
        self.assertEqual(orbit.t0, space.TIME_OFFSET_20_DEGREES)

    def test_lisa_custom_t0_overrides_default(self):
        orbit = space_orbit.LisaEqualArmOrbit(t0=0.0)
        self.assertEqual(orbit.t0, 0.0)
        default_orbit = space_orbit.LisaEqualArmOrbit()
        self.assertGreater(
            numpy.max(numpy.abs(
                orbit.compute_position([1e7])
                - default_orbit.compute_position([1e7]))),
            1e6)

    def test_taiji_matches_independent_fixture_implementation(self):
        production = space_orbit.TaijiEqualArmOrbit(kappa0=PHI0_REAL)
        fixture = AnalyticTaijiOrbit()
        times = numpy.random.uniform(0.0, 3.15e7, size=20)
        self.assertLess(
            numpy.max(numpy.abs(production.compute_position(times)
                                - fixture.compute_position(times))), 1e-3,
            'TaijiEqualArmOrbit does not match the independent '
            'AnalyticTaijiOrbit test fixture')

    def test_taiji_default_kappa0_anchors_to_real_earth(self):
        self.assertAlmostEqual(
            space_orbit.TaijiEqualArmOrbit().kappa0, PHI0_REAL, places=9)

    def test_taiji_custom_kappa0_overrides_default(self):
        orbit = space_orbit.TaijiEqualArmOrbit(kappa0=0.0)
        self.assertEqual(orbit.kappa0, 0.0)
        default_orbit = space_orbit.TaijiEqualArmOrbit()
        # a different kappa0 must give a different position at a fixed
        # time (i.e. the parameter is not silently ignored)
        self.assertGreater(
            numpy.max(numpy.abs(
                orbit.compute_position([1e7])
                - default_orbit.compute_position([1e7]))),
            1e6)

    def test_arm_lengths_match_design_values(self):
        cases = [
            (space_orbit.LisaEqualArmOrbit(), ARMLENGTH),
            (space_orbit.TaijiEqualArmOrbit(), TAIJI_ARMLENGTH),
        ]
        times = numpy.random.uniform(0.0, 3.15e7, size=20)
        for orbit, expected_armlength in cases:
            for d in _arm_lengths(orbit, times):
                self.assertLess(
                    numpy.max(numpy.abs(d - expected_armlength)), 1.0)

    def test_velocity_matches_finite_difference_of_position(self):
        """`compute_velocity` is an exact analytic derivative, not a
        finite-difference approximation -- but a central finite difference
        of `compute_position`, at a small enough step, must still agree
        with it to high precision.
        """
        dt = 1.0  # s
        for orbit in (space_orbit.LisaEqualArmOrbit(),
                      space_orbit.TaijiEqualArmOrbit()):
            t = numpy.random.uniform(1e5, 3.1e7, size=10)
            finite_diff_vel = (orbit.compute_position(t + dt)
                               - orbit.compute_position(t - dt)) / (2 * dt)
            self.assertLess(
                numpy.max(numpy.abs(
                    orbit.compute_velocity(t) - finite_diff_vel)), 1e-3)

    def test_acceleration_matches_finite_difference_of_velocity(self):
        """Same rationale, one derivative order up."""
        dt = 1.0  # s
        for orbit in (space_orbit.LisaEqualArmOrbit(),
                      space_orbit.TaijiEqualArmOrbit()):
            t = numpy.random.uniform(1e5, 3.1e7, size=10)
            finite_diff_acc = (orbit.compute_velocity(t + dt)
                               - orbit.compute_velocity(t - dt)) / (2 * dt)
            self.assertLess(
                numpy.max(numpy.abs(
                    orbit.compute_acceleration(t) - finite_diff_acc)), 1e-6)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestEqualArmAnalyticOrbits))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
