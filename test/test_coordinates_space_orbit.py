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
Regression tests for `pycbc.coordinates.space_orbit`.

The core claim being tested is backward compatibility: the generic,
orbit-provider-based machinery in `space_orbit` must reproduce the existing,
hard-coded circular-orbit LISA functions in `pycbc.coordinates.space` exactly
(to floating point precision) when fed the same analytic circular orbit that
those functions assume. The reference orbit used here is the "Equal arm
analytic orbit" defined in the LISA Data Challenge manual
(LISA-LCST-SGS-MAN-001, Sec. 8.1.1, Eq. 48-52); this is the same closed-form
per-spacecraft trajectory implemented by `lisaorbits.EqualArmlengthOrbits`,
reimplemented here directly (not imported from `lisaorbits`) so this test has
no optional dependency.

Beyond the LISA special case, this file also checks that `constellation_frame`
generalizes correctly to constellations with genuinely different geometry:
Taiji (heliocentric, like LISA but leading the Earth by 20 degrees instead of
trailing it, with a different arm length/eccentricity) and TianQin (a fast,
rigidly-rotating triangle around a geocentric guiding center, in a plane
fixed in inertial space rather than precessing over a year). The Taiji
fixture uses the same first-order-in-eccentricity Keplerian expansion as the
LDC manual's Eq. 48-52 above (Rubbo, Cornish & Poujade 2004, Phys. Rev. D 69,
082003), just with Taiji's own arm length/eccentricity and a +20 degree lead
angle instead of LISA's trailing offset. The TianQin fixture follows Hu et al
2018 (Class. Quantum Grav. 35, 095008), simplified to a pure circular
heliocentric guiding center coincident with the Earth (consistent with how
the LISA/Taiji guiding center above is also treated as the zero-eccentricity
limit of the same family of orbits).
"""
import numpy
import unittest
from astropy.constants import au

from pycbc.coordinates import space
from pycbc.coordinates import space_orbit
from utils import simple_exit


seed = 8202
numpy.random.seed(seed)

OMEGA_0 = 1.99098659277e-7  # LISA/Taiji-like orbital angular frequency [rad/s]
ARMLENGTH = 2.5e9  # matches pycbc.detector.space._space_detectors['LISA']
SEMI_MAJOR_AXIS = au.value
ECCENTRICITY = ARMLENGTH / (2 * SEMI_MAJOR_AXIS * numpy.sqrt(3))
T0 = space.TIME_OFFSET_20_DEGREES


class AnalyticEqualArmOrbit:
    """Reference implementation of the LDC manual's 'Equal arm analytic
    orbit' (LISA-LCST-SGS-MAN-001, Sec. 8.1.1, Eq. 48-52). Used only to
    validate `space_orbit.constellation_frame` and related functions against
    the existing analytic functions in `pycbc.coordinates.space`, which
    assume this same orbit but only track the (eccentricity-independent)
    guiding center.
    """
    def compute_position(self, t, sc=(1, 2, 3)):
        t = numpy.atleast_1d(numpy.asarray(t, dtype=float))
        sc = numpy.atleast_1d(sc)
        alpha = OMEGA_0 * (t + T0)
        out = numpy.empty((len(t), len(sc), 3))
        for k, n in enumerate(sc):
            beta_n = (n - 1) * 2 * numpy.pi / 3.0
            out[:, k, 0] = SEMI_MAJOR_AXIS * numpy.cos(alpha) \
                + SEMI_MAJOR_AXIS * ECCENTRICITY * (
                    numpy.sin(alpha) * numpy.cos(alpha) * numpy.sin(beta_n)
                    - (1 + numpy.sin(alpha) ** 2) * numpy.cos(beta_n))
            out[:, k, 1] = SEMI_MAJOR_AXIS * numpy.sin(alpha) \
                + SEMI_MAJOR_AXIS * ECCENTRICITY * (
                    numpy.sin(alpha) * numpy.cos(alpha) * numpy.cos(beta_n)
                    - (1 + numpy.cos(alpha) ** 2) * numpy.sin(beta_n))
            out[:, k, 2] = -numpy.sqrt(3) * SEMI_MAJOR_AXIS * ECCENTRICITY \
                * numpy.cos(alpha - beta_n)
        return out


# Taiji shares the LISA fixture's functional form (the same
# first-order-in-eccentricity Keplerian expansion as the LDC manual's
# Eq. 48-52 above, per Rubbo, Cornish & Poujade 2004, Phys. Rev. D 69,
# 082003), differing only in arm length/eccentricity and in leading (rather
# than trailing) the Earth-like guiding center by 20 degrees.
TAIJI_ARMLENGTH = 3.0e9  # matches pycbc.detector.space._space_detectors
TAIJI_ECCENTRICITY = TAIJI_ARMLENGTH / (2 * SEMI_MAJOR_AXIS * numpy.sqrt(3))
TAIJI_LEAD_ANGLE = numpy.deg2rad(20.0)


class AnalyticTaijiOrbit:
    """Reference implementation of the Taiji heliocentric orbit (first order
    in eccentricity, per Rubbo, Cornish & Poujade 2004, Phys. Rev. D 69,
    082003). Used only to check that `constellation_frame`/`NumericOrbits`
    handle a constellation genuinely different from the LISA fixture above
    (different arm length, eccentricity, and heliocentric lead angle).
    """
    def compute_position(self, t, sc=(1, 2, 3)):
        t = numpy.atleast_1d(numpy.asarray(t, dtype=float))
        sc = numpy.atleast_1d(sc)
        alpha = OMEGA_0 * t + TAIJI_LEAD_ANGLE
        out = numpy.empty((len(t), len(sc), 3))
        for k, n in enumerate(sc):
            beta_n = (n - 1) * 2 * numpy.pi / 3.0
            out[:, k, 0] = SEMI_MAJOR_AXIS * numpy.cos(alpha) \
                + SEMI_MAJOR_AXIS * TAIJI_ECCENTRICITY * (
                    numpy.sin(alpha) * numpy.cos(alpha) * numpy.sin(beta_n)
                    - (1 + numpy.sin(alpha) ** 2) * numpy.cos(beta_n))
            out[:, k, 1] = SEMI_MAJOR_AXIS * numpy.sin(alpha) \
                + SEMI_MAJOR_AXIS * TAIJI_ECCENTRICITY * (
                    numpy.sin(alpha) * numpy.cos(alpha) * numpy.cos(beta_n)
                    - (1 + numpy.cos(alpha) ** 2) * numpy.sin(beta_n))
            out[:, k, 2] = -numpy.sqrt(3) * SEMI_MAJOR_AXIS \
                * TAIJI_ECCENTRICITY * numpy.cos(alpha - beta_n)
        return out


# TianQin: a fast-rotating rigid triangle (per Hu et al 2018, Class.
# Quantum Grav. 35, 095008) whose plane is fixed in inertial space (pointing
# towards the calibration source RX J0806.3+1527), around a guiding center
# that coincides with the Earth. The guiding center itself is simplified
# here to a pure circular heliocentric orbit, matching the zero-eccentricity
# treatment of the LISA/Taiji guiding center used elsewhere in this file
# (the Earth's real orbital eccentricity is irrelevant to what this test
# checks).
TIANQIN_ARMLENGTH = 1.7e8  # matches pycbc.detector.space._space_detectors
TIANQIN_LAMBDA_S = numpy.deg2rad(120.5)  # direction to RX J0806.3+1527
TIANQIN_BETA_S = numpy.deg2rad(-4.7)
TIANQIN_F_SC = 1.0 / (3.65 * 86400.0)  # ~3.65 day rotation period


class AnalyticTianQinOrbit:
    """Reference implementation of the TianQin geocentric orbit (per Hu et
    al 2018, Class. Quantum Grav. 35, 095008). Used to check that
    `constellation_frame` correctly derives a *non-precessing* constellation
    plane, in contrast to LISA/Taiji above.
    """
    def compute_position(self, t, sc=(1, 2, 3), initial_orbit_phase=0.0):
        t = numpy.atleast_1d(numpy.asarray(t, dtype=float))
        sc = numpy.atleast_1d(sc)
        earth = SEMI_MAJOR_AXIS * numpy.stack([
            numpy.cos(OMEGA_0 * t), numpy.sin(OMEGA_0 * t),
            numpy.zeros_like(t)], axis=-1)
        out = numpy.empty((len(t), len(sc), 3))
        for k, n in enumerate(sc):
            kappa_n = 2 * numpy.pi / 3.0 * (n - 1) + initial_orbit_phase
            phase = 2 * numpy.pi * TIANQIN_F_SC * t + kappa_n
            xn = (TIANQIN_ARMLENGTH / numpy.sqrt(3)) * (
                numpy.sin(TIANQIN_BETA_S) * numpy.cos(TIANQIN_LAMBDA_S)
                * numpy.sin(phase)
                + numpy.sin(TIANQIN_LAMBDA_S) * numpy.cos(phase))
            yn = (TIANQIN_ARMLENGTH / numpy.sqrt(3)) * (
                numpy.sin(TIANQIN_BETA_S) * numpy.sin(TIANQIN_LAMBDA_S)
                * numpy.sin(phase)
                - numpy.cos(TIANQIN_LAMBDA_S) * numpy.cos(phase))
            zn = -(TIANQIN_ARMLENGTH / numpy.sqrt(3)) \
                * numpy.cos(TIANQIN_BETA_S) * numpy.sin(phase)
            out[:, k, 0] = earth[:, 0] + xn
            out[:, k, 1] = earth[:, 1] + yn
            out[:, k, 2] = earth[:, 2] + zn
        return out


def _arm_lengths(orbit, t):
    pos = orbit.compute_position(t)
    r1, r2, r3 = pos[:, 0], pos[:, 1], pos[:, 2]
    return (numpy.linalg.norm(r1 - r2, axis=-1),
            numpy.linalg.norm(r2 - r3, axis=-1),
            numpy.linalg.norm(r1 - r3, axis=-1))


class TestConstellationFrame(unittest.TestCase):
    def setUp(self):
        self.orbit = AnalyticEqualArmOrbit()
        self.precision = 1e-6  # relative to ~1 AU length scale / O(1) rotation
        self.times = numpy.random.uniform(0.0, 3.15e7, size=50)

    def test_centroid_matches_lisa_position_ssb(self):
        """The constellation centroid derived from the 3 spacecraft
        positions must match the existing (eccentricity-independent)
        guiding-center formula `lisa_position_ssb`, since the average of the
        Equal Arm orbit's eccentricity terms over the 3 spacecraft is
        exactly zero at every order.
        """
        centroid, _ = space_orbit.constellation_frame(self.times, self.orbit)
        expected = numpy.array([
            space.lisa_position_ssb(t, T0)[0].flatten().astype(float)
            for t in self.times
        ])
        maxdiff = numpy.max(numpy.abs(centroid - expected))
        self.assertLess(maxdiff, self.precision * SEMI_MAJOR_AXIS,
                        f'constellation centroid does not match '
                        f'lisa_position_ssb; max abs diff {maxdiff}')

    def test_rotation_matches_rotation_matrix_ssb_to_lisa(self):
        """The rotation matrix derived from the 3 spacecraft positions must
        match the existing analytic `rotation_matrix_ssb_to_lisa`.
        """
        _, rotation = space_orbit.constellation_frame(self.times, self.orbit)
        for i, t in enumerate(self.times):
            alpha = OMEGA_0 * (t + T0)
            expected = space.rotation_matrix_ssb_to_lisa(alpha)
            maxdiff = numpy.max(numpy.abs(rotation[i] - expected))
            self.assertLess(maxdiff, self.precision,
                            f'rotation matrix at t={t} does not match '
                            f'rotation_matrix_ssb_to_lisa; max abs diff '
                            f'{maxdiff}')

    def test_t_detector_from_ssb_matches_t_lisa_from_ssb(self):
        lam, beta = 1.234, 0.456
        k_ssb = space.localization_to_propagation_vector(
            lam, beta, use_astropy=False).flatten()
        for t_ssb in self.times[:10]:
            expected = space.t_lisa_from_ssb(t_ssb, lam, beta, T0)
            derived = space_orbit.t_detector_from_ssb(
                t_ssb, k_ssb, self.orbit)
            self.assertAlmostEqual(derived, expected, places=3,
                msg=f't_detector_from_ssb does not match t_lisa_from_ssb '
                    f'at t_ssb={t_ssb}')

    def test_t_ssb_from_t_detector_matches_t_ssb_from_t_lisa(self):
        lam, beta = 1.234, 0.456
        k_ssb = space.localization_to_propagation_vector(
            lam, beta, use_astropy=False).flatten()
        for t_lisa in self.times[:10]:
            expected = space.t_ssb_from_t_lisa(t_lisa, lam, beta, T0)
            derived = space_orbit.t_ssb_from_t_detector(
                t_lisa, k_ssb, self.orbit)
            self.assertAlmostEqual(derived, expected, places=3,
                msg=f't_ssb_from_t_detector does not match '
                    f't_ssb_from_t_lisa at t_lisa={t_lisa}')


class TestNumericOrbits(unittest.TestCase):
    def setUp(self):
        self.orbit = AnalyticEqualArmOrbit()
        # coarse grid, similar density to a typical lisaorbits orbit file
        self.t_grid = numpy.linspace(0.0, 3.15e7, 400)
        self.positions = self.orbit.compute_position(self.t_grid)
        self.numeric = space_orbit.NumericOrbits(self.t_grid, self.positions)
        self.t_query = numpy.random.uniform(1e5, 3.1e7, size=20)

    def test_position_interpolation_accuracy(self):
        """Interpolated positions must reproduce the analytic orbit to a
        small fraction of the ~1 AU length scale, at off-grid query times.
        """
        interpolated = self.numeric.compute_position(self.t_query)
        exact = self.orbit.compute_position(self.t_query)
        maxdiff = numpy.max(numpy.abs(interpolated - exact))
        self.assertLess(maxdiff, 1e-6 * SEMI_MAJOR_AXIS,
                        f'NumericOrbits position interpolation error too '
                        f'large: {maxdiff} m')

    def test_constellation_frame_consistent_with_numeric_orbits(self):
        """Feeding `constellation_frame` a `NumericOrbits` instance built
        from a dense sampling of the analytic orbit must reproduce the
        result of feeding it the analytic orbit provider directly.
        """
        centroid_num, rotation_num = space_orbit.constellation_frame(
            self.t_query, self.numeric)
        centroid_ana, rotation_ana = space_orbit.constellation_frame(
            self.t_query, self.orbit)
        self.assertLess(
            numpy.max(numpy.abs(centroid_num - centroid_ana)),
            1e-6 * SEMI_MAJOR_AXIS)
        self.assertLess(
            numpy.max(numpy.abs(rotation_num - rotation_ana)), 1e-6)

    def test_rejects_mismatched_shapes(self):
        with self.assertRaises(ValueError):
            space_orbit.NumericOrbits(self.t_grid, self.positions[:-1])
        with self.assertRaises(ValueError):
            space_orbit.NumericOrbits(
                self.t_grid, self.positions,
                velocities=self.positions[:, :, :2])


class TestTaijiOrbit(unittest.TestCase):
    def setUp(self):
        self.orbit = AnalyticTaijiOrbit()
        self.times = numpy.random.uniform(0.0, 3.15e7, size=50)

    def test_arm_length_matches_design_value(self):
        for d in _arm_lengths(self.orbit, self.times):
            self.assertLess(
                numpy.max(numpy.abs(d - TAIJI_ARMLENGTH)), 1.0,
                'Taiji arm length deviates from the 3e6 km design value '
                'by more than 1 m')

    def test_leads_earth_like_guiding_center_by_20_degrees(self):
        """Taiji's guiding center must lead a coincident, zero-eccentricity
        circular orbit (the same idealized Earth proxy used to define the
        Taiji fixture's own zeroth-order term) by exactly 20 degrees, at
        every time (alpha'' = alpha - beta + 20 deg).
        """
        centroid, _ = space_orbit.constellation_frame(self.times, self.orbit)
        earth_like = SEMI_MAJOR_AXIS * numpy.stack([
            numpy.cos(OMEGA_0 * self.times), numpy.sin(OMEGA_0 * self.times),
            numpy.zeros_like(self.times)], axis=-1)
        cos_angle = numpy.sum(centroid * earth_like, axis=-1) / (
            numpy.linalg.norm(centroid, axis=-1)
            * numpy.linalg.norm(earth_like, axis=-1))
        angle_deg = numpy.rad2deg(numpy.arccos(numpy.clip(cos_angle, -1, 1)))
        self.assertLess(numpy.max(numpy.abs(angle_deg - 20.0)), 1e-6)

    def test_rotation_orthonormal(self):
        _, rotation = space_orbit.constellation_frame(self.times, self.orbit)
        for r in rotation:
            self.assertLess(
                numpy.max(numpy.abs(r @ r.T - numpy.eye(3))), 1e-8)

    def test_round_trip_time_delay(self):
        lam, beta = 0.7, -0.3
        k_ssb = space.localization_to_propagation_vector(
            lam, beta, use_astropy=False).flatten()
        for t_ssb in self.times[:10]:
            t_det = space_orbit.t_detector_from_ssb(t_ssb, k_ssb, self.orbit)
            t_ssb_roundtrip = space_orbit.t_ssb_from_t_detector(
                t_det, k_ssb, self.orbit)
            self.assertAlmostEqual(t_ssb_roundtrip, t_ssb, places=3)


class TestTianQinOrbit(unittest.TestCase):
    def setUp(self):
        self.orbit = AnalyticTianQinOrbit()
        self.times = numpy.random.uniform(0.0, 3.15e7, size=50)

    def test_arm_length_matches_design_value(self):
        for d in _arm_lengths(self.orbit, self.times):
            self.assertLess(
                numpy.max(numpy.abs(d - TIANQIN_ARMLENGTH)), 1.0,
                'TianQin arm length deviates from the 1.7e5 km design '
                'value by more than 1 m')

    def test_centroid_matches_earthlike_guiding_center(self):
        """TianQin's 3 spacecraft are equally spaced in phase around their
        guiding center, so their average must cancel the fast-rotation term
        exactly, leaving just the (here, circular) Earth-like guiding
        center.
        """
        centroid, _ = space_orbit.constellation_frame(self.times, self.orbit)
        earth_like = SEMI_MAJOR_AXIS * numpy.stack([
            numpy.cos(OMEGA_0 * self.times), numpy.sin(OMEGA_0 * self.times),
            numpy.zeros_like(self.times)], axis=-1)
        self.assertLess(numpy.max(numpy.abs(centroid - earth_like)), 1e-3)

    def test_constellation_plane_normal_is_time_independent(self):
        """Unlike LISA/Taiji (whose constellation plane precesses over a
        year), TianQin's plane is fixed in inertial space (pointing towards
        RX J0806.3+1527): the normal derived by `constellation_frame` must
        be the same at every sampled time, even though the in-plane axes
        rotate rapidly (~3.65 day period).
        """
        year_times = numpy.linspace(0.0, 3.15e7, 50)
        _, rotation = space_orbit.constellation_frame(year_times, self.orbit)
        normal = rotation[:, 2, :]
        x_axis = rotation[:, 0, :]
        self.assertLess(
            numpy.max(numpy.abs(normal - normal[0])), 1e-6,
            'TianQin constellation plane normal is not time-independent')
        # the in-plane axis, in contrast, must rotate substantially over a
        # year (~100 fast rotations), i.e. this is not a degenerate check
        self.assertGreater(numpy.max(numpy.abs(x_axis - x_axis[0])), 0.5)

    def test_rotation_orthonormal(self):
        _, rotation = space_orbit.constellation_frame(self.times, self.orbit)
        for r in rotation:
            self.assertLess(
                numpy.max(numpy.abs(r @ r.T - numpy.eye(3))), 1e-8)

    def test_round_trip_time_delay(self):
        lam, beta = 2.1, 0.5
        k_ssb = space.localization_to_propagation_vector(
            lam, beta, use_astropy=False).flatten()
        for t_ssb in self.times[:10]:
            t_det = space_orbit.t_detector_from_ssb(t_ssb, k_ssb, self.orbit)
            t_ssb_roundtrip = space_orbit.t_ssb_from_t_detector(
                t_det, k_ssb, self.orbit)
            self.assertAlmostEqual(t_ssb_roundtrip, t_ssb, places=3)


class TestNumericOrbitsGeneralizesBeyondLisa(unittest.TestCase):
    """`NumericOrbits` must interpolate Taiji- and TianQin-shaped orbits
    accurately too, not just the LISA special case in `TestNumericOrbits`
    above. TianQin's much shorter (~3.65 day) rotation period requires a
    denser sampling grid than LISA/Taiji's yearly one to interpolate well.
    """
    def test_taiji_interpolation_accuracy(self):
        orbit = AnalyticTaijiOrbit()
        t_grid = numpy.linspace(0.0, 3.15e7, 400)
        numeric = space_orbit.NumericOrbits(
            t_grid, orbit.compute_position(t_grid))
        t_query = numpy.random.uniform(1e5, 3.1e7, size=20)
        interpolated = numeric.compute_position(t_query)
        exact = orbit.compute_position(t_query)
        self.assertLess(numpy.max(numpy.abs(interpolated - exact)),
                        1e-6 * SEMI_MAJOR_AXIS)

    def test_tianqin_interpolation_accuracy(self):
        orbit = AnalyticTianQinOrbit()
        # dense enough to resolve the ~3.65 day rotation period
        t_grid = numpy.linspace(0.0, 3.15e7, 20000)
        numeric = space_orbit.NumericOrbits(
            t_grid, orbit.compute_position(t_grid))
        t_query = numpy.random.uniform(1e5, 3.1e7, size=20)
        interpolated = numeric.compute_position(t_query)
        exact = orbit.compute_position(t_query)
        self.assertLess(numpy.max(numpy.abs(interpolated - exact)),
                        1e-6 * TIANQIN_ARMLENGTH)


class TestOptionalLisaorbitsDuckTyping(unittest.TestCase):
    """If `lisaorbits` happens to be installed in the test environment,
    verify that a real `lisaorbits.Orbits` instance can be passed directly
    to `constellation_frame` with no adapter code (duck typing), and that
    the result is physically sane (finite, non-degenerate). This test is
    skipped entirely if `lisaorbits` is not installed; `space_orbit` itself
    never imports it.
    """
    def test_lisaorbits_instance_accepted_directly(self):
        try:
            import lisaorbits
        except ImportError:
            self.skipTest('lisaorbits not installed; skipping duck-typing '
                          'cross-check')
        orbit = lisaorbits.EqualArmlengthOrbits()
        times = numpy.array([1e6, 1e7])
        centroid, rotation = space_orbit.constellation_frame(times, orbit)
        self.assertEqual(centroid.shape, (2, 3))
        self.assertEqual(rotation.shape, (2, 3, 3))
        self.assertTrue(numpy.all(numpy.isfinite(centroid)))
        self.assertTrue(numpy.all(numpy.isfinite(rotation)))
        # rotation matrices must be orthogonal
        for r in rotation:
            self.assertLess(
                numpy.max(numpy.abs(r @ r.T - numpy.eye(3))), 1e-8)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestConstellationFrame))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestNumericOrbits))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestTaijiOrbit))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestTianQinOrbit))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestNumericOrbitsGeneralizesBeyondLisa))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestOptionalLisaorbitsDuckTyping))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
