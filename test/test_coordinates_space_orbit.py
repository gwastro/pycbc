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
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestOptionalLisaorbitsDuckTyping))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
