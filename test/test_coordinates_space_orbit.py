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
Regression tests for `pycbc.coordinates.space_orbit.NumericOrbits`.

`NumericOrbits` interpolates a tabulated constellation orbit (times and
spacecraft positions, optionally velocities) with a scipy B-spline, and
reads/writes pycbc's own HDF5 orbit format. The reference orbit used here is
the "Equal arm analytic orbit" of the LISA Data Challenge manual
(LISA-LCST-SGS-MAN-001, Sec. 8.1.1, Eq. 48-52), reimplemented directly in
this file, so these tests compare two independent implementations.
"""
import os
import tempfile
import h5py
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

    def test_rejects_mismatched_shapes(self):
        with self.assertRaises(ValueError):
            space_orbit.NumericOrbits(self.t_grid, self.positions[:-1])
        with self.assertRaises(ValueError):
            space_orbit.NumericOrbits(
                self.t_grid, self.positions,
                velocities=self.positions[:, :, :2])

    def test_velocity_is_the_derivative_of_the_position_spline(self):
        """With no velocities supplied, `compute_velocity` is the analytic
        derivative of the position spline; it must agree with a central
        finite difference of `compute_position` on the same instance.
        """
        dt = 1.0
        fd = (self.numeric.compute_position(self.t_query + dt)
              - self.numeric.compute_position(self.t_query - dt)) / (2 * dt)
        analytic = self.numeric.compute_velocity(self.t_query)
        self.assertLess(numpy.max(numpy.abs(analytic - fd)), 1e-3)

    def test_acceleration_is_the_second_derivative(self):
        """`compute_acceleration` is the second analytic derivative of the
        same spline, and must agree with a central finite difference of
        `compute_velocity`.
        """
        dt = 1.0
        fd = (self.numeric.compute_velocity(self.t_query + dt)
              - self.numeric.compute_velocity(self.t_query - dt)) / (2 * dt)
        analytic = self.numeric.compute_acceleration(self.t_query)
        self.assertLess(numpy.max(numpy.abs(analytic - fd)), 1e-6)

    def test_explicit_velocities_are_used_verbatim(self):
        """When velocities are supplied they must be interpolated in their
        own right, not re-derived from the positions.
        """
        marker = numpy.full_like(self.positions, 1234.5)
        numeric = space_orbit.NumericOrbits(
            self.t_grid, self.positions, velocities=marker)
        got = numeric.compute_velocity(self.t_query)
        self.assertLess(numpy.max(numpy.abs(got - 1234.5)), 1e-6)


class TestNumericOrbitsFileIO(unittest.TestCase):
    """`to_file`, paired with `from_file`, is what lets an orbit built in
    memory be saved once and then referenced by a PE config's `orbit-file`
    option.
    """
    def setUp(self):
        self.exact_orbit = AnalyticEqualArmOrbit()
        self.t_grid = numpy.arange(0.0, 3.15e7, 86400.0)
        self.tmpdir = tempfile.mkdtemp()

    def test_to_file_round_trip_without_velocities(self):
        """With no explicit velocities at construction, the `velocities`
        dataset must be omitted (not written as zeros or anything else), so
        a reloaded instance re-derives velocities from the position spline
        the same way the original did.
        """
        orbit = space_orbit.NumericOrbits(
            self.t_grid, self.exact_orbit.compute_position(self.t_grid))
        path = os.path.join(self.tmpdir, 'roundtrip_no_vel.hdf5')
        orbit.to_file(path)

        with h5py.File(path, 'r') as f:
            self.assertIn('t', f)
            self.assertIn('positions', f)
            self.assertNotIn('velocities', f)

        reloaded = space_orbit.NumericOrbits.from_file(path)
        query_t = numpy.linspace(self.t_grid[5], self.t_grid[-5], 30)
        self.assertLess(numpy.max(numpy.abs(
            reloaded.compute_position(query_t)
            - orbit.compute_position(query_t))), 1e-6)
        self.assertLess(numpy.max(numpy.abs(
            reloaded.compute_velocity(query_t)
            - orbit.compute_velocity(query_t))), 1e-6)
        self.assertLess(numpy.max(numpy.abs(
            reloaded.compute_acceleration(query_t)
            - orbit.compute_acceleration(query_t))), 1e-6)

    def test_to_file_round_trip_with_explicit_velocities(self):
        positions = self.exact_orbit.compute_position(self.t_grid)
        dt = 1.0
        velocities = (self.exact_orbit.compute_position(self.t_grid + dt)
                      - self.exact_orbit.compute_position(self.t_grid - dt)
                      ) / (2 * dt)
        orbit = space_orbit.NumericOrbits(
            self.t_grid, positions, velocities=velocities)
        path = os.path.join(self.tmpdir, 'roundtrip_with_vel.hdf5')
        orbit.to_file(path)

        with h5py.File(path, 'r') as f:
            self.assertIn('velocities', f)

        reloaded = space_orbit.NumericOrbits.from_file(path)
        query_t = numpy.linspace(self.t_grid[5], self.t_grid[-5], 30)
        self.assertLess(numpy.max(numpy.abs(
            reloaded.compute_velocity(query_t)
            - orbit.compute_velocity(query_t))), 1e-6)

    def test_to_file_group_option_round_trip(self):
        orbit = space_orbit.NumericOrbits(
            self.t_grid, self.exact_orbit.compute_position(self.t_grid))
        path = os.path.join(self.tmpdir, 'grouped.hdf5')
        orbit.to_file(path, group='lisa')
        orbit.to_file(path, group='taiji', mode='a')

        for group in ('lisa', 'taiji'):
            reloaded = space_orbit.NumericOrbits.from_file(path, group=group)
            query_t = numpy.linspace(self.t_grid[5], self.t_grid[-5], 30)
            self.assertLess(numpy.max(numpy.abs(
                reloaded.compute_position(query_t)
                - orbit.compute_position(query_t))), 1e-6)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestNumericOrbits))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestNumericOrbitsFileIO))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
