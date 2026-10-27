# Copyright (C) 2026 Shichao Wu, Alex Nitz, Alex Correia
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

"""Tests for the spacecraft orbit interface."""

import unittest

import numpy

from pycbc.detector import BaseOrbit
from pycbc.detector.orbits import BaseOrbit as ModuleBaseOrbit


class PolynomialOrbit(BaseOrbit):
    """Polynomial orbit with x=j+t**2, y=j*t and z=-j."""

    def _evaluate(self, t, sc, derivative):
        times = self._prepare_times(t)[:, None]
        labels = (self._sc_indices(sc) + 1)[None, :]
        times, labels = numpy.broadcast_arrays(times, labels)
        zeros = numpy.zeros_like(times)
        if derivative == 0:
            components = labels + times**2, labels * times, -labels
        elif derivative == 1:
            components = 2 * times, labels, zeros
        else:
            components = zeros + 2, zeros, zeros
        return numpy.stack(components, axis=-1)

    def compute_position(self, t, sc=None):
        return self._evaluate(t, sc, 0)

    def compute_velocity(self, t, sc=None):
        return self._evaluate(t, sc, 1)

    def compute_acceleration(self, t, sc=None):
        return self._evaluate(t, sc, 2)


class TestBaseOrbit(unittest.TestCase):
    def setUp(self):
        self.orbit = PolynomialOrbit()

    def test_public_export(self):
        self.assertIs(BaseOrbit, ModuleBaseOrbit)

    def test_base_is_abstract(self):
        with self.assertRaises(TypeError):
            BaseOrbit()

    def test_each_method_is_required(self):
        methods = ('compute_position', 'compute_velocity',
                   'compute_acceleration')
        for missing in methods:
            with self.subTest(missing=missing):
                namespace = {name: getattr(PolynomialOrbit, name)
                             for name in methods if name != missing}
                incomplete = type('IncompleteOrbit', (BaseOrbit,), namespace)
                with self.assertRaises(TypeError):
                    incomplete()

    def test_analytic_metadata(self):
        self.assertIsNone(self.orbit.t_interp)
        self.assertEqual(self.orbit.num_sc, 3)

    def test_interpolation_metadata_is_copied(self):
        times = numpy.array([-5., 0., 2., 10.])
        orbit = PolynomialOrbit(t_interp=times)
        times[0] = 100.
        numpy.testing.assert_array_equal(orbit.t_interp, [-5., 0., 2., 10.])
        self.assertEqual(orbit.t_interp.ndim, 1)

    def test_invalid_interpolation_times(self):
        for times in ([], [1.], [2., 1.], [1., 1.], [[1., 2.]],
                      [0., numpy.nan], [0., numpy.inf], [0., 1j],
                      ['0', '1'], True):
            with self.subTest(times=times), self.assertRaises(ValueError):
                PolynomialOrbit(t_interp=times)

    def test_spacecraft_count(self):
        orbit = PolynomialOrbit(num_sc=numpy.int64(4))
        self.assertEqual(orbit.compute_position(0.).shape, (1, 4, 3))
        numpy.testing.assert_array_equal(orbit._sc_indices(None), [0, 1, 2, 3])

    def test_invalid_spacecraft_count(self):
        for count in (0, -1, 3., True, numpy.bool_(True), '3', None):
            with self.subTest(count=count), self.assertRaises(ValueError):
                PolynomialOrbit(num_sc=count)

    def test_known_position_velocity_acceleration(self):
        expected = ([[[7., 6., -3.]]], [[[4., 3., 0.]]], [[[2., 0., 0.]]])
        for name, values in zip(('compute_position', 'compute_velocity',
                                 'compute_acceleration'), expected):
            with self.subTest(method=name):
                numpy.testing.assert_array_equal(
                    getattr(self.orbit, name)(2., 3), values)

    def test_scalar_and_vector_shapes(self):
        for method in (self.orbit.compute_position,
                       self.orbit.compute_velocity,
                       self.orbit.compute_acceleration):
            for times, labels, shape in ((1., None, (1, 3, 3)),
                                         (1., 2, (1, 1, 3)),
                                         ([1., 2.], [3, 1], (2, 2, 3)),
                                         ([1., 2.], 1, (2, 1, 3))):
                with self.subTest(method=method.__name__, shape=shape):
                    self.assertEqual(method(times, labels).shape, shape)

    def test_order_and_repeats(self):
        times = [2., -1., 2., 0.]
        labels = [3, 1, 3]
        for method in (self.orbit.compute_position,
                       self.orbit.compute_velocity,
                       self.orbit.compute_acceleration):
            result = method(times, labels)
            for i, time in enumerate(times):
                for j, label in enumerate(labels):
                    numpy.testing.assert_array_equal(
                        result[i, j], method(time, label)[0, 0])

    def test_empty_time_axis(self):
        for method in (self.orbit.compute_position,
                       self.orbit.compute_velocity,
                       self.orbit.compute_acceleration):
            self.assertEqual(method([], [1, 2]).shape, (0, 2, 3))

    def test_empty_spacecraft_axis(self):
        labels = numpy.array([], dtype=int)
        for method in (self.orbit.compute_position,
                       self.orbit.compute_velocity,
                       self.orbit.compute_acceleration):
            self.assertEqual(method([0., 1.], labels).shape, (2, 0, 3))

    def test_invalid_evaluation_times(self):
        for times in (numpy.nan, numpy.inf, -numpy.inf, 1j, [[1., 2.]],
                      ['1'], None, True):
            for method in (self.orbit.compute_position,
                           self.orbit.compute_velocity,
                           self.orbit.compute_acceleration):
                with self.subTest(times=times, method=method.__name__):
                    with self.assertRaises(ValueError):
                        method(times)

    def test_invalid_spacecraft_labels(self):
        for labels in (0, -1, 4, 1., 1.5, True, [1, 4], [[1, 2]],
                       numpy.nan, '1', 1j):
            for method in (self.orbit.compute_position,
                           self.orbit.compute_velocity,
                           self.orbit.compute_acceleration):
                with self.subTest(labels=labels, method=method.__name__):
                    with self.assertRaises(ValueError):
                        method(0., labels)

    def test_unsigned_spacecraft_labels(self):
        numpy.testing.assert_array_equal(
            self.orbit._sc_indices(numpy.array([3, 1], dtype=numpy.uint64)),
            [2, 0])


if __name__ == '__main__':
    unittest.main()
