# Copyright (C) 2026  Alex Nitz
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
"""Tests of the wrapper around the distance marginalization interpolant.

The distance-marginalized likelihood is interpolated over a fixed range of
signal to noise ratio; outside it the value is dropped to zero, which
biases the result, so the run has to say so. And what the wrapper returns
for a single point has to be recognizable to its caller as a single point,
or the answer is marginalized a second time over nothing.
"""

import logging
import unittest

import numpy
from utils import simple_exit

from pycbc.inference.models.tools import (marginalize_likelihood,
                                          setup_distance_marg_interpolant)

VSAMPLES = 1000


class TestDistanceInterpolantWrapper(unittest.TestCase):

    def interpolant(self, snr_range=(5, 40)):
        dist_rescale = numpy.linspace(0.5, 2.0, 60)
        weights = numpy.ones(60) / 60
        return setup_distance_marg_interpolant(
            (dist_rescale, weights), snr_range=snr_range, density=(60, 60))

    def test_warns_once_and_only_outside_the_range(self):
        interp = self.interpolant()
        with self.assertNoLogs(level='WARNING'):
            # a signal to noise ratio of about 20, well inside (5, 40)
            interp(200.0, 100.0)
            interp(numpy.array([150.0, 300.0]), numpy.array([80.0, 120.0]))

        with self.assertLogs(level='WARNING') as caught:
            interp(1e9, 1e12)
            interp(1e9, 1e12)
            interp(numpy.array([1e9]), numpy.array([1e12]))
        # exactly one warning despite three excursions
        self.assertEqual(len(caught.records), 1)
        self.assertIn('snr_range', caught.records[0].getMessage().lower()
                      .replace('signal to noise ratio', 'snr_range'))

    def test_out_of_range_still_returns_minus_inf(self):
        """The warning must not change the value: still -inf out of range."""
        interp = self.interpolant()
        logging.disable(logging.CRITICAL)
        try:
            self.assertEqual(interp(1e9, 1e12), -numpy.inf)
            v = interp(numpy.array([200.0, 1e9]), numpy.array([100.0, 1e12]))
        finally:
            logging.disable(logging.NOTSET)
        self.assertTrue(numpy.isfinite(v[0]))
        self.assertEqual(v[1], -numpy.inf)

    def test_a_scalar_query_is_not_marginalized_a_second_time(self):
        """A single point has nothing to marginalize over.

        marginalize_likelihood decides whether its input is a vector of
        drawn points by asking whether it is a float, and a
        zero-dimensional array is not one. Handed that, it folds the value
        through the vector-marginalization weight and the likelihood comes
        out low by exactly log(marginalize_vector_samples) -- 6.9 nats at
        the default. So the wrapper has to return a float for a scalar
        query, and an array only for an array query.
        """
        interp = self.interpolant()
        one = interp(200.0, 100.0)
        self.assertIsInstance(one, float)
        self.assertEqual(
            marginalize_likelihood(200.0, 100.0, logw=-numpy.log(VSAMPLES),
                                   interpolator=interp, distance=True,
                                   phase=True),
            one)

        many = interp(numpy.array([200.0, 150.0]),
                      numpy.array([100.0, 80.0]))
        self.assertEqual(numpy.ndim(many), 1)
        self.assertEqual(len(many), 2)

    def test_a_zero_dimensional_query_is_also_a_single_point(self):
        """numpy.asarray(5.0) is one point, but it is not a float.

        The wrapper picks between the two interpolators on the length of
        the query, and a zero-dimensional array has no length, so deciding
        on isinstance(x, float) rather than its dimension raises TypeError
        before either of them is reached.
        """
        interp = self.interpolant()
        one = interp(numpy.asarray(200.0), numpy.asarray(100.0))
        self.assertIsInstance(one, float)
        self.assertEqual(one, interp(200.0, 100.0))


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestDistanceInterpolantWrapper))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
