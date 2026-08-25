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
"""Tests the two routes through the distance marginalization interpolant.

The grid is evenly spaced in the log of each coordinate, so the cell a
point falls in can be calculated instead of searched for, which is what
the vector route does. It is a different spline from the one FITPACK
fits, so what has to hold is that it interpolates the same grid to the
same accuracy, and that a query long enough to take it gets the answer a
short one would.
"""

import unittest
import warnings

import numpy
from scipy import ndimage

from utils import simple_exit

from pycbc.inference.models.tools import (setup_distance_marg_interpolant,
                                          marginalize_likelihood)

# a stand-in for dist_ref/dist_locs and its weights; the interpolant does
# not care where the distance grid came from
SAMPLES = 500
RESCALE = numpy.linspace(1.0, 0.2, SAMPLES)
DIST_MARG = (RESCALE, numpy.ones(SAMPLES) / SAMPLES)
# narrow enough that a small grid resolves it, so the test measures the
# evaluation rather than the density
SNR_RANGE = (5, 15)
DENSITY = (100, 100)
# the wrapper sends anything shorter than this to the spline
MANY = 100


class TestMargDistanceInterp(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # staticmethod, or looking it up on an instance would bind it and
        # pass the test case in as the first inner product
        cls.interp = staticmethod(setup_distance_marg_interpolant(
            DIST_MARG, phase=True, snr_range=SNR_RANGE, density=DENSITY))
        smin, smax = RESCALE.min(), RESCALE.max()
        lo = SNR_RANGE[0] ** 2.0 / smax
        hi = SNR_RANGE[1] ** 2.0 / smin
        numpy.random.seed(11)
        # inside the grid, which is where the wrapper returns a number
        cls.sh = numpy.exp(numpy.random.uniform(
            numpy.log(lo * 1.2), numpy.log(hi * 0.8), 200))
        cls.hh = numpy.exp(numpy.random.uniform(
            numpy.log(lo * lo / smax * 1.2), numpy.log(hi / smin * 0.8), 200))
        cls.truth = numpy.array(
            [marginalize_likelihood(a, b, distance=DIST_MARG, phase=True)
             for a, b in zip(cls.sh, cls.hh)])

    def scalar_route(self):
        """One point at a time, which is under the length that switches"""
        return numpy.array([self.interp(float(a), float(b))
                            for a, b in zip(self.sh, self.hh)])

    def test_a_long_query_takes_the_calculated_index(self):
        """The route being claimed has to be the one a long query runs.

        Without this the tests below would pass on a wrapper that quietly
        sent everything to the spline.
        """
        calls = []
        real = ndimage.map_coordinates

        def counted(*args, **kwargs):
            calls.append(len(args[1][0]))
            return real(*args, **kwargs)

        ndimage.map_coordinates = counted
        try:
            self.interp(self.sh[:MANY - 1], self.hh[:MANY - 1])
            self.assertEqual(calls, [], "a short query left the spline")
            self.interp(self.sh, self.hh)
            self.assertEqual(calls, [len(self.sh)],
                             "a long query did not use the calculated index")
        finally:
            ndimage.map_coordinates = real

    def test_the_two_routes_agree(self):
        """Both interpolate the same grid, so they answer the same.

        They are not the same spline and do not agree to the last bit;
        they have to agree to well inside the error either makes against
        the function being interpolated, which is asserted below.
        """
        gap = numpy.abs(self.interp(self.sh, self.hh) - self.scalar_route())
        worst = numpy.abs(self.scalar_route() - self.truth).max()
        self.assertLess(
            gap.max(), worst,
            "the routes differ by %s, more than the %s either is from the "
            "function it interpolates" % (gap.max(), worst))

    def test_neither_route_is_further_from_the_integral(self):
        """The change is meant to cost nothing in accuracy.

        Both are measured against the integral they approximate, computed
        without any interpolation, and the calculated-index route may not
        be the worse of the two.
        """
        vector = numpy.abs(self.interp(self.sh, self.hh) - self.truth)
        scalar = numpy.abs(self.scalar_route() - self.truth)
        self.assertLess(
            vector.max(), 1e-2,
            "interpolation error %s is too large for this grid to say "
            "anything" % vector.max())
        self.assertLessEqual(
            vector.max(), 2.0 * scalar.max(),
            "calculated index gave %s against the spline's %s"
            % (vector.max(), scalar.max()))

    def test_a_query_below_the_grid_is_rejected_quietly(self):
        """sh is not positive when the phase is not marginalized over.

        Such a point is out of range and the caller replaces it with
        -inf, but the index is calculated from a logarithm, which would
        warn and produce a nan on the way there. The route has to hold
        the point inside the grid first.
        """
        sh = numpy.concatenate([self.sh, [-3.0, -1e-8]])
        hh = numpy.concatenate([self.hh, [self.hh[0], self.hh[1]]])
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            v = self.interp(sh, hh)
        self.assertTrue(numpy.isneginf(v[-2:]).all(),
                        "a negative inner product gave %s" % v[-2:])
        self.assertFalse(numpy.isnan(v).any(), "the route produced a nan")


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestMargDistanceInterp))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
