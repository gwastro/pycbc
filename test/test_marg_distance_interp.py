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
"""Tests of the distance marginalization interpolant.

The grid is evenly spaced in the log of each coordinate, so the cell a
point falls in is calculated instead of searched for. That has to leave
the answer as close to the integral it stands in for as before; outside
the interpolated range the value is dropped to zero, which biases the
result, so the run has to say so; and what comes back for a single point
has to be recognizable to the caller as a single point, or the answer is
marginalized a second time over nothing.
"""

import logging
import unittest

import numpy
from scipy.special import i0
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
VSAMPLES = 1000
# modest enough that the reference below can be summed without logsumexp
# and the Bessel function taken without its exponential scaling
PLAIN = [(30., 100.), (40., 200.), (25., 60.), (60., 400.), (35., 90.)]


def integrate_over_distance(sh, hh, phase):
    """The marginalization written out, in plain arithmetic.

    Nothing here comes from the code under test: the amplitude goes as the
    inverse distance, so at a rescaling r the inner products are r*sh and
    r**2*hh, and marginalizing is the weighted mean of the likelihood ratio
    over the grid. Marginalizing the overall phase as well replaces
    exp(r*sh) with the modified Bessel function of the first kind.
    """
    r = RESCALE
    if phase:
        lr = i0(r * abs(sh)) * numpy.exp(-0.5 * r * r * hh)
    else:
        lr = numpy.exp(r * sh.real - 0.5 * r * r * hh)
    return float(numpy.log((DIST_MARG[1] * lr).sum()))


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
        # without phase marginalization the inner product keeps its sign
        # rather than its magnitude, so the queries are not the same set
        cls.nophase = staticmethod(setup_distance_marg_interpolant(
            DIST_MARG, phase=False, snr_range=SNR_RANGE, density=DENSITY))
        cls.nophase_truth = numpy.array(
            [marginalize_likelihood(a, b, distance=DIST_MARG, phase=False)
             for a, b in zip(cls.sh, cls.hh)])

    def test_the_marginalization_is_the_integral_it_claims(self):
        """Anchor both routes to an independently written integral.

        The interpolant is checked against marginalize_likelihood
        everywhere else, which is the same module: if the brute force path
        were wrong they would agree and both be wrong. This is the only
        assertion here that does not use it as its reference.
        """
        for sh, hh in PLAIN:
            for phase in (True, False):
                want = integrate_over_distance(numpy.complex128(sh), hh,
                                               phase)
                brute = marginalize_likelihood(numpy.complex128(sh), hh,
                                               distance=DIST_MARG,
                                               phase=phase)
                self.assertAlmostEqual(
                    brute, want, places=10,
                    msg="brute force distance marginalization is %s, the "
                        "integral is %s (sh=%s hh=%s phase=%s)"
                        % (brute, want, sh, hh, phase))
                interp = self.interp if phase else self.nophase
                self.assertLess(
                    abs(interp(sh, hh) - want), 0.01,
                    "the interpolant is %s, the integral is %s (sh=%s "
                    "hh=%s phase=%s)" % (interp(sh, hh), want, sh, hh, phase))

    def test_it_still_stands_in_for_the_integral(self):
        """The interpolant has to answer what it interpolates.

        A tenth of a nat is far inside what a 100x100 grid resolves, and
        far outside what a change of evaluation route could move it by, so
        this fails on a wrong index and passes on a coarse grid.
        """
        for name, interp, truth in (
                ('phase', self.interp, self.truth),
                ('no phase', self.nophase, self.nophase_truth)):
            worst = numpy.abs(interp(self.sh, self.hh) - truth).max()
            self.assertLess(worst, 0.1, "%s: off the marginalized likelihood "
                            "by %s nats" % (name, worst))

    def test_a_single_point_is_not_marginalized_a_second_time(self):
        """A single point has nothing to marginalize over.

        marginalize_likelihood decides whether its input is a vector of
        drawn points by asking whether it is a float, and the interpolant
        hands back a length-one array. Handed that, it folds the value
        through the vector-marginalization weight and the likelihood comes
        out low by exactly log(marginalize_vector_samples).
        """
        one = self.interp(float(self.sh[0]), float(self.hh[0]))
        self.assertIsInstance(one, float)
        self.assertEqual(
            marginalize_likelihood(float(self.sh[0]), float(self.hh[0]),
                                   logw=-numpy.log(VSAMPLES),
                                   interpolator=self.interp, distance=True,
                                   phase=True),
            one)
        # an array query comes back as one value per point, and a single
        # point has to be the value the array route gives for it -- the
        # scalar case is a different return, not a different calculation
        many = self.interp(self.sh[:3], self.hh[:3])
        self.assertEqual(numpy.shape(many), (3,))
        self.assertEqual(one, many[0])

    def test_a_zero_dimensional_query_is_also_a_single_point(self):
        """numpy.asarray(5.0) is one point, but it is not a float."""
        one = self.interp(numpy.asarray(self.sh[0]),
                          numpy.asarray(self.hh[0]))
        self.assertIsInstance(one, float)
        self.assertEqual(one, self.interp(float(self.sh[0]),
                                          float(self.hh[0])))

    def test_out_of_range_is_minus_inf_and_said_once(self):
        """The warning must not change the value: still -inf out of range.

        On its own interpolant, since the warning is said once per
        interpolant and any other test that leaves the grid spends it.
        """
        interp = setup_distance_marg_interpolant(
            DIST_MARG, phase=True, snr_range=SNR_RANGE, density=(20, 20))
        with self.assertNoLogs(level='WARNING'):
            interp(float(self.sh[0]), float(self.hh[0]))

        with self.assertLogs(level='WARNING') as caught:
            self.assertEqual(interp(1e9, 1e12), -numpy.inf)
            v = interp(numpy.array([self.sh[0], 1e9]),
                       numpy.array([self.hh[0], 1e12]))
        # exactly one warning despite two excursions
        self.assertEqual(len(caught.records), 1)
        self.assertIn('snr_range', caught.records[0].getMessage())
        self.assertTrue(numpy.isfinite(v[0]))
        self.assertEqual(v[1], -numpy.inf)

    def test_a_negative_inner_product_is_not_a_nan(self):
        """sh keeps its sign when the phase is not marginalized over.

        The index is the log of the inner product, so a negative one has
        to be held inside the grid before the log rather than after, or it
        is a nan that the bounds mask never gets to overwrite.
        """
        logging.disable(logging.CRITICAL)
        try:
            # raise rather than warn, so taking the log of a negative fails
            # here even though the bounds mask would overwrite the nan
            with numpy.errstate(invalid='raise'):
                v = self.nophase(numpy.full(5, -50.0), numpy.full(5, 100.0))
        finally:
            logging.disable(logging.NOTSET)
        self.assertFalse(numpy.isnan(v).any())
        self.assertTrue(numpy.all(v == -numpy.inf))


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestMargDistanceInterp))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
