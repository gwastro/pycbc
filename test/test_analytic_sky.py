# Copyright (C) 2026 Alexander Nitz
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
"""Unit tests for the analytic delay-to-sky extrinsic draw."""

import unittest

import numpy
from utils import simple_exit

from pycbc.detector import Detector
from pycbc.distributions import CosAngle, JointDistribution, Uniform
from pycbc.inference.models.tools import DistMarg
from pycbc.types import TimeSeries


class TestDrawInCells(unittest.TestCase):
    """The per-sample window, at the bounds the sky draw actually reaches.

    Each sample gets its own bounds, because the second detector's cell
    has to lie within one light crossing of wherever the first was drawn.
    Under peak locking the locked region is narrower than that crossing,
    so the bounds are usually the whole series; a window holding no cell
    at all needs two regions further apart than the crossing, which takes
    a mis-locked detector.
    """

    def setUp(self):
        numpy.random.seed(4)
        self.n = 64
        # a peak in the middle, so restricting the bounds has to move the draw
        self.logw = -0.02 * (numpy.arange(self.n) - 32.0) ** 2.0
        self.rng = numpy.random.default_rng(9)

    def draw(self, lo, hi, size=2000):
        return DistMarg._draw_in_cells(
            self.logw, numpy.full(size, lo), numpy.full(size, hi), self.rng
        )

    def test_a_window_holding_no_cell_is_rejected(self):
        """Nothing to draw, so the sample must carry no weight."""
        idx, mass = self.draw(20, 20)
        self.assertTrue(numpy.isneginf(mass).all(), "an empty window returned a weight")
        self.assertTrue(
            (idx >= 0).all() and (idx < self.n).all(),
            "a rejected sample still has to return a valid cell",
        )

    def test_the_draw_stays_inside_the_window(self):
        idx, mass = self.draw(10, 20)
        self.assertGreaterEqual(idx.min(), 10)
        self.assertLess(idx.max(), 20)
        self.assertTrue(numpy.isfinite(mass).all())

    def test_the_whole_series_is_the_usual_case(self):
        """Peak locking leaves the bounds at the ends most of the time."""
        idx, mass = self.draw(0, self.n)
        self.assertTrue(numpy.isfinite(mass).all())
        self.assertGreater(len(numpy.unique(idx)), 5)

    def test_the_mass_is_the_probability_within_the_window(self):
        """What it returns has to be the log probability it drew with."""
        lo, hi = 10, 20
        idx, mass = self.draw(lo, hi, size=20000)
        p = numpy.exp(self.logw[lo:hi] - self.logw.max())
        want = numpy.log(p / p.sum())
        self.assertLess(numpy.abs(mass - want[idx - lo]).max(), 1e-12)


class TestEllipseSlice(unittest.TestCase):
    """The segment of dt13 that a fixed dt12 leaves physical.

    The delays live inside the ellipse dt^T S dt <= 1, whose interior is
    the set of delay pairs a real sky direction can produce. Fixing dt12
    cuts a segment out of it, and what comes back has to be that segment:
    its ends on the ellipse, its interior inside, and nothing at all when
    the dt12 asked for is already outside.
    """

    def setUp(self):
        # a positive definite S, as the detector geometry always gives
        self.smat = numpy.array([[9.0e3, 2.0e3], [2.0e3, 7.0e3]])

    def quad(self, dt12, dt13):
        s = self.smat
        return s[0, 0] * dt12**2 + 2 * s[0, 1] * dt12 * dt13 + s[1, 1] * dt13**2

    def test_the_ends_sit_on_the_ellipse(self):
        dt12 = numpy.linspace(-0.008, 0.008, 40)
        mid, half, bad = DistMarg._ellipse_slice(self.smat, dt12)
        good = ~bad
        self.assertGreater(good.sum(), 10, "nothing to test")
        for end in (mid[good] - half[good], mid[good] + half[good]):
            self.assertLess(
                numpy.abs(self.quad(dt12[good], end) - 1.0).max(),
                1e-9,
                "a segment end is not on the ellipse",
            )

    def test_the_middle_is_inside(self):
        dt12 = numpy.linspace(-0.008, 0.008, 40)
        mid, _, bad = DistMarg._ellipse_slice(self.smat, dt12)
        good = ~bad
        self.assertLessEqual(self.quad(dt12[good], mid[good]).max(), 1.0)

    def test_a_dt12_outside_the_ellipse_has_no_slice(self):
        """Past the widest delay the baseline allows, nothing is physical."""
        widest = (1.0 / numpy.linalg.inv(self.smat)[0, 0]) ** -0.5
        dt12 = numpy.array([widest * 1.01, widest * 2.0, -widest * 1.5])
        _, _, bad = DistMarg._ellipse_slice(self.smat, dt12)
        self.assertTrue(bad.all(), "a dt12 outside the ellipse got a slice")
        inside = numpy.array([0.0, widest * 0.5])
        _, _, bad2 = DistMarg._ellipse_slice(self.smat, inside)
        self.assertFalse(bad2.any(), "a dt12 inside the ellipse got none")


class TestAnalyticSkyDraw(unittest.TestCase):
    """The draw itself, tilted by made-up SNR series.

    Whatever tilts the draw, its weights are prior over proposal, so they
    average to one over the draws and reweight them back to the prior.
    """

    TC = 1187008882.4
    NUM = 200000

    def draw(self, ifos, tilt=1.0, sky=None, seed=0):
        numpy.random.seed(seed)
        marg = DistMarg()
        marg.data = dict.fromkeys(ifos)
        marg._current_params = None
        if sky is None:
            sky = Uniform(ra=(0, 2 * numpy.pi)), CosAngle(dec=None)
        prior = JointDistribution(
            ["tc", "ra", "dec"],
            Uniform(tc=(self.TC - 0.05, self.TC + 0.05)),
            *sky,
        )
        marg.setup_marginalization(
            ["tc", "ra", "dec"],
            prior=prior,
            marginalize_vector_params="tc,ra,dec",
            marginalize_vector_samples=self.NUM,
            marginalize_sky_analytic=True,
        )
        # a peak in each detector, at times that need not be consistent,
        # on top of noise
        times = numpy.arange(-0.2, 0.2, 1 / 4096.0)
        snrs = {}
        for k, ifo in enumerate(ifos):
            peak = numpy.exp(-(((times - 0.01 * (k - 1)) / 0.004) ** 2))
            loglr = tilt * (peak + 0.2 * numpy.random.random(len(times)))
            snrs[ifo] = TimeSeries(
                (2 * loglr) ** 0.5 + 0j, delta_t=1 / 4096.0, epoch=self.TC - 0.2
            )
        ok = marg.analytic_sky_draw(snrs, ifos, self.NUM)
        return marg, ok

    def test_the_weights_average_to_one(self):
        for ifos in (
            ["H1"],
            ["H1", "L1"],
            ["H1", "L1", "V1"],
            ["H1", "L1", "V1", "K1"],
        ):
            with self.subTest(ifos=ifos):
                marg, ok = self.draw(ifos)
                self.assertTrue(ok)
                logw = marg.marginalize_vector_params["logw_partial"]
                self.assertFalse(numpy.isnan(logw).any())
                self.assertAlmostEqual(numpy.exp(logw).mean(), 1.0, delta=0.01)

    def test_the_weighted_draws_follow_the_prior(self):
        marg, _ = self.draw(["H1", "L1", "V1"])
        vec = marg.marginalize_vector_params
        weight = numpy.exp(vec["logw_partial"])
        sindec = numpy.sin(vec["dec"])
        for value, want in (
            (sindec, 0.0),
            (sindec**2, 1 / 3.0),
            (numpy.cos(vec["dec"]) * numpy.cos(vec["ra"]), 0.0),
            ((vec["tc"] - self.TC) / 0.05, 0.0),
        ):
            self.assertAlmostEqual((weight * value).mean(), want, delta=0.01)

    def test_the_stored_factors_are_those_of_the_sky_returned(self):
        """The models use the stored factors or recompute them from the sky"""
        ifos = ["H1", "L1", "V1"]
        marg, _ = self.draw(ifos)
        vec = marg.marginalize_vector_params
        fplus, _, delay = marg.precalc_antenna_factors
        for ifo in ifos:
            det = Detector(ifo, reference_time=self.TC)
            for i in numpy.flatnonzero(numpy.isfinite(vec["logw_partial"]))[:50]:
                ra, dec, tc = vec["ra"][i], vec["dec"][i], vec["tc"][i]
                # up to the earth's rotation across the tc prior
                self.assertAlmostEqual(
                    det.time_delay_from_earth_center(ra, dec, tc),
                    delay[ifo][i],
                    delta=1e-7,
                )
                fp, _ = det.antenna_pattern(ra, dec, 0, tc)
                self.assertAlmostEqual(fp, fplus[ifo][i], delta=1e-4)

    def test_a_sky_prior_that_is_not_isotropic_is_followed(self):
        """The draw is isotropic, so the weights carry the prior's shape."""
        sky = Uniform(ra=(1.0, 2.0)), Uniform(dec=(-0.5, 0.8))
        for ifos in (["H1", "L1"], ["H1", "L1", "V1"]):
            with self.subTest(ifos=ifos):
                marg, ok = self.draw(ifos, sky=sky)
                self.assertTrue(ok)
                vec = marg.marginalize_vector_params
                weight = numpy.exp(vec["logw_partial"])
                self.assertAlmostEqual(weight.mean(), 1.0, delta=0.03)
                self.assertAlmostEqual((weight * vec["ra"]).mean(), 1.5, delta=0.04)
                self.assertAlmostEqual((weight * vec["dec"]).mean(), 0.15, delta=0.03)

    def test_colocated_detectors_leave_it_to_the_map(self):
        _, ok = self.draw(["H1", "H2"])
        self.assertFalse(ok)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestDrawInCells))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestEllipseSlice))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestAnalyticSkyDraw))

if __name__ == "__main__":
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
