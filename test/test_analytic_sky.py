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
"""Tests for the analytic delay-to-sky extrinsic draw.

The helpers it is built from, the draw itself on made-up SNR series, and
the draw through a model, against the pregenerated sky map.
"""

import copy
import unittest

import numpy
from utils import simple_exit

from pycbc.detector import Detector
from pycbc.distributions import CosAngle, JointDistribution, SinAngle, Uniform
from pycbc.inference import models
from pycbc.inference.models.tools import DistMarg
from pycbc.noise import noise_from_psd
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.types import TimeSeries
from pycbc.waveform import get_td_waveform


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
        """A drawn cell outside the window would be an unphysical delay."""
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
        """The ends are the most extreme dt13 a real direction can give."""
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
        """Between the ends, every dt13 is a real direction."""
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

    def draw(self, ifos, tilt=1.0, sky=None, seed=0, tc_half=0.05, frame=None):
        numpy.random.seed(seed)
        marg = DistMarg()
        marg.data = dict.fromkeys(ifos)
        marg._current_params = None
        if frame is not None:
            marg.static_params = {"tc_ref_frame": frame}
        if sky is None:
            sky = Uniform(ra=(0, 2 * numpy.pi)), CosAngle(dec=None)
        prior = JointDistribution(
            ["tc", "ra", "dec"],
            Uniform(tc=(self.TC - tc_half, self.TC + tc_half)),
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
        """The weights are prior over proposal, so their mean is one.

        A wrong normalization or Jacobian anywhere in the draw shows here.
        """
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
        """Moments of the isotropic sky and of the uniform tc prior.

        Weight that is right on average but in the wrong places, such as
        on the wrong mirror root or with the earth turned the wrong way,
        shows here.
        """
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

    def test_tc_in_a_detector_frame(self):
        """With tc in Hanford's frame, the prior and the roots are in it too.

        A tc prior narrower than the light travel time makes the choice
        between the two mirror roots depend on the frame.
        """
        for half in (0.05, 0.004):
            with self.subTest(tc_half=half):
                marg, ok = self.draw(["H1", "L1", "V1"], tc_half=half, frame="H1")
                self.assertTrue(ok)
                vec = marg.marginalize_vector_params
                weight = numpy.exp(vec["logw_partial"])
                # few draws survive the narrow prior, so the spread is 0.01
                self.assertAlmostEqual(weight.mean(), 1.0, delta=0.04)
                self.assertAlmostEqual(
                    (weight * (vec["tc"] - self.TC) / half).mean(), 0.0, delta=0.02
                )

    def test_colocated_detectors_leave_it_to_the_map(self):
        """With no delay to invert, the draw hands over to the map."""
        _, ok = self.draw(["H1", "H2"])
        self.assertFalse(ok)


class TestThroughAModel(unittest.TestCase):
    """The draw as a model uses it, against the pregenerated sky map.

    A BNS in simulated noise in H1, L1 and V1, marginalized over tc, ra
    and dec by RelativeTimeDom. Three detectors, since with two the second
    delay is a free azimuth and with one there is no delay to invert.
    Both draws estimate the same integral, so the means of five calls
    must agree; each call scatters by 0.2 to 0.3 at this number of draws.
    """

    TC = 1187008882.42
    FLOW, SEGLEN, SRATE = 25.0, 32, 2048
    IFOS = ["H1", "L1", "V1"]
    VARIABLE = ["distance", "inclination", "tc", "ra", "dec"]
    POINT = {"distance": 40.0, "inclination": 0.5}
    SKY = {"ra": 1.7, "dec": -0.4, "polarization": 0.3}

    @classmethod
    def make_data(cls, distance):
        hp, hc = get_td_waveform(
            approximant="TaylorF2",
            f_lower=cls.FLOW,
            delta_t=1.0 / cls.SRATE,
            mass1=1.4,
            mass2=1.35,
            distance=distance,
            inclination=0.5,
            coa_phase=1.1,
        )
        hp.start_time += cls.TC
        hc.start_time += cls.TC
        data = {}
        for n, ifo in enumerate(cls.IFOS):
            noise = noise_from_psd(
                int(cls.SEGLEN * cls.SRATE), 1.0 / cls.SRATE, cls.psds[ifo], seed=11 + n
            )
            noise._epoch = cls.TC - cls.SEGLEN / 2
            sig = Detector(ifo).project_wave(
                hp, hc, cls.SKY["ra"], cls.SKY["dec"], cls.SKY["polarization"]
            )
            data[ifo] = noise.add_into(sig).to_frequencyseries()
        return data

    @classmethod
    def setUpClass(cls):
        psd = aLIGOZeroDetHighPower(
            int(cls.SRATE * cls.SEGLEN / 2) + 1, 1.0 / cls.SEGLEN, cls.FLOW
        )
        cls.psds = dict.fromkeys(cls.IFOS, psd)
        cls.data = cls.make_data(40.0)
        cls.static = {
            "mass1": 1.4,
            "mass2": 1.35,
            "f_lower": cls.FLOW,
            "approximant": "TaylorF2",
            "polarization": cls.SKY["polarization"],
        }
        cls.fiducial = {"mass1": 1.4, "tc": cls.TC, **cls.SKY}

    def loglr(self, analytic, data=None, point=None, static=None):
        """The mean of five marginalized likelihoods, and the model."""
        data = self.data if data is None else data
        numpy.random.seed(5)
        prior = JointDistribution(
            list(self.VARIABLE),
            SinAngle(inclination=None),
            Uniform(distance=(10, 200)),
            Uniform(tc=(self.TC - 0.1, self.TC + 0.1)),
            Uniform(ra=(0, 2 * numpy.pi)),
            CosAngle(dec=None),
        )
        model = models.RelativeTimeDom(
            list(self.VARIABLE),
            {i: copy.deepcopy(data[i]) for i in self.IFOS},
            low_frequency_cutoff=dict.fromkeys(self.IFOS, self.FLOW),
            psds=self.psds,
            static_params=self.static if static is None else static,
            prior=prior,
            fiducial_params=self.fiducial,
            epsilon=0.1,
            marginalize_phase=True,
            marginalize_vector_params="tc,ra,dec",
            marginalize_vector_samples=2000,
            sample_rate=4096,
            marginalize_sky_analytic=analytic,
        )
        values = []
        for _ in range(5):
            model.update(**(self.POINT if point is None else point))
            values.append(model.loglr)
        return numpy.array(values), model

    def test_it_agrees_with_the_map(self):
        """The analytic draw is used, and gives the map's answer."""
        analytic, model = self.loglr(True)
        # only the map builds its table of sky points
        self.assertFalse(hasattr(model, "tinfo"))
        sky_map, model = self.loglr(False)
        self.assertTrue(hasattr(model, "tinfo"))
        self.assertAlmostEqual(analytic.mean(), sky_map.mean(), delta=0.6)

    def test_tc_in_a_detector_frame(self):
        """The draw returns tc in the frame it is sampled in.

        Returned at the geocentre, it would be read as Hanford's time and
        give 630.0 here, against 636.5 from the map.
        """
        static = dict(self.static, tc_ref_frame="H1")
        analytic, _ = self.loglr(True, static=static)
        sky_map, _ = self.loglr(False, static=static)
        self.assertAlmostEqual(analytic.mean(), sky_map.mean(), delta=0.6)

    def test_the_map_survives_a_detector_without_signal(self):
        """At twice the distance Virgo has no peak above its noise.

        Its drawn times spread over the whole window, so few delays match
        a sky point. The map used to fall back to the same prior draws on
        every call (108.2 here); it now draws from the louder detectors.
        """
        data = self.make_data(80.0)
        point = dict(self.POINT, distance=80.0)
        analytic, _ = self.loglr(True, data=data, point=point)
        sky_map, _ = self.loglr(False, data=data, point=point)
        self.assertGreater(sky_map.std(), 0.0)
        self.assertAlmostEqual(analytic.mean(), sky_map.mean(), delta=0.6)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestDrawInCells))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestEllipseSlice))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestAnalyticSkyDraw))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestThroughAModel))

if __name__ == "__main__":
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
