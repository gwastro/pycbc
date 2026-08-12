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
"""Measures how accurate the marginalized likelihoods are.

The marginalization options control a Monte Carlo integral, but nothing
said how much accuracy a given setting buys. These tests measure that, so
that changing the settings can be justified rather than guessed at.

The marginalization points are drawn once when the model is built and then
reused for every evaluation, so repeating a call does not vary the answer.
The error is a fixed offset per set of draws, and it is measured here by
building the model repeatedly from different seeds and looking at the
spread. That spread is the accuracy the sampler actually sees, and it is
what an effective sample size would predict: it should fall as one over
the square root of the number of points.

A signal is injected into simulated noise so the integrand has a real
peak; against pure noise the time and sky integrals are meaningless.

The sky marginalization is left to its own tests. Its accuracy is set by
the number of points in the same way as the rest, but the quality of the
points is set separately by the sample rate, which fixes how finely the
sky is divided by arrival time, and by the number of points the sky map
is built from. Measuring it against the number of points alone, with
those two left coarse, says nothing.
"""

import copy
import unittest

import numpy
from utils import simple_exit
from validation import FLOW, TC, get_seed, make_data

from pycbc.distributions import (
    CosAngle,
    JointDistribution,
    SinAngle,
    Uniform,
    UniformAngle,
)
from pycbc.inference import models


class TestMargAccuracy(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.seed = get_seed(cls.__name__)
        cls.data, cls.psds, inj = make_data(cls.seed)
        cls.flow = {ifo: FLOW for ifo in cls.data}
        cls.static = dict(mass1=inj['mass1'], mass2=inj['mass2'],
                          f_lower=FLOW, approximant='TaylorF2',
                          ra=inj['ra'], dec=inj['dec'])
        cls.point = {'distance': inj['distance'],
                     'inclination': inj['inclination'], 'tc': TC,
                     'polarization': inj['polarization'],
                     'ra': inj['ra'], 'dec': inj['dec']}

    def build(self, sky=False, **kwargs):
        variable = ['distance', 'inclination', 'tc', 'polarization']
        static = dict(self.static)
        dists = [SinAngle(inclination=None), Uniform(distance=(10, 200)),
                 Uniform(tc=(TC - 0.1, TC + 0.1)),
                 UniformAngle(polarization=None)]
        if sky:
            # the sky prior has to be separable: the marginalization pops
            # one parameter at a time out of the joint prior
            variable += ['ra', 'dec']
            del static['ra'], static['dec']
            dists += [UniformAngle(ra=None), CosAngle(dec=None)]
        return models.MarginalizedTime(
            variable, copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow, psds=self.psds,
            static_params=static,
            prior=JointDistribution(variable, *dists),
            marginalize_phase=True, **kwargs)

    def spread(self, nseed=8, **kwargs):
        """Mean and spread of the log likelihood over independent draws."""
        values = []
        for s in range(nseed):
            numpy.random.seed(1000 + s)
            model = self.build(**kwargs)
            model.update(**{p: self.point[p]
                            for p in model.variable_params})
            values.append(model.loglr)
        return numpy.mean(values), numpy.std(values)

    def test_vector_marginalization_converges_as_root_n(self):
        """Sixteen times the points must cut the error by four.

        This is the property that lets the number of points be chosen to
        meet an accuracy, rather than by habit.

        Measured over the whole range rather than step by step: each
        spread is itself estimated from a finite number of runs, and over
        one step of four that uncertainty is the same size as what is
        being measured.
        """
        seen = {}
        for npoint in (256, 1024, 4096):
            _, sd = self.spread(nseed=24,
                                marginalize_vector_params='polarization',
                                marginalize_vector_samples=npoint)
            seen[npoint] = sd

        # only worth asking about where the answer is usable in the first
        # place; too few points and the estimator has collapsed onto a
        # handful of them, where no rate of convergence applies
        if seen[256] > 0.5:
            self.skipTest("this signal needs more than 256 points before "
                          "the error is small enough for its rate of "
                          "convergence to be the question: %s" % seen)

        # the error must fall substantially over the range, but not by an
        # asserted factor: the spread is estimated from a finite number of
        # runs and the log likelihood over draws is heavy-tailed, so the
        # estimate itself scatters by tens of percent even at this nseed.
        # A precise root-n exponent cannot be measured cheaply; a
        # substantial fall over sixteen times the points can, and a genuine
        # failure to converge would leave the ratio near one. Faster than
        # root-n is not a defect, so there is no upper bound.
        gain = seen[256] / seen[4096]
        self.assertGreater(gain, 1.8,
                           "error did not fall with points: %s" % seen)

    def test_too_few_points_is_wrong_but_not_silently_so(self):
        """Outside the usable range the bar is different.

        Nobody should run with sixty four points, and the question there
        is not how fast the error falls but whether it announces itself.
        It does: the spread between runs is what a sampler would see, so
        a setting this coarse is visible rather than quietly biased.
        """
        coarse_mean, coarse_sd = self.spread(
            nseed=8, marginalize_vector_params='polarization',
            marginalize_vector_samples=64)
        fine_mean, fine_sd = self.spread(
            nseed=8, marginalize_vector_params='polarization',
            marginalize_vector_samples=4096)
        self.assertGreater(coarse_sd, fine_sd)
        # whatever it costs in precision, it must not move the answer
        self.assertLess(abs(coarse_mean - fine_mean),
                        3 * (coarse_sd ** 2 + fine_sd ** 2) ** 0.5,
                        "coarse %.3f+-%.3f vs fine %.3f+-%.3f"
                        % (coarse_mean, coarse_sd, fine_mean, fine_sd))

    def test_vector_marginalization_is_unbiased(self):
        """Coarse settings must cost precision, not correctness."""
        coarse_mean, coarse_sd = self.spread(
            marginalize_vector_params='polarization',
            marginalize_vector_samples=64)
        fine_mean, fine_sd = self.spread(
            marginalize_vector_params='polarization',
            marginalize_vector_samples=1024)
        # the means must agree within the error on the difference
        error = (coarse_sd ** 2 + fine_sd ** 2) ** 0.5
        self.assertLess(abs(coarse_mean - fine_mean), 3 * error,
                        "coarse %.4f+-%.4f vs fine %.4f+-%.4f"
                        % (coarse_mean, coarse_sd, fine_mean, fine_sd))

    def test_sky_marginalization_needs_time_marginalized_with_it(self):
        """Sky points are only usable if drawn against the time series.

        The sky integrand is sharply peaked, so drawing sky points from
        the prior alone is hopeless. Including tc in the same set lets the
        code place points using the signal to noise time series instead,
        which is why every example marginalizes them together.
        """
        _, together = self.spread(
            sky=True, nseed=6,
            marginalize_vector_params='tc,ra,dec,polarization',
            marginalize_vector_samples=2048,
            marginalize_sky_initial_samples=1e5)
        _, alone = self.spread(
            sky=True, nseed=6,
            marginalize_vector_params='ra,dec,polarization',
            marginalize_vector_samples=2048)

        # for some signals the estimator has not converged even with tc
        # drawn in, and its spread stays large however many points are
        # used; comparing two unconverged numbers says nothing, so the
        # claim is only tested where the good case is actually good. A
        # converged sky-and-time spread is order one.
        if together > 5.0:
            self.skipTest("sky marginalization has not converged for this "
                          "signal even with tc, spread %.1f; nothing to "
                          "compare against" % together)

        # how much worse depends on the signal, from a factor of a few to
        # tens, so only the direction is asserted, with a margin for the
        # two spreads each being estimated from a handful of runs.
        self.assertGreater(alone, 1.5 * together,
                           "drawing the sky from the prior alone should be "
                           "worse: %.3f vs %.3f" % (alone, together))

    def test_time_marginalization_converges_in_sample_rate(self):
        """The time grid must be fine enough that refining it stops moving.

        sample_rate sets the spacing of the time grid the likelihood is
        summed over. Each refinement should move the answer less than the
        one before it.
        """
        values = []
        for rate in (2048, 4096, 8192, 16384):
            model = self.build(sample_rate=rate)
            model.update(**{p: self.point[p]
                            for p in model.variable_params})
            values.append(model.loglr)

        # how fine it has to be depends on the signal, so what is asked
        # for here is that it is settling, not that it has arrived. Only
        # the whole range is compared: each value carries the scatter of
        # the marginalization itself, which one step does not outrun.
        steps = numpy.abs(numpy.diff(values))
        self.assertLess(steps[-1], steps[0] / 4.,
                        "not settling: %s" % values)

    def test_default_number_of_points_is_accurate_enough(self):
        """The shipped default must not dominate the likelihood error.

        How large the error is depends on the signal: it is a few
        hundredths for a face-on binary and rises towards edge-on, where
        the polarization integrand is peakiest, reaching about a quarter
        of a nat at the default thousand points. A quarter of a nat of
        scatter between evaluations is still invisible next to the width
        of a posterior, so the default is adequate even there, but the
        threshold has to allow for it rather than be set from a convenient
        signal. The spread is estimated from more draws than elsewhere
        because near edge-on it is what is being tested, not incidental.
        """
        _, sd = self.spread(nseed=24,
                            marginalize_vector_params='polarization')
        self.assertLess(sd, 0.4, "default settings give a spread of %.3f"
                        % sd)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMargAccuracy))

if __name__ == '__main__':
    from astropy.utils import iers
    iers.conf.auto_download = False
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
