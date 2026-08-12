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
"""

import copy
import unittest

import numpy
from utils import simple_exit

from pycbc.detector import Detector
from pycbc.distributions import (
    CosAngle,
    JointDistribution,
    SinAngle,
    Uniform,
    UniformAngle,
)
from pycbc.inference import models
from pycbc.noise import noise_from_psd
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.waveform import get_td_waveform

TC = 1187008882.42840
FLOW, SEGLEN, SRATE = 25., 32, 2048
INJ = dict(mass1=1.4, mass2=1.35, distance=60., inclination=0.5,
           ra=1.7, dec=-0.4, polarization=0.3, coa_phase=1.1)


class TestMargAccuracy(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        flen = int(SRATE * SEGLEN / 2) + 1
        psd = aLIGOZeroDetHighPower(flen, 1. / SEGLEN, FLOW)
        hp, hc = get_td_waveform(
            approximant='TaylorF2', f_lower=FLOW, delta_t=1. / SRATE,
            **{k: INJ[k] for k in ('mass1', 'mass2', 'distance',
                                   'inclination', 'coa_phase')})
        peak = float(hp.sample_times[numpy.argmax(abs(hp.data) ** 2
                                                  + abs(hc.data) ** 2)])
        cls.data, cls.psds = {}, {}
        seed = 3
        for ifo in ['H1', 'L1', 'V1']:
            ts = noise_from_psd(int(SEGLEN * SRATE), 1. / SRATE, psd,
                                seed=seed)
            seed += 101
            ts._epoch = TC - SEGLEN / 2
            signal = Detector(ifo).project_wave(
                hp, hc, INJ['ra'], INJ['dec'], INJ['polarization'])
            signal.start_time += TC - peak
            cls.data[ifo] = ts.add_into(signal).to_frequencyseries()
            cls.psds[ifo] = psd
        cls.flow = {ifo: FLOW for ifo in cls.data}
        cls.static = dict(mass1=INJ['mass1'], mass2=INJ['mass2'],
                          f_lower=FLOW, approximant='TaylorF2',
                          ra=INJ['ra'], dec=INJ['dec'])
        cls.point = {'distance': INJ['distance'],
                     'inclination': INJ['inclination'], 'tc': TC,
                     'polarization': INJ['polarization'],
                     'ra': INJ['ra'], 'dec': INJ['dec']}

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
        """Quadrupling the points must halve the error.

        This is the property that lets the number of points be chosen to
        meet an accuracy, rather than by habit.
        """
        seen = {}
        for npoint in (64, 256, 1024):
            _, sd = self.spread(marginalize_vector_params='polarization',
                                marginalize_vector_samples=npoint)
            seen[npoint] = sd

        self.assertLess(seen[1024], seen[64],
                        "error did not fall with points: %s" % seen)
        for coarse, fine in ((64, 256), (256, 1024)):
            gain = seen[coarse] / seen[fine]
            self.assertGreater(gain, 1.4, "expected about 2: %s" % seen)
            self.assertLess(gain, 3.0, "expected about 2: %s" % seen)

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
            marginalize_vector_samples=512,
            marginalize_sky_initial_samples=1e5)
        _, alone = self.spread(
            sky=True, nseed=6,
            marginalize_vector_params='ra,dec,polarization',
            marginalize_vector_samples=512)

        self.assertLess(together, 2.0,
                        "sky marginalization should be usable when drawn "
                        "with tc, got a spread of %.3f" % together)
        self.assertGreater(alone, 10 * together,
                           "drawing the sky from the prior alone should be "
                           "much worse: %.3f vs %.3f" % (alone, together))

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

        steps = numpy.abs(numpy.diff(values))
        for coarse, fine in zip(steps, steps[1:], strict=False):
            self.assertLess(fine, coarse,
                            "refining the time grid should settle: %s"
                            % values)
        self.assertLess(steps[-1], 0.1,
                        "not converged by 16384 Hz: %s" % values)

    def test_default_number_of_points_is_accurate_enough(self):
        """The shipped default must not dominate the likelihood error.

        An error well below 0.1 in the log likelihood is invisible next to
        the width of a posterior, so the default should sit under that.
        """
        _, sd = self.spread(marginalize_vector_params='polarization')
        self.assertLess(sd, 0.1, "default settings give a spread of %.3f"
                        % sd)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMargAccuracy))

if __name__ == '__main__':
    from astropy.utils import iers
    iers.conf.auto_download = False
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
