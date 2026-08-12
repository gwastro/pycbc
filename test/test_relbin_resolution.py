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
"""Tests that the relative binning resolution check tracks the actual error.

These build models directly with parameters chosen to stress the
approximation, rather than running a sampler, so the whole file runs in
well under a minute on cached data.
"""

import unittest

import numpy
from utils import simple_exit

from pycbc.distributions import JointDistribution, SinAngle, Uniform
from pycbc.inference import models
from pycbc.noise import noise_from_psd
from pycbc.psd import aLIGOZeroDetHighPower


class TestRelbinResolution(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # a short simulated signal is enough, and avoids downloading data
        tc = 1187008882.42840
        flow, seglen, srate = 25., 16, 2048
        delta_f = 1. / seglen
        flen = int(srate * seglen / 2) + 1
        psd = aLIGOZeroDetHighPower(flen, delta_f, flow)
        cls.psds, cls.data = {}, {}
        seed = 7
        for ifo in ['H1', 'L1']:
            ts = noise_from_psd(int(seglen * srate), 1. / srate, psd,
                                seed=seed)
            ts._epoch = tc - seglen / 2
            seed += 11
            cls.data[ifo] = ts.to_frequencyseries()
            cls.psds[ifo] = psd
        cls.flow = {'H1': flow, 'L1': flow}
        cls.static = {'mass1': 1.3757, 'mass2': 1.3757, 'f_lower': flow,
                      'approximant': 'TaylorF2', 'polarization': 0.,
                      'ra': 3.44615914, 'dec': -0.40808407, 'tc': tc}
        cls.variable = ['distance', 'inclination']
        cls.prior = JointDistribution(cls.variable, SinAngle(inclination=None),
                                      Uniform(distance=(10, 100)))
        cls.q = {'distance': 42.0, 'inclination': 2.5}

    def model(self, epsilon, static=None):
        return models.Relative(
            list(self.variable), {k: v.copy() for k, v in self.data.items()},
            low_frequency_cutoff=self.flow, psds=self.psds,
            static_params=static or self.static, prior=self.prior,
            fiducial_params={'mass1': 1.3756}, epsilon=epsilon)

    def test_error_falls_with_resolution(self):
        """Adding bins must resolve the ratio better, and at second order.

        The bin width scales with epsilon and the error of a linear
        interpolation goes as the width squared, so halving epsilon should
        cut the reported error by roughly four.
        """
        values = [self.model(eps).check_bin_resolution()
                  for eps in (1.0, 0.5, 0.25)]
        self.assertTrue(values[0] > values[1] > values[2],
                        "error should fall as bins are added: %s" % values)
        for coarse, fine in zip(values, values[1:], strict=False):
            self.assertGreater(coarse / fine, 2.,
                               "expected close to second order: %s" % values)

    def test_error_tracks_the_likelihood_error(self):
        """The reported error must predict the error it stands in for.

        Both are measured against a much finer binning of the same model,
        so this needs no external reference.
        """
        exact = self.model(0.005)
        exact.update(**self.q)
        reference = exact.loglr

        seen = []
        for epsilon in (1.0, 0.25, 0.05):
            model = self.model(epsilon)
            model.update(**self.q)
            seen.append((model.check_bin_resolution(),
                         abs(model.loglr - reference)))

        # the diagnostic is only useful if the two fall together
        for (ecoarse, lcoarse), (efine, lfine) in zip(seen, seen[1:], strict=False):
            self.assertGreater(ecoarse, efine, "diagnostic: %s" % seen)
            self.assertGreater(lcoarse, lfine, "likelihood: %s" % seen)

    def test_flags_a_poorly_placed_fiducial(self):
        """Bins that suit the fiducial need not suit the signal."""
        offset = dict(self.static)
        offset['mass1'], offset['mass2'] = 1.39, 1.36
        good = self.model(0.5).check_bin_resolution()
        bad = self.model(0.5, static=offset).check_bin_resolution()
        self.assertGreater(bad, good * 10)
        self.assertGreater(bad, 1e-3)

    def test_check_leaves_the_model_alone(self):
        """A diagnostic must not move the point the model is sitting at."""
        model = self.model(0.5)
        model.update(**self.q)
        before = model.loglr
        model.check_bin_resolution()
        self.assertEqual(model.current_params['distance'], self.q['distance'])
        self.assertEqual(before, model.loglr)

    def test_check_is_cheap(self):
        """It must cost a handful of waveforms, not a full analysis."""
        model = self.model(0.5)
        value = model.check_bin_resolution(ndraw=3)
        self.assertTrue(numpy.isfinite(value))
        self.assertGreaterEqual(value, 0.)


suite = unittest.TestSuite()
suite.addTest(
    unittest.TestLoader().loadTestsFromTestCase(TestRelbinResolution))

if __name__ == '__main__':
    from astropy.utils import iers
    iers.conf.auto_download = False
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
