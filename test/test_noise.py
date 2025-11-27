# Copyright (C) 2019  Alex Nitz
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

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""
These are the unittests for noise generation
"""
import unittest
import pycbc.psd
from utils import simple_exit
from pycbc.noise.reproduceable import noise_from_string as rep_noise_from_string
from pycbc.fft.fftw import set_measure_level
from pycbc.noise.reproduceable import normal
from pycbc.noise.gaussian import noise_from_string as gaussian_noise_from_string

set_measure_level(0)

class TestNoiseReproduceable(unittest.TestCase):
    def setUp(self,*args):
        self.ts = rep_noise_from_string(
            'aLIGOZeroDetHighPower',
            100,
            200,
            sample_rate=1024,
            seed=0,
            low_frequency_cutoff=1.0,
            filter_duration=64
        )

    def test_consistent_result(self):
        # This just checks that the result hasn't changed. If it has
        # you should find out why
        summ = self.ts.sum()
        comp = 4.265258573533567e-18

        diff = abs(summ - comp)
        self.assertLess(diff, 1e-30)

        std = self.ts.data.std()
        comp = 1.1461606806180472e-20
        diff = abs(std - comp)
        self.assertLess(diff, 1e-30)


    def test_noise_psd(self):
        p = self.ts.psd(4)
        p2 = pycbc.psd.from_string('aLIGOZeroDetHighPower', len(p),
                                  p.delta_f, 1.0)
        kmin = int(1.0 / p.delta_f)
        kmax = int(500 / p.delta_f)
        ratio = p[kmin:kmax] / p2[kmin:kmax]
        ave = ratio.numpy().mean()
        self.assertAlmostEqual(ave, 1, 1)

    def test_noise_reproducible(self):
        ts1 = normal(20, 30, sample_rate=16384, seed=87693)
        ts2 = normal(25, 35, sample_rate=16384, seed=87693)
        self.assertEqual(ts1.time_slice(25, 30), ts2.time_slice(25, 30))

class TestNoiseGaussian(unittest.TestCase):
    def setUp(self,*args):
        self.ts = gaussian_noise_from_string(
            'aLIGOZeroDetHighPower',
            102400,
            delta_t=1./1024,
            seed=1865,
            low_frequency_cutoff=1.0,
        )
        self.ts._epoch = 100

    def test_noise_psd(self):
        p = self.ts.psd(4)
        p2 = pycbc.psd.from_string('aLIGOZeroDetHighPower', len(p),
                                  p.delta_f, 1.0)
        kmin = int(1.0 / p.delta_f)
        kmax = int(500 / p.delta_f)
        ratio = p[kmin:kmax] / p2[kmin:kmax]
        ave = ratio.numpy().mean()
        self.assertAlmostEqual(ave, 1, 1)

    def test_consistent_result(self):
        # This just checks that the result hasn't changed. If it has
        # you should find out why
        summ = self.ts.sum()
        comp = -3.4666738998970245e-33

        diff = abs(summ - comp)
        self.assertLess(diff, 1e-30)

        std = self.ts.data.std()
        comp = 1.3670623471855284e-20

        diff = abs(std - comp)
        self.assertLess(diff, 1e-30)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestNoiseReproduceable))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestNoiseGaussian))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
