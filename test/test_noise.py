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
import numpy
import pycbc.psd
from utils import simple_exit
from pycbc.noise.reproduceable import noise_from_string
from pycbc.fft.fftw import set_measure_level
from hashlib import md5
set_measure_level(0)

class TestNoise(unittest.TestCase):
    def setUp(self,*args):
        self.ts = noise_from_string('aLIGOZeroDetHighPower',
                                    0, 100,
                                    sample_rate=1024,
                                    seed=0,
                                    low_frequency_cutoff=1.0,
                                    filter_duration=64)

    def test_consistent_result(self):
        # This just checks that the result hasn't changed. If it has
        # you should find out why
        summ = self.ts.sum()
        comp = 4.597515648402546e-19
        diff = abs(summ - comp)
        self.assertTrue(diff < 1e-30)

    def test_noise_psd(self):
        p = self.ts.psd(4)
        p2 = pycbc.psd.from_string('aLIGOZeroDetHighPower', len(p),
                                  p.delta_f, 1.0)
        kmin = int(1.0 / p.delta_f)
        kmax = int(500 / p.delta_f)
        ratio = p[kmin:kmax] / p2[kmin:kmax]
        ave = ratio.numpy().mean()
        self.assertAlmostEqual(ave, 1, 2)

suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestNoise))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
