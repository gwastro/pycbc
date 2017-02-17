# Copyright (C) 2012  Alex Nitz, Josh Willis
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
These are the unittests for the pycbc.filter.matchedfilter module
"""
import sys
import unittest
from pycbc.types import *
from pycbc.filter import *
from pycbc.scheme import *
from utils import parse_args_all_schemes, simple_exit
from numpy.random import uniform
import scipy.signal
from pycbc.filter.resample import lfilter

_scheme, _context = parse_args_all_schemes("Resampling")

class TestUtils(unittest.TestCase):
    def setUp(self,*args):
        self.scheme = _scheme
        self.context = _context
        self.delta_t = 1.0 / 4096
        self.target_delta_t = 1.0 / 1024
        self.a = TimeSeries([1,2,3,4], delta_t=self.delta_t, dtype=float32)
        self.b = TimeSeries([1,2,3,4], delta_t=self.delta_t, dtype=float64)
        self.c = TimeSeries([1,2,3,4], delta_t=self.delta_t, dtype=complex64)
        self.d = Array([1,2,3,4], dtype=float32)

    if _scheme == 'cpu':
        def test_resample_float32(self):
            ra = resample_to_delta_t(self.a, self.target_delta_t)
            self.assertAlmostEqual(ra[0], 0.00696246)
            ra = resample_to_delta_t(self.a, self.delta_t)
            self.assertAlmostEqual(ra[0], 1)

        def test_resample_float64(self):
            rb = resample_to_delta_t(self.b, self.target_delta_t)
            self.assertAlmostEqual(rb[0], 0.00696246)
            rb = resample_to_delta_t(self.b, self.delta_t)
            self.assertAlmostEqual(rb[0], 1)

    def test_resample_errors(self):
        self.assertRaises(TypeError, resample_to_delta_t, self.c, self.target_delta_t)
        self.assertRaises(TypeError, resample_to_delta_t, self.d, self.target_delta_t)

        if self.scheme != 'cpu':
            with self.context:
                self.assertRaises(TypeError, resample_to_delta_t, self.a, self.target_delta_t)

    def test_lfilter(self):
        "Check our hand written lfilter"
        c = uniform(-10, 10, size=1024)
        ts = uniform(-1, 1, size=4323)
        
        ref = scipy.signal.lfilter(c, 1.0, ts)
        test = lfilter(c, ts)

        # These only agree where there is no fft wraparound
        # so excluded corrupted region from test
        ref = ref[len(c):]
        test = test[len(c):]

        maxreldiff =  ((ref - test) / ref).max()
        
        self.assertTrue(maxreldiff < 1e-7)

suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestUtils))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
