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
These are the unittests for the pycbc.waveform module
"""
import sys
import pycbc
import unittest
import numpy
from pycbc.types import *
from pycbc.scheme import *
from pycbc.events import *
import pycbc.fft
from utils import parse_args_all_schemes, simple_exit

_scheme, _context = parse_args_all_schemes("Threshold")

from pycbc.events.threshold_cpu import threshold_numpy
trusted_threshold = threshold_numpy

class TestThreshold(unittest.TestCase):
    def setUp(self,*args):
        self.context = _context
        self.scheme = _scheme
        r = numpy.random.uniform(low=-1, high=1.0, size=2**20)
        i = numpy.random.uniform(low=-1, high=1.0, size=2**20)
        v = r + i*1.0j
        self.series = Array(v, dtype=complex64)
        self.threshold = 1.3
        self.locs, self.vals = trusted_threshold(self.series, self.threshold)
        self.tolerance = 1e-6
        print(len(self.locs), len(self.vals))

    def test_threshold(self):
        with self.context:
            locs, vals = threshold(self.series, self.threshold)
            self.assertTrue((locs == self.locs).all())
            self.assertTrue((vals == self.vals).all())
            print(len(locs), len(vals))
suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestThreshold))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
