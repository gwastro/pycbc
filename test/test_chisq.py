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
from utils import parse_args_all_schemes, simple_exit

_scheme, _context = parse_args_all_schemes("correlate")

from pycbc.vetoes.chisq_cpu import chisq_accum_bin_numpy
from pycbc.vetoes import chisq_accum_bin
trusted_accum = chisq_accum_bin_numpy

class TestChisq(unittest.TestCase):
    def setUp(self,*args):
        self.context = _context
        self.scheme = _scheme
        self.tolerance = 1e-6
        xr = numpy.random.uniform(low=-1, high=1.0, size=2**20)
        xi = numpy.random.uniform(low=-1, high=1.0, size=2**20)
        self.x = Array(xr + xi * 1.0j, dtype=complex64)
        self.z = zeros(2**20, dtype=float32)
        for i in range(0, 4):
            trusted_accum(self.z, self.x)
        
    def test_accum(self):
        with self.context:
            z = zeros(2**20, dtype=float32)
            for i in range(0, 4):
                chisq_accum_bin(z, self.x)
            self.assertTrue(self.z.almost_equal_elem(z, self.tolerance))
            
suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestChisq))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
