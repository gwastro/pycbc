# Copyright (C) 2016 Christopher M. Biwer
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
These are the unittests for the pycbc.inference module
"""
import sys
import pycbc
import unittest
import numpy
from pycbc import inference
from utils import parse_args_cpu_only, simple_exit

# tests only need to happen on the CPU
parse_args_cpu_only("Prior")

class TestPriorEvaluator(unittest.TestCase):
    def setUp(self, *args):

        # create a list to hold distributions
        self._dists = []

        # loop over all distributions
        for name,distclass in inference.distributions.iteritems():

            # Guassian distribution
            if name == inference.distributions.Gaussian.name:
                dist = distributions[name](["mass1", "mass2"], low=[10.0,10.0],
                                          high=[50.0,50.0], mean=[25.0,25.0],
                                          var=[10.0,10.0])

            # subclasses of Uniform
            else:
                self._dists.append( distclass(mass1=(0.0,1.0), mass2=(0.0,1.0)) )

    def test_accum(self):
        pass
 
suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestPriorEvaluator))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
