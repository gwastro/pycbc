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
from pycbc.pnutils import *
from pycbc.scheme import *
from utils import parse_args_cpu_only

# We only need CPU tests
parse_args_cpu_only("PN Utilities")

class TestUtils(unittest.TestCase):
    def test_mass1_mass2_to_tau0_tau3(self):
        result = mass1_mass2_to_tau0_tau3(3.0,5.0,15.0)
        answer = (8.9886054117032569e-08, 0.00068132339501069089)
        self.assertAlmostEqual(result[0]/answer[0],1,places=6)
        self.assertAlmostEqual(result[1]/answer[1],1,places=6)

    def test_tau0_tau3_to_mtotal_eta(self):
        result = tau0_tau3_to_mtotal_eta(10,5,15.0)
        answer = [0.00052771444173289135, 0.019562392458119287]
        self.assertAlmostEqual(result[0]/answer[0],1,places=6)
        self.assertAlmostEqual(result[1]/answer[1],1,places=6)

    def test_tau0_tau3_to_mass1_mass2(self):
        result = tau0_tau3_to_mass1_mass2(10,5,15.0)
        answer =  [0.00051718082839390816, 1.0533613338983225e-05]
        self.assertAlmostEqual(result[0]/answer[0],1,places=6)
        self.assertAlmostEqual(result[1]/answer[1],1,places=6)

    def test_mass1_mass2_to_mtotal_eta(self):
        result = mass1_mass2_to_mtotal_eta(5,10)
        answer = [15.0, 0.22222222222222221]
        self.assertAlmostEqual(result[0]/answer[0],1,places=6)
        self.assertAlmostEqual(result[1]/answer[1],1,places=6)

    def test_mass1_mass2_to_mchirp_eta(self):
        result = mass1_mass2_to_mchirp_eta(5,10)
        answer = [6.0836434189320574, 0.22222222222222224]
        self.assertAlmostEqual(result[0]/answer[0],1,places=6)
        self.assertAlmostEqual(result[1]/answer[1],1,places=6)

suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestUtils))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
