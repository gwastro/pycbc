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
import numpy
from pycbc.pnutils import *
from pycbc.scheme import *
from utils import parse_args_cpu_only, simple_exit

# We only need CPU tests
parse_args_cpu_only("PN Utilities")

class TestUtils(unittest.TestCase):
    def test_mass1_mass2_to_tau0_tau3(self):
        result = mass1_mass2_to_tau0_tau3(3.0,5.0,15.0)
        answer = (63.039052988077955, 2.353532999897545)
        self.assertAlmostEqual(result[0]/answer[0],1,places=6)
        self.assertAlmostEqual(result[1]/answer[1],1,places=6)

    def test_tau0_tau3_to_mtotal_eta(self):
        result = tau0_tau3_to_mtotal_eta(93.84928959285253,2.9198487498891126,20.0)
        answer = [5., 4.*1./5./5.]
        self.assertAlmostEqual(result[0]/answer[0],1,places=6)
        self.assertAlmostEqual(result[1]/answer[1],1,places=6)

    def test_tau0_tau3_to_mass1_mass2(self):
        result = tau0_tau3_to_mass1_mass2(12.410035910174642,0.9266455525603574,30.0)
        answer =  [6., 2.]
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

    def test_mass1_mass2_spin1z_spin2z_to_beta_sigma_gamma(self):
        # with no spin
        result = mass1_mass2_spin1z_spin2z_to_beta_sigma_gamma(1.4, 1.4,
                                                               0., 0.)
        for i in xrange(3):
            self.assertAlmostEqual(result[i], 0, places=6)
        # with spin
        result = mass1_mass2_spin1z_spin2z_to_beta_sigma_gamma(10., 1.4,
                                                               0.9, 0.1)
        answer = [7.208723197, 3.251802285, 243.2697314]
        for r, a in zip(result, answer):
            self.assertAlmostEqual(r / a, 1, places=6)
        result = mass1_mass2_spin1z_spin2z_to_beta_sigma_gamma(5., 5.,
                                                               0.5, -0.7)
        answer = [-0.7833333333, 0.07250000000, -24.59479718]
        for r, a in zip(result, answer):
            self.assertAlmostEqual(r / a, 1, places=6)
        # using array arguments
        mass1 = numpy.array([1.4, 10., 5., 5.])
        mass2 = numpy.array([1.4, 1.4, 5., 5.])
        spin1 = numpy.array([0., 0.9, 0.5, -0.7])
        spin2 = numpy.array([0., 0.1, -0.7, 0.5])
        answer = numpy.array([
            [0., 0., 0.],
            [7.208723197, 3.251802285, 243.2697314],
            [-0.7833333333, 0.07250000000, -24.59479718],
            [-0.7833333333, 0.07250000000, -24.59479718]
        ]).T
        result = mass1_mass2_spin1z_spin2z_to_beta_sigma_gamma(mass1, mass2,
                                                               spin1, spin2)
        for error in (result - answer).ravel():
            self.assertAlmostEqual(error, 0, places=6)

suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestUtils))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
