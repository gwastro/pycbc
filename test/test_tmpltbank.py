# Copyright (C) 2013 Ian Harry
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
These are the unittests for the pycbc.tmpltbank module
"""

import os
import numpy
import pycbc.tmpltbank
import pycbc.psd
from pycbc.types import Array
import difflib
import sys

import unittest
from utils import parse_args_cpu_only, simple_exit

# This will return whatever is appropriate, depending on whether this
# particular instance of the unittest was called for CPU, CUDA, or OpenCL
parse_args_cpu_only("Template bank module")

class TmpltbankTestClass(unittest.TestCase):
    def setUp(self):
        # Where are my data files?
        if os.path.isfile('test/data/ZERO_DET_high_P.txt'):
            self.dataDir = 'test/data/'
        elif os.path.isfile('data/ZERO_DET_high_P.txt'):
            self.dataDir = 'data/'
        else:
            self.assertTrue(False, msg="Cannot find data files!")

        self.deltaF = 0.1
        self.f_low = 15
        self.f_upper = 2000
        self.sampleRate = 4096

        self.segLen = 1./self.deltaF
        self.psdSize = int(self.segLen * self.sampleRate) / 2. + 1

        self.psd = pycbc.psd.from_txt('%sZERO_DET_high_P.txt' %(self.dataDir),\
                self.psdSize, self.deltaF, self.f_low, is_asd_file=True)

        self.evals, self.evecs, self.gs, self.moments = \
              pycbc.tmpltbank.determine_eigen_directions(\
                    self.psd, 'taylorF4_45PN', 70, self.f_low, self.f_upper, \
                    self.deltaF)

    def test_eigen_directions(self):
        evalsStock = Array(numpy.loadtxt('%sstockEvals.dat'%(self.dataDir)))
        evecsStock = Array(numpy.loadtxt('%sstockEvecs.dat'%(self.dataDir)))
        evalsCurr = Array(self.evals[self.f_upper])
        evecsCurr = Array(self.evecs[self.f_upper])
        errMsg = "pycbc.tmpltbank.determine_eigen_directions has failed "
        errMsg += "sanity check."
        self.assertTrue(evalsStock.almost_equal_elem(evalsCurr, tol=1E-4),\
                        msg=errMsg)
        self.assertTrue(evecsStock.almost_equal_elem(evecsCurr, tol=1E-4),\
                        msg=errMsg)

    def tearDown(self):
        pass

suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TmpltbankTestClass))
 
if  __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)

