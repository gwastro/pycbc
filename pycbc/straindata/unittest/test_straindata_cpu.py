#!/usr/bin/python
#
# Copyright (C) 2011 Karsten Wiesner
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
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
# pycbc should follow <http://www.python.org/dev/peps/pep-0008/>
#
# vim: tabstop=4:softtabstop=4:shiftwidth=4:expandtab

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

"""
unittest for the StrainData class
"""

import sys
# preliminary hard coded path to packages 
sys.path.append('/Users/kawies/dev/src/pycbc')
sys.path.append('/Users/kawies/dev/src/pycbc/pycbc/straindata')

from straindata_cpu import StrainData as DUT_StrainData
import unittest

#import random

class TestStrainDataCPU(unittest.TestCase):

    # Test Fixture ---------------------------

    def setUp(self):
        print "setUp called"
        
        self.length =   128
        self.segments = 2
        
        self.dut= DUT_StrainData(self.segments, self.length, "H1")

    def tearDown(self):
        print "tearDown called"

    # Sub tests -------------------------------

    def test_1_initial_timeseries(self):
        """
        Test proper initialization, datatype and access of the initial 
        double precision time series for strain data
        """

        self.assertTrue(repr(self.dut.strain_time_series).find("datavectorcpu.real_vector_double_t") >= 0,
        " Wrong type of datavector for straindata at the initial timepoint")
        
        print self.dut.strain_time_series

        print repr(self.dut.strain_time_series)
         



    def test_2_access_datavectors(self):

        self.dut.strain_time_series[0]= 1.234
        self.assertEqual(self.dut.strain_time_series[0], 1.234)
        
        #self.dut.strain_time_series[self.length-1]= 1.234
        #self.assertEqual(self.dut.strain_time_series[self.length-1], 1.234)

        #self.dut.strain_freq_series[0]= 1.23456789
        #self.assertEqual(self.dut.strain_freq_series[0], 1.23456789)
        
        #self.dut.strain_freq_series[self.length-1]= 1.23456789
        #self.assertEqual(self.dut.strain_freq_series[self.length-1], 1.23456789)

                
        # make sure the shuffled sequence does not lose any elements
        #random.shuffle(self.seq)
        #self.seq.sort()
        #self.assertEqual(self.seq, range(10))

        # should raise an exception for an immutable sequence
        #self.assertRaises(TypeError, random.shuffle, (1,2,3))

    def test_3_attributes_datavectors(self):
        pass
        #element = random.choice(self.seq)
        #self.assertTrue(element in self.seq)

    def test_4_sample(self):
        pass
        #with self.assertRaises(ValueError):
        #    random.sample(self.seq, 20)
        #for element in random.sample(self.seq, 5):
        #    self.assertTrue(element in self.seq)


#if __name__ == '__main__':
#    unittest.main()
 
# automate the process of creating a test suite    
          
suite = unittest.TestLoader().loadTestsFromTestCase(TestStrainDataCPU)
unittest.TextTestRunner(verbosity=2).run(suite)

# (order in which the various test cases will be run is determined by sorting 
# the test function names with the built-in cmp() function)