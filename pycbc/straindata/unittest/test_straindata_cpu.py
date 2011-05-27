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

# preliminary hard coded path to packages 
import sys
sys.path.append('/Users/kawies/dev/src/pycbc')
sys.path.append('/Users/kawies/dev/src/pycbc/pycbc/straindata')

from straindata_cpu import StrainData as DUT_StrainData
import unittest
import random

class TestStrainDataCPU(unittest.TestCase):

    # Test Fixture ---------------------------

    def setUp(self):
        
        self.length =   5
        self.segments = 3
        self.interferometer = "H1"
        self.dut= DUT_StrainData(self.segments, self.length, self.interferometer)

        random.seed()
        self.digits_to_check_single_against_double = 7

    def tearDown(self):
        pass
        
    # Sub tests -------------------------------

    def test_1_initial_timeseries(self):
        """
        Test proper datatype, initialization and access of the initial 
        double precision time series of strain data s(t)
        """

        # check type
        self.assertTrue(repr(self.dut.strain_time_series).find("datavectorcpu.real_vector_double_t") >= 0,
        " Wrong type of datavector for straindata at the initial phase")

        # check correct initialization
        for i in range(self.length):
            self.assertEquals(self.dut.strain_time_series[i], 0.0,
            "strain_time_series not initialized by 0.0 at index: {0}".format(i))
        
        # access test
        for i in range(self.length):
            tmp= random.uniform(-1,1)
            self.dut.strain_time_series[i] = tmp
            self.assertEquals(self.dut.strain_time_series[i], tmp,
            "strain_time_series access test failed at index: {0}".format(i))

    def test_2_convert(self):
        """
        Test the conversion of the initial 
        double precision time series of strain data s(t) to a single precision
        time series
        """
        pass


    def test_3_segmenting(self):
        """
        Test proper datatype, initialization and access of the segmented 
        single precision frequency series of strain data stilde(f)
        """

        # iterate over segments
        for stilde in self.dut:
            # check type
            self.assertTrue(repr(stilde).find("datavectorcpu.complex_vector_single_t") >= 0,
            " Wrong type of datavector for stilde")

            # check correct initialization
            for i in range(self.length):
                self.assertEquals(stilde[i], 0.0,
                "strain_freq_series not initialized by 0.0 at index: {0}".format(i))
        
            # access test
            for j in range(self.length * 2): # complex!
                tmp= random.uniform(-1,1)
                stilde[j] = tmp
                self.assertAlmostEquals(stilde[j], tmp, 
                self.digits_to_check_single_against_double,
                "stilde access test failed at index: {0}".format(j))

        # 2nd run to check if the iterator was resetted correctly
        for stilde in self.dut:
                   
            for j in range(self.length * 2): # complex!
                tmp= random.uniform(-1,1)
                stilde[j] = tmp
                self.assertAlmostEquals(stilde[j], tmp, 
                self.digits_to_check_single_against_double,
                "stilde access test failed at index: {0}".format(j))
        
        #element = random.choice(self.seq)
        #self.assertTrue(element in self.seq)
        
        #with self.assertRaises(ValueError):
        #    random.sample(self.seq, 20)
        #for element in random.sample(self.seq, 5):
        #    self.assertTrue(element in self.seq)

        # make sure the shuffled sequence does not lose any elements
        #random.shuffle(self.seq)
        #self.seq.sort()
        #self.assertEqual(self.seq, range(10))

        # should raise an exception for an immutable sequence
        #self.assertRaises(TypeError, random.shuffle, (1,2,3))


#if __name__ == '__main__':
#    unittest.main()
 
# automate the process of creating a test suite    
          
suite = unittest.TestLoader().loadTestsFromTestCase(TestStrainDataCPU)
unittest.TextTestRunner(verbosity=2).run(suite)

# (order in which the various test cases will be run is determined by sorting 
# the test function names with the built-in cmp() function)