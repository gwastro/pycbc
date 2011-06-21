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

from straindata_cpu import StrainDataCpu as DUT_StrainData

import unittest
import random

class TestStrainDataCPU(unittest.TestCase):

    # Test Fixture ---------------------------

    def setUp(self):
        
        self.search_time = 2048
        self.sample_freq = 4096
        self.length =      self.search_time * self.sample_freq
        self.segments = 15
        self.gps_start_time= 871147532
        self.gps_end_time= self.gps_start_time + self.search_time
        self.interferometer = "H1"
        
        self.dut= DUT_StrainData(self.gps_start_time, self.gps_end_time, 
                                 self.segments, self.sample_freq, 
                                 self.interferometer)


        print
        print "setup strain w/ segment length: {0}".format(self.dut.segments_length)

        random.seed()
        self.digits_to_check_single_against_double = 7

    def tearDown(self):
        pass
        
    # Sub tests -------------------------------

    def test_1_initial_timeseries(self): #######################################
        """
        Test proper datatype, initialization and access of the initial 
        double precision time series of strain data s(t)
        """

        # check type
        self.assertTrue(repr(self.dut.time_series).
        find("datavectorcpu.real_vector_double_t") >= 0,
        " Wrong type of datavector for straindata at the initial phase")


        # check correct initialization
        for i in range(self.length):
            self.assertEquals(self.dut.time_series[i], 0.0,
            "time_series not initialized by 0.0 at index: {0}".format(i))
        
        # check datavectors for throwing exceptions if 
        # trying to access out of bounds
        def out_of_bounds_access(self, vector_to_test):
            vector_to_test[self.length] = 2.0
            tmp = vector_to_test[self.length]
        
        self.assertRaises(ValueError, out_of_bounds_access, self, 
                          self.dut.time_series)

        self.assertRaises(ValueError, self.dut.time_series.set_start, 
                          self.length)
        
        #try:
        #    self.dut.time_series[self.length] = 2.0
        #except ValueError:
        #    print 'Catched ValueError'
            
        # access test
        for i in range(self.length):
            tmp= random.uniform(-1,1)
            self.dut.time_series[i] = tmp
            self.assertEquals(self.dut.time_series[i], tmp, 
            "time_series access test failed at index: {0}".format(i))

    def test_2_convert_to_single_preci(self): ##################################
        """
        Test the conversion of the initial 
        double precision time series of strain data s(t) to a single precision
        time series
        """
        
        # initialize datavector by random numbers and save in temporarily list 
        tmp_series = []
        for i in range(self.length):
            tmp = random.uniform(-1,1)
            tmp_series.append(tmp)
            self.dut.time_series[i] = tmp
        
        self.dut.convert_to_single_preci()
        
        # check type (straindata converted to single precision)
        self.assertTrue(repr(self.dut.time_series).
        find("datavectorcpu.real_vector_single_t") >= 0,
        " Wrong type of datavector for straindata after convert() call")
        
        # check data integrity after conversion
        for i in range(self.length):
                self.assertAlmostEquals(self.dut.time_series[i], 
                                        tmp_series[i], 
                self.digits_to_check_single_against_double,
                "straindata access test failed after convert() at index: {0}".
                format(i))

    def test_3_segmenting(self): ###############################################
        """
        Test proper datatype, initialization and access of the segmented 
        single precision frequency series of strain data stilde(f)
        """

        # iterate over segments
        for stilde in self.dut:
            # check type
            self.assertTrue(repr(stilde).
            find("datavectorcpu.complex_vector_single_t") >= 0,
            " Wrong type of datavector for stilde")
            
            # check correct initialization
            for i in range(self.dut.segments_length):
                self.assertEquals(stilde[i], (0+0j),
                "frequency_series not initialized by (0+0j) at index: {0}".
                format(i))

            # check datavectors for throwing exceptions if 
            # trying to access out of bounds
            def out_of_bounds_access(self, vector_to_test):
                vector_to_test[self.length] = 2.0
                tmp = vector_to_test[self.length]
        
            self.assertRaises(ValueError, out_of_bounds_access, self, 
                                stilde)

            self.assertRaises(ValueError, stilde.set_start, 
                                self.length)
       
            # access test
            for i in range(self.dut.segments_length):
                tmp= complex(random.uniform(-1,1), random.uniform(-1,1))
                stilde[i] = tmp
                self.assertAlmostEquals(stilde[i], tmp, 
                self.digits_to_check_single_against_double,
                "stilde access test failed at index: {0}".format(i))

        # 2nd run to check if the iterator was resetted correctly
        for stilde in self.dut:
                   
            for i in range(self.dut.segments_length):
                tmp= complex(random.uniform(-1,1), random.uniform(-1,1))
                stilde[i] = tmp
                self.assertAlmostEquals(stilde[i], tmp, 
                self.digits_to_check_single_against_double,
                "stilde access test failed at index: {0}".format(i))
        

# automate the process of creating a test suite    
suite = unittest.TestLoader().loadTestsFromTestCase(TestStrainDataCPU)
unittest.TextTestRunner(verbosity=2).run(suite)

# (order in which the various test cases will be run is determined by sorting 
# the test function names with the built-in cmp() function)