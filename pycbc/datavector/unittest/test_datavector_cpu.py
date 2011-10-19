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
unittest for the datavector cpu class
"""

from pycbc.datavector.clayer.cpu import real_vector_double_cpu_t
from pycbc.datavector.clayer.cpu import real_vector_single_cpu_t
from pycbc.datavector.clayer.cpu import complex_vector_single_cpu_t
from pycbc.datavector.clayer.cpu import complex_vector_double_cpu_t


import unittest
import random

class TestDatavectorCPU(unittest.TestCase):

    # Test Fixture ---------------------------

    def setUp(self):
        
        self.length = 128 # typ design spec: 4096
        self.delta_x = 1.0
        
        self.dut_real_single= real_vector_single_cpu_t(self.length, 
                                                       self.delta_x)
        self.dut_real_double= real_vector_double_cpu_t(self.length, 
                                                       self.delta_x)
        self.dut_complex_single= complex_vector_single_cpu_t(self.length, 
                                                             self.delta_x)
        self.dut_complex_double= complex_vector_double_cpu_t(self.length, 
                                                             self.delta_x)

        random.seed()
        self.digits_to_check_single_against_double = 7

    def tearDown(self):
        pass
        
    # Sub tests -------------------------------

    def test_1_real_single_datavector(self): ###################################
        """
        Test proper datatype, initialization and access of the  
        single precision datavector
        """

        # check type
        self.assertTrue(repr(self.dut_real_single).
        find("pycbc.datavector.clayer.cpu.real_vector_single_cpu_t") >= 0,
        " Wrong type of datavector after instanciation")

        # check correct initialization
        for i in range(self.length):
            self.assertEquals(self.dut_real_single[i], 0.0,
            "datavector not initialized by 0.0 at index: {0}".format(i))
        
        # check datavectors for throwing exceptions if 
        # trying to access out of bounds
        def out_of_bounds_access(self, vector_to_test):
            vector_to_test[self.length] = 2.0
            tmp = vector_to_test[self.length]
        
        self.assertRaises(ValueError, out_of_bounds_access, self, 
                          self.dut_real_single)

        self.assertRaises(ValueError, self.dut_real_single.set_start, 
                          self.length)

        # access test
        for i in range(self.length):
            tmp= random.uniform(-1,1)
            self.dut_real_single[i] = tmp
            self.assertAlmostEquals(self.dut_real_single[i], tmp, 
            self.digits_to_check_single_against_double,
            "datavector access test failed at index: {0}".format(i))

    def test_2_real_double_datavector(self): ###################################
        """
        Test proper datatype, initialization and access of the  
        double precision datavector
        """

        # check type
        self.assertTrue(repr(self.dut_real_double).
        find("pycbc.datavector.clayer.cpu.real_vector_double_cpu_t") >= 0,
        " Wrong type of datavector after instanciation")

        # check correct initialization
        for i in range(self.length):
            self.assertEquals(self.dut_real_double[i], 0.0,
            "datavector not initialized by 0.0 at index: {0}".format(i))
        
        # check datavectors for throwing exceptions if 
        # trying to access out of bounds
        def out_of_bounds_access(self, vector_to_test):
            vector_to_test[self.length] = 2.0
            tmp = vector_to_test[self.length]
        
        self.assertRaises(ValueError, out_of_bounds_access, self, 
                          self.dut_real_double)

        self.assertRaises(ValueError, self.dut_real_double.set_start, 
                          self.length)

        # access test
        for i in range(self.length):
            tmp= random.uniform(-1,1)
            self.dut_real_double[i] = tmp
            self.assertEquals(self.dut_real_double[i], tmp, 
            "datavector access test failed at index: {0}".format(i))

    def test_3_complex_single_datavector(self): ################################
        """
        Test proper datatype, initialization and access of the  
        single precision complex datavector
        """

        # check type
        self.assertTrue(repr(self.dut_complex_single).
        find("pycbc.datavector.clayer.cpu.complex_vector_single_cpu_t") >= 0,
            " Wrong type of datavector for dut_complex_single")
            
        # check correct initialization
        for i in range(self.length):
            self.assertEquals(self.dut_complex_single[i], (0+0j),
            "dut_complex_single not initialized by (0+0j) at index: {0}".
            format(i))

        # check datavectors for throwing exceptions if 
        # trying to access out of bounds
        def out_of_bounds_access(self, vector_to_test):
            vector_to_test[self.length] = (2.0+3.0j)
            tmp = vector_to_test[self.length]
        
            self.assertRaises(ValueError, out_of_bounds_access, self, 
                                self.dut_complex_single)

            self.assertRaises(ValueError, self.dut_complex_single.set_start, 
                                self.length)
       
        # access test
        for i in range(self.length):
                tmp= complex(random.uniform(-1,1), random.uniform(-1,1))
                self.dut_complex_single[i] = tmp
                self.assertAlmostEquals(self.dut_complex_single[i], tmp, 
                self.digits_to_check_single_against_double,
                "dut_complex_single access test failed at index: {0}".format(i))

    def test_4_complex_double_datavector(self): ################################
        """
        Test proper datatype, initialization and access of the  
        double precision complex datavector
        """

        # check type
        self.assertTrue(repr(self.dut_complex_double).
        find("pycbc.datavector.clayer.cpu.complex_vector_double_cpu_t") >= 0,
            " Wrong type of datavector for dut_complex_double")
            
        # check correct initialization
        for i in range(self.length):
            self.assertEquals(self.dut_complex_double[i], (0+0j),
            "dut_complex_double not initialized by (0+0j) at index: {0}".
            format(i))

        # check datavectors for throwing exceptions if 
        # trying to access out of bounds
        def out_of_bounds_access(self, vector_to_test):
            vector_to_test[self.length] = (2.0+3.0j)
            tmp = vector_to_test[self.length]
        
            self.assertRaises(ValueError, out_of_bounds_access, self, 
                                self.dut_complex_double)

            self.assertRaises(ValueError, self.dut_complex_double.set_start, 
                                self.length)
       
        # access test
        for i in range(self.length):
                tmp= complex(random.uniform(-1,1), random.uniform(-1,1))
                self.dut_complex_double[i] = tmp
                self.assertEquals(self.dut_complex_double[i], tmp, 
                "dut_complex_double access test failed at index: {0}".format(i))

# automate the process of creating a test suite    
suite = unittest.TestLoader().loadTestsFromTestCase(TestDatavectorCPU)
unittest.TextTestRunner(verbosity=2).run(suite)

# (order in which the various test cases will be run is determined by sorting 
# the test function names with the built-in cmp() function)
