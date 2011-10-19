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
unittest for the datavector opencl class
"""

from pycbc.datavector.clayer.cpu import real_vector_double_cpu_t
from pycbc.datavector.clayer.cpu import real_vector_single_cpu_t
from pycbc.datavector.clayer.cpu import complex_vector_single_cpu_t
from pycbc.datavector.clayer.cpu import complex_vector_double_cpu_t

from pycbc.datavector.clayer.opencl import real_vector_double_opencl_t
from pycbc.datavector.clayer.opencl import real_vector_single_opencl_t
from pycbc.datavector.clayer.opencl import complex_vector_opencl_cpu_t
from pycbc.datavector.clayer.opencl import complex_vector_opencl_cpu_t

import unittest
import random

class TestDatavectorOpenCl(unittest.TestCase):

    # Test Fixture ---------------------------

    def setUp(self):
        
        self.length = 128 # typ design spec: 4096
        self.delta_x = 1.0
        
        self.dut_real_single= real_vector_single_opencl_t(self.length, 
                                                          self.delta_x)
        self.dut_real_double= real_vector_double_opencl_t(self.length, 
                                                          self.delta_x)
        self.dut_complex_single= complex_vector_single_opencl_t(self.length, 
                                                                self.delta_x)
        self.dut_complex_double= complex_vector_double_opencl_t(self.length, 
                                                                self.delta_x)
        
        #print                                                        
        #print self.dut_real_single
        #print repr(self.dut_real_single)
        #print self.dut_real_double
        #print repr(self.dut_real_double)
        #print self.dut_complex_single
        #print repr(self.dut_complex_single)
        #print self.dut_complex_double
        #print repr(self.dut_complex_double)

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
        find("pycbc.datavector.clayer.opencl.real_vector_single_opencl_t") >= 0,
        " Wrong type of datavector after instanciation")
        
        reference_vector= real_vector_single_cpu_t(self.length, self.delta_x)

        # check correct initialization
        transfer_real_vector_single_to_cpu(reference_vector, 
                                           self.dut_real_single)
        for i in range(self.length):
            self.assertEquals(reference_vector[i], 0.0,
            "datavector not initialized by 0.0 at index: {0}".format(i))

        # data integrity test:
        
        # generate data
        for i in range(self.length):
            reference_vector[i]= random.uniform(-1,1)

        # transfer to GPU
        transfer_real_vector_single_from_cpu(self.dut_real_single, 
                                             reference_vector)
        
        test_vector= real_vector_single_cpu_t(self.length, self.delta_x)

        # fetch from GPU
        transfer_real_vector_single_to_cpu(test_vector, self.dut_real_single)

        # check element by element
        for i in range(self.length):
            self.assertEquals(reference_vector[i], test_vector[i],
            "data integrity failed at index: {0}".
            format(i))

    def test_2_real_double_datavector(self): ###################################
        """
        Test proper datatype, initialization and access of the  
        double precision datavector
        """

        # check type
        self.assertTrue(repr(self.dut_real_double).
        find("pycbc.datavector.clayer.opencl.real_vector_double_opencl_t") >= 0,
        " Wrong type of datavector after instanciation")
        
        reference_vector= real_vector_double_cpu_t(self.length, self.delta_x)

        # check correct initialization
        transfer_real_vector_double_to_cpu(reference_vector, 
                                           self.dut_real_double)
        for i in range(self.length):
            self.assertEquals(reference_vector[i], 0.0,
            "datavector not initialized by 0.0 at index: {0}".format(i))

        # data integrity test:
        
        # generate data
        for i in range(self.length):
            reference_vector[i]= random.uniform(-1,1)

        # transfer to GPU
        transfer_real_vector_double_from_cpu(self.dut_real_double, 
                                             reference_vector)
        
        test_vector= real_vector_double_cpu_t(self.length, self.delta_x)

        # fetch from GPU
        transfer_real_vector_double_to_cpu(test_vector, self.dut_real_double)

        # check element by element
        for i in range(self.length):
            self.assertEquals(reference_vector[i], test_vector[i],
            "data integrity failed at index: {0}".
            format(i))


    def test_3_complex_single_datavector(self): ################################
        """
        Test proper datatype, initialization and access of the  
        single precision complex datavector
        """

        # check type
        self.assertTrue(repr(self.dut_complex_single).
        find("pycbc.datavector.clayer.opencl.complex_vector_single_opencl_t") >= 0,
            " Wrong type of datavector for dut_complex_single")
            
        reference_vector= complex_vector_single_cpu_t(self.length, self.delta_x)

        # check correct initialization
        transfer_complex_vector_single_to_cpu(reference_vector, 
                                              self.dut_complex_single)
        for i in range(self.length):
            self.assertEquals(reference_vector[i], (0+0j),
            "datavector not initialized by 0+0j at index: {0}".format(i))

        # data integrity test:
        
        # generate data
        for i in range(self.length):
            reference_vector[i]= complex(random.uniform(-1,1), 
                                         random.uniform(-1,1))

        # transfer to GPU
        transfer_complex_vector_single_from_cpu(self.dut_complex_single, 
                                                reference_vector)
        
        test_vector= complex_vector_single_cpu_t(self.length, self.delta_x)

        # fetch from GPU
        transfer_complex_vector_single_to_cpu(test_vector, 
                                              self.dut_complex_single)

        # check element by element
        for i in range(self.length):
            self.assertEquals(reference_vector[i], test_vector[i],
            "data integrity failed at index: {0}".
            format(i))
            

    def test_4_complex_double_datavector(self): ################################
        """
        Test proper datatype, initialization and access of the  
        double precision complex datavector
        """

        # check type
        self.assertTrue(repr(self.dut_complex_double).
        find("pycbc.datavector.clayer.opencl.complex_vector_double_opencl_t") >= 0,
            " Wrong type of datavector for dut_complex_double")


        reference_vector= complex_vector_double_cpu_t(self.length, self.delta_x)

        # check correct initialization
        transfer_complex_vector_double_to_cpu(reference_vector, 
                                              self.dut_complex_double)
        for i in range(self.length):
            self.assertEquals(reference_vector[i], (0+0j),
            "datavector not initialized by 0+0j at index: {0}".format(i))

        # data integrity test:
        
        # generate data
        for i in range(self.length):
            reference_vector[i]= complex(random.uniform(-1,1), 
                                         random.uniform(-1,1))

        # transfer to GPU
        transfer_complex_vector_double_from_cpu(self.dut_complex_double, 
                                                reference_vector)
        
        test_vector= complex_vector_double_cpu_t(self.length, self.delta_x)

        # fetch from GPU
        transfer_complex_vector_double_to_cpu(test_vector, 
                                              self.dut_complex_double)

        # check element by element
        for i in range(self.length):
            self.assertEquals(reference_vector[i], test_vector[i],
            "data integrity failed at index: {0}".
            format(i))
            

# automate the process of creating a test suite    
suite = unittest.TestLoader().loadTestsFromTestCase(TestDatavectorOpenCl)
unittest.TextTestRunner(verbosity=2).run(suite)

# (order in which the various test cases will be run is determined by sorting 
# the test function names with the built-in cmp() function)
