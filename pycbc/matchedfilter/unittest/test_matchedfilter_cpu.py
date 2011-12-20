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
# ==============================================================================
#
#                                   Preamble
#
# ==============================================================================
#

"""
unittest for the MatchedFilterBase class
"""


from pycbc.clayer.pycbccpu import cpu_context_t
from pycbc.datavector.cpu import complex_vector_single_cpu_t
from pycbc.matchedfilter.cpu import MatchedFilterCpu as DUT_MatchedFilterCpu

import unittest
import random
import logging

class TestMatchedFilterCPU(unittest.TestCase):

    # Test Fixture ---------------------------
    logging.basicConfig(level=logging.DEBUG,
                    format='%(name)s %(asctime)s %(levelname)s %(message)s',
                    filename='test_matchedfilter_cpu.log',
                    filemode='w')

    logger= logging.getLogger('pycbc.matchedfilter.unittest.test_matchedfilter_cpu')

    start_message = 'Starting pycbc matchedfilter unittest ...'
    logger.debug(start_message)

    def setUp(self):
        
        self.context = cpu_context_t(1)
        
        self.strain_len = 10
        self.strain_dx = 1

        self.stilde = complex_vector_single_cpu_t(self.context, self.strain_len,
                                                                self.strain_dx)

        self.htilde = complex_vector_single_cpu_t(self.context, self.strain_len,
                                                                self.strain_dx)
                                                                
        self.dut= DUT_MatchedFilterCpu(self.context, self.strain_len, 
                                                     self.strain_dx)

        print "setup matchedfilter w/ len={0}, dx={1}".format(self.strain_len, 
                                                              self.strain_dx)

        random.seed()
        
    def tearDown(self):
        pass
        
    # Sub tests -------------------------------

    def test_1_init(self): #######################################
        """
        Testing initalization 
        """
        print repr(self.dut)
        #print repr(self.dut._devicecontext)

        for i in range(self.strain_len):
            self.stilde[i]= random.uniform(-1,1)
            self.htilde[i]= random.uniform(-1,1)

        print self.stilde
        for i in range(self.strain_len):
            print self.stilde[i]

        print self.htilde
        for i in range(self.strain_len):
            print self.htilde[i]
        
        snr = self.dut.perform_generate_snr(self.stilde, self.htilde) # mf perform!
        
        print snr
        for i in range(self.strain_len):
            print snr[i]
            
        print "found maximum snr: {0}".format(self.dut.perform_find_max(snr))

        


# automate the process of creating a test suite    
suite = unittest.TestLoader().loadTestsFromTestCase(TestMatchedFilterCPU)
unittest.TextTestRunner(verbosity=2).run(suite)

# (order in which the various test cases will be run is determined by sorting 
# the test function names with the built-in cmp() function)
