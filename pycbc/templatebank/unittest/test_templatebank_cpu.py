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
unittest for the TemplateBank class
"""


from pycbc.clayer.cpu import cpu_context_t

from pycbc.templatebank.cpu import TemplateBankCpu as DUT_TemplateBank

import unittest
import random
import logging

class TestTemplateBankCPU(unittest.TestCase):

    # Test Fixture ---------------------------
    logging.basicConfig(level=logging.DEBUG,
                    format='%(name)s %(asctime)s %(levelname)s %(message)s',
                    filename='test_templatebank_cpu.log',
                    filemode='w')

    logger= logging.getLogger('pycbc.templatebank.unittest.test_templatebank_cpu')

    start_message = 'Starting pycbc templatebank unittest ...'
    logger.debug(start_message)

    def setUp(self):
        
        self.context = cpu_context_t(1)
        
        self.wave_len = 256
        self.wave_dx = 1
        
        self.dut= DUT_TemplateBank(self.context, 
                                   self.wave_len, 
                                   self.wave_dx) #, 
                                   #'H1-TMPLTBANK-871154847-2048.xml.gz')

        print "setup templatebank w/ approximation model: {0}".format(self.dut.approximation_model)

        random.seed()
        
        
    def tearDown(self):
        pass
        
    # Sub tests -------------------------------

    def test_1_init_and_iterator(self): #######################################
        """
        Test proper instanciation and the iterator over templates
        """
        print repr(self.dut)
        
        print repr(self.dut._devicecontext)

        for template in self.dut:
            print repr(template)
            
        print repr(self.dut.pre_condition_vector_t)
        
        #templatebank has instanciated a waveformgenerator
        #which could have generated a precondition vector
        print self.dut.precondition_data
        for i in range(self.wave_len):
            print self.dut.precondition_data[i]

        print "start generate waveform filters"    
        for template in self.dut:
            print template
            htilde = self.dut.waveform_generator.perform_generate_waveform_filter(template)
                
            print htilde
            
        print "end generate waveform filters"
        
        # check just the last one
        for i in range(self.wave_len):
            print htilde[i]


# automate the process of creating a test suite    
suite = unittest.TestLoader().loadTestsFromTestCase(TestTemplateBankCPU)
unittest.TextTestRunner(verbosity=2).run(suite)

# (order in which the various test cases will be run is determined by sorting 
# the test function names with the built-in cmp() function)