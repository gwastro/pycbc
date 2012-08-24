
# Copyright (C) 2012  Alex Nitz
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
import pycbc
import unittest
from pycbc.types import *
from pycbc.scheme import *
from pycbc.filter import *
from pycbc.waveform import *
import pycbc.fft
import numpy 
import base_test

import optparse
from optparse import OptionParser

_parser = OptionParser()

def _check_scheme(option, opt_str, scheme, parser):
    if scheme=='cuda' and not pycbc.HAVE_CUDA:
        raise optparse.OptionValueError("CUDA not found")

    if scheme=='opencl' and not pycbc.HAVE_OPENCL:
        raise optparse.OptionValueError("OpenCL not found")
    setattr (parser.values, option.dest, scheme)

_parser.add_option('--scheme','-s', action='callback', type = 'choice', 
                    choices = ('cpu','cuda','opencl'), 
                    default = 'cpu', dest = 'scheme', callback = _check_scheme,
                    help = 'specifies processing scheme, can be cpu [default], cuda, or opencl')

_parser.add_option('--device-num','-d', action='store', type = 'int', 
                    dest = 'devicenum', default=0,
                    help = 'specifies a GPU device to use for CUDA or OpenCL, 0 by default')

(_opt_list, _args) = _parser.parse_args()

#Changing the optvalues to a dict makes them easier to read
_options = vars(_opt_list)

if _options['scheme'] == 'cpu':
    context = DefaultScheme()
if _options['scheme'] == 'cuda':
    context = CUDAScheme(device_num=_options['devicenum'])
if _options['scheme'] == 'opencl':
    context = OpenCLScheme(device_num=_options['devicenum'])

class TestWaveform(base_test.function_base,unittest.TestCase):
    def setUp(self,*args): 
        self.context = context 

    def test_generation(self):
        with self.context:
            for waveform in ['IMRPhenomB','SpinTaylorT4','SEOBNRv1','EOBNRv2']:
                hc,hp = get_td_waveform(approximant=waveform,mass1=20,mass2=20,delta_t=1.0/4096,f_lower=40)
                self.assertTrue(len(hc)> 0)    
            for waveform in ['TaylorF2', 'IMRPhenomB']:
                htilde = get_fd_waveform(approximant=waveform,mass1=20,mass2=20,delta_f=1.0/256,f_lower=40)     
                self.assertTrue(len(htilde)> 0)                                                         

    def test_errors(self):
        func = get_fd_waveform
        self.assertRaises(ValueError,func,approximant="BLAH")
        self.assertRaises(ValueError,func,approximant="TaylorF2",mass1=3)
        self.assertRaises(ValueError,func,approximant="TaylorF2",mass1=3,mass2=3)
        self.assertRaises(ValueError,func,approximant="TaylorF2",mass1=3,mass2=3,phase_order=7)
        self.assertRaises(ValueError,func,approximant="TaylorF2",mass1=3,mass2=3,phase_order=7)
        self.assertRaises(ValueError,func,approximant="TaylorF2",mass1=3)

        for func in [get_fd_waveform,get_td_waveform]:
            self.assertRaises(ValueError,func,approximant="BLAH")
            self.assertRaises(ValueError,func,approximant="IMRPhenomB",mass1=3)
            self.assertRaises(ValueError,func,approximant="IMRPhenomB",mass1=3,mass2=3)
            self.assertRaises(ValueError,func,approximant="IMRPhenomB",mass1=3,mass2=3,phase_order=7)
            self.assertRaises(ValueError,func,approximant="IMRPhenomB",mass1=3,mass2=3,phase_order=7)
            self.assertRaises(ValueError,func,approximant="IMRPhenomB",mass1=3)

    if _options['scheme'] == 'cuda':
        def test_taylorf2(self): 
            for order in [0,7]:
                h = get_fd_waveform(approximant= "TaylorF2", mass1=1,mass2=1,phase_order=7,amplitude_order=7,delta_f = 1.0/1024,f_lower=15.0) 
                with self.context:
                    s = get_fd_waveform(approximant= "TaylorF2", mass1=1,mass2=1,phase_order=7,amplitude_order=7,delta_f = 1.0/1024,f_lower=15.0)
                    o,i = match(h,s)
                    self.assertAlmostEqual(1,o,places=4)

            for p_order in [0,7]:
                for a_order in [0,7]:
                    h = get_fd_waveform(approximant= "TaylorF2", mass1=1,mass2=1,phase_order=7,amplitude_order=7,delta_f = 1.0/1024,f_lower=15.0) 
                    with self.context:
                        s = get_fd_waveform(approximant= "TaylorF2", mass1=1,mass2=1,phase_order=7,amplitude_order=7,delta_f = 1.0/1024,f_lower=15.0)     
                        diff = ((h-s)).sum()
                    self.assertTrue(diff<10e-32)

    
suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestWaveform))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
        
    NotImpErrors = 0
    for error in results.errors:
        for errormsg in error:
            if type(errormsg) is str:
                if 'NotImplemented' in errormsg:
                    NotImpErrors +=1
                    break
    if results.wasSuccessful():
        sys.exit(0)
    elif len(results.failures)==0 and len(results.errors)==NotImpErrors:
        sys.exit(1)
    else:
        sys.exit(2)
