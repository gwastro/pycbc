
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
import unittest
from pycbc.types import *
from pycbc.filter import *
from pycbc.scheme import *
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

class TestUtils(unittest.TestCase):
    def setUp(self,*args): 
        self.context = context
        self.delta_t = 1.0 / 4096
        self.target_delta_t = 1.0 / 1024
        self.a = TimeSeries([1,2,3,4], delta_t=self.delta_t, dtype=float32)
        self.b = TimeSeries([1,2,3,4], delta_t=self.delta_t, dtype=float64)
        self.c = TimeSeries([1,2,3,4], delta_t=self.delta_t, dtype=complex64)
        self.d = Array([1,2,3,4], dtype=float32)          

    if type(context) is DefaultScheme:
        def test_resample_float32(self):
            ra = resample_to_delta_t(self.a, self.target_delta_t)
            self.assertAlmostEqual(ra[0], 0.00696246)
            ra = resample_to_delta_t(self.a, self.delta_t)
            self.assertAlmostEqual(ra[0], 1)
        
        def test_resample_float64(self):
            rb = resample_to_delta_t(self.b, self.target_delta_t)
            self.assertAlmostEqual(rb[0], 0.00696246)
            rb = resample_to_delta_t(self.b, self.delta_t)
            self.assertAlmostEqual(rb[0], 1)

    def test_resample_errors(self):
        self.assertRaises(TypeError, resample_to_delta_t, self.c, self.target_delta_t)
        self.assertRaises(TypeError, resample_to_delta_t, self.d, self.target_delta_t)

        if type(self.context) is CUDAScheme:
            with self.context:
                self.assertRaises(TypeError, resample_to_delta_t, self.a, self.target_delta_t)

                                                      
    
suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestUtils))

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
