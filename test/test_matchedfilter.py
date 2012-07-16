
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
import pycbc.fft
import numpy 

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

class TestMatchedFilter(unittest.TestCase):
    def setUp(self):
        from math import sin
        # Use sine wave as test signal
        data = numpy.sin(numpy.arange(0,100,100/(4096.0*64)))
        self.filt = TimeSeries(data,dtype=float32,delta_t=1.0/4096)
        self.filt2 = (self.filt*1)
        self.filt2[0:len(self.filt2)/2].fill(0)  
        self.filt_offset = TimeSeries(numpy.roll(data,4096*32), dtype=float32,
                                      delta_t=1.0/4096)
                                      
        self.filtD = TimeSeries(data,dtype=float64,delta_t=1.0/4096)
        self.filt2D = (self.filtD*1)
        self.filt2D[0:len(self.filt2D)/2].fill(0)  
        self.filt_offsetD = TimeSeries(numpy.roll(data,4096*32), dtype=float64,
                                      delta_t=1.0/4096)                             
         

    def test_ave_snr_noise(self):
        #Test that the average snr in noise is 2
        from numpy.random import normal

        noise = normal(0.0,2,4096*64)
        nplus= TimeSeries(noise,dtype=float32,delta_t=1.0/4096) 
        ntilde = get_frequencyseries(nplus)
        # Calculate a Faux psd for normalization, replace with better algorithm
        psd = (ntilde).squared_norm()  / float(len(nplus)) * nplus.delta_t *2.0

        snr,norm = matchedfilter(self.filt,nplus,psd=psd)   
        ave = snr.squared_norm().sum() * norm * norm /len(snr)
        self.assertAlmostEqual(2,ave,places=5)
        
        
        noise = normal(0.0,2,4096*64)
        nplus= TimeSeries(noise,dtype=float64,delta_t=1.0/4096) 
        ntilde = get_frequencyseries(nplus)
        # Calculate a Faux psd for normalization, replace with better algorithm
        psd = (ntilde).squared_norm()  / float(len(nplus)) * nplus.delta_t *2.0

        snr,norm = matchedfilter(self.filtD,nplus,psd=psd)   
        ave = snr.squared_norm().sum() * norm * norm /len(snr)
        self.assertAlmostEqual(2,ave,places=5)
        
    def test_perfect_match(self):
        o = match(self.filt,self.filt)
        self.assertAlmostEqual(1,o,places=4)
        
        o = match(self.filtD,self.filtD)
        self.assertAlmostEqual(1,o,places=4)
        
    def test_perfect_match_offset(self):
        o = match(self.filt,self.filt_offset)
        self.assertAlmostEqual(1,o,places=4)
        
        o = match(self.filtD,self.filt_offsetD)
        self.assertAlmostEqual(1,o,places=4)
        
    def test_imperfect_match(self):
        f = get_frequencyseries(self.filt)
        f2 = get_frequencyseries(self.filt2)
        o = match(self.filt,self.filt2)
        self.assertAlmostEqual(sqrt(0.5),o,places=3)
        
        f = get_frequencyseries(self.filtD)
        f2 = get_frequencyseries(self.filt2D)
        o = match(self.filtD,self.filt2D)
        self.assertAlmostEqual(sqrt(0.5),o,places=3)
        
    def test_errors(self):

        #Check that an incompatible data and filter produce an error
        self.assertRaises(ValueError,match,self.filt,self.filt[0:5])
        
        #Check that an incompatible psd produces an error
        self.assertRaises(TypeError,match,self.filt,self.filt,psd=self.filt) 
        psd = FrequencySeries(zeros(len(self.filt)/2+1),delta_f=100000)
        self.assertRaises(TypeError,match,self.filt,self.filt,psd=psd)
        
        #Check that only TimeSeries or FrequencySeries are accepted
        self.assertRaises(TypeError,match,zeros(10),zeros(10))


scheme = None

if _options['scheme'] == 'cpu':
    scheme = DefaultScheme()
if _options['scheme'] == 'cuda':
    scheme = CUDAScheme(device_num=_options['devicenum'])
if _options['scheme'] == 'opencl':
    scheme = OpenCLScheme(device_num=_options['devicenum'])
    
suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMatchedFilter))

if __name__ == '__main__':
    with scheme:
        results = unittest.TextTestRunner(verbosity=2).run(suite)
        
    NotImpErrors = 0
    for error in results.errors:
        for errormsg in error:
            if type(errormsg) is str:
                if 'NotImplented' in errormsg:
                    NotImpErrors +=1
                    break
    if results.wasSuccessful():
        sys.exit(0)
    elif len(results.failures)==0 and len(results.errors)==NotImpErrors:
        sys.exit(1)
    else:
        sys.exit(2)
