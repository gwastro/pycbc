# Copyright (C) 2012  Andrew Miller, Alex Nitz
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
'''
These are the unittests for the pycbc frame/cache reading functions
'''


import pycbc
import unittest
from pycbc.types import *
from pycbc.scheme import *
import pycbc.frame
import numpy
import base_test
import tempfile
import lal
from pylal import frutils, Fr
import os
import subprocess
import sys

import optparse
from optparse import OptionParser

_parser = OptionParser()

def _check_scheme(option, opt_str, scheme, parser):
    if scheme=='cuda' and not pycbc.HAVE_CUDA:
        raise optparse.OptionValueError("CUDA not found")

    if scheme=='opencl' and not pycbc.HAVE_OPENCL:
        raise optparse.OptionValueError("OpenCL not found")
    setattr (parser.values, option.dest, scheme)

_parser.add_option('--scheme','-s', action='callback', type = 'choice', choices = ('cpu','cuda','opencl'), 
                    default = 'cpu', dest = 'scheme', callback = _check_scheme,
                    help = 'specifies processing scheme, can be cpu [default], cuda, or opencl')

_parser.add_option('--device-num','-d', action='store', type = 'int', dest = 'devicenum', default=0,
                    help = 'specifies a GPU device to use for CUDA or OpenCL, 0 by default')

(_opt_list, _args) = _parser.parse_args()

#Changing the optvalues to a dict makes them easier to read
_options = vars(_opt_list)

# ********************GENERIC ARRAY TESTS ***********************

# Although these tests don't require CPU/GPU swapping of inputs, we will
# still inherit from the function_base so that we can use the check function.
class FrameTestBase(base_test.function_base):

    def checkCurrentState(self, inputs, results, places):
        super(FrameTestBase,self).checkCurrentState(inputs, results, places)
        for a in inputs:
            if isinstance(a,pycbc.types.Array):
                self.assertEqual(a.delta_t, self.delta_t)
                self.assertTrue(a.dtype == self.dtype)

    def setUp(self):
        if _options['scheme'] == 'cpu':
            self.scheme = type(None)
        elif _options['scheme'] == 'cuda':
            self.scheme = pycbc.scheme.CUDAScheme
        elif _options['scheme'] == 'opencl':
            self.scheme = pycbc.scheme.OpenCLScheme 
        # Number of decimal places to compare for single precision
        if self.dtype == numpy.float32 or self.dtype == numpy.complex64:
            self.places = 5
        # Number of decimal places to compare for double precision
        else:
            self.places = 13
            
        self.size = pow(2,12)

        if self.dtype == numpy.float32 or self.dtype == numpy.float64:
            self.kind = 'real'
        else:
            self.kind = 'complex'

        self.data1 = numpy.array(numpy.random.rand(self.size), dtype=self.dtype)
        self.data2 = numpy.array(numpy.random.rand(self.size), dtype=self.dtype)
        
        # If the dtype is complex, we should throw in some complex values as well
        if self.kind == 'complex':
            self.data1 += numpy.random.rand(self.size) * 1j
            self.data2 += numpy.random.rand(self.size) * 1j        
            
        self.delta_t = .5
        self.epoch = lal.LIGOTimeGPS(123456,0)
        
    def test_frame(self):
        
        # This is a file in the temp directory that will be deleted when it is garbage collected
        frmfile = tempfile.NamedTemporaryFile()  
        filename = frmfile.name + ".gwf"
        
        # Now we will create a frame file, specifiying that it is a timeseries
        Fr.frputvect(filename,[{'name':'channel1', 'data':self.data1, 'start':int(self.epoch), 'dx':self.delta_t,'type':1},
                                {'name':'channel2', 'data':self.data2, 'start':int(self.epoch), 'dx':self.delta_t,'type':1}])
        
        with self.context:
            if _options['scheme'] == 'cpu':
                # Reading just one channel first
                ts1 = pycbc.frame.read_frame(filename,'channel1')
                # Chacking all values
                self.checkCurrentState((ts1,),(self.data1,),self.places)
                # Now checking the start time
                self.assertTrue(ts1.start_time == self.epoch)
                # And the duration
                self.assertTrue(ts1.end_time-ts1.start_time == self.size*self.delta_t)
                
                # Now reading multiple channels
                ts2 = pycbc.frame.read_frame(filename,['channel1','channel2'])
                # We should get back a list
                self.assertTrue(type(ts2) is list)
                self.checkCurrentState(ts2, (self.data1,self.data2), self.places)
                self.assertTrue(ts2[0].start_time == self.epoch)
                self.assertTrue(ts2[1].start_time == self.epoch)
                self.assertTrue(ts2[0].end_time-ts2[0].start_time == self.size*self.delta_t)
                self.assertTrue(ts2[1].end_time-ts2[1].start_time == self.size*self.delta_t)
                
                # These are the times and indices for the segment we will try to read
                start = self.epoch+10
                end = self.epoch+50
                startind = int(10/self.delta_t)
                endind = int(50/self.delta_t)
                
                # Now reading in a specific segment with an integer
                ts3 = pycbc.frame.read_frame(filename, 'channel1',
                                             start_time=int(start),
                                             end_time=int(end))
                
                # The same, but with a LIGOTimeGPS for the start and end times
                ts4 = pycbc.frame.read_frame(filename, 'channel1',
                                             start_time=start, end_time=end)

                # Now we will check those two TimeSeries
                self.checkCurrentState((ts3,ts4), (self.data1[startind:endind],self.data1[startind:endind]), self.places)
                self.assertTrue(40 - (float(ts3.end_time)-float(ts3.start_time)) < self.delta_t)
                self.assertTrue(ts3.start_time == start)
                
                self.assertTrue(40 - (float(ts4.end_time)-float(ts4.start_time)) < self.delta_t)
                self.assertTrue(ts4.start_time == start)

                # And now some cases that should raise errors

                # There must be a span grater than 0
                self.assertRaises(ValueError, pycbc.frame.read_frame, filename,
                                  'channel1', start_time=self.epoch,
                                  end_time=self.epoch)
                # The start must be before the end
                self.assertRaises(ValueError, pycbc.frame.read_frame, filename,
                                  'channel1', start_time=self.epoch+1,
                                  end_time=self.epoch)
                # Non integer times should also raise an error
                badtime = lal.LIGOTimeGPS(int(self.epoch)+5,1000)
                
                self.assertRaises(ValueError, pycbc.frame.read_frame, filename,
                                  'channel1', start_time=self.epoch,
                                  end_time=badtime)
                self.assertRaises(ValueError, pycbc.frame.read_frame, filename,
                                  'channel1', start_time=float(self.epoch),
                                  end_time=float(badtime))
    

if _options['scheme']=='cpu':
    context = pycbc.scheme.CPUScheme()
elif _options['scheme']=='cuda':
    context = pycbc.scheme.CUDAScheme(device_num=_options['devicenum'])
elif _options['scheme']=='opencl':
    context = pycbc.scheme.OpenCLScheme(device_num=_options['devicenum'])
        
TestClasses = []
types = [numpy.float32, numpy.float64, numpy.complex64, numpy.complex128]

for ty in types:
    klass = type('{0}_{1}_Test'.format(_options['scheme'],ty.__name__),
                    (FrameTestBase,unittest.TestCase),
                    {'dtype': ty, 'context' : context})
    TestClasses.append(klass)

if __name__ == '__main__':
    suite = unittest.TestSuite()
    for klass in TestClasses:
        suite.addTest(unittest.TestLoader().loadTestsFromTestCase(klass))
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
        
