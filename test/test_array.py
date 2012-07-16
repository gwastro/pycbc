# Copyright (C) 2012  Alex Nitz, Andrew Miller
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
These are the unittests for the pycbc array type
'''


import pycbc
import unittest
from pycbc.types import *
from pycbc.scheme import *
import numpy
import base_test
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

# By importing the current schemes array type, it will make it easier to check the  array types later
if _options['scheme'] == 'cuda':
    import pycuda
    import pycuda.gpuarray
    from pycuda.gpuarray import GPUArray as SchemeArray
elif _options['scheme'] == 'opencl':
    import pyopencl
    import pyopencl.array
    from pyopencl.array import Array as SchemeArray
elif _options['scheme'] == 'cpu':
    from numpy import ndarray as SchemeArray


class ArrayTestBase(base_test.array_base):
    def setUp(self):
    
        # We need to check for correct creation from all dtypes, 
        # and errors from incorrect operations so the other precision of 
        # odtype needs to be available as well
        self.other_precision = {numpy.complex64 : numpy.complex128,
                            numpy.complex128 : numpy.complex64,
                            numpy.float32 : numpy.float64,
                            numpy.float64 : numpy.float32}
        
        # Number of decimal places to compare for single precision
        if self.dtype == numpy.float32 or self.dtype == numpy.complex64:
            self.places = 5
        # Number of decimal places to compare for double precision
        else:
            self.places = 13
            
        # We will also need to check whether dtype and odtype are real or complex,
        # so that we can test non-zero imaginary parts. 
        if self.dtype == numpy.float32 or self.dtype == numpy.float64:
            self.kind = 'real'
        else:
            self.kind = 'complex'
        if self.odtype == numpy.float32 or self.odtype == numpy.float64:
            self.okind = 'real'
        else:
            self.okind = 'complex'
        # Now that the kinds are set, we need to call our parent method to set up all the
        # inputs and answers for our functions
        self.setNumbers()
        
        # Here two test arrays are created.
        # We will use multiples of each to check things are on the cpu/gpu when they should be.
        self.a1 = Array(self.alist, dtype=self.dtype)
        self.a2 = Array(self.alist, dtype=self.dtype)
        self.a3 = Array(self.alist, dtype=self.dtype)

        self.b1 = Array(self.blist, dtype=self.odtype)
        self.b2 = Array(self.blist, dtype=self.odtype)
        
        # And now we'll make a copy of the scalar. We make a copy so we can still check it later.
        self.s = self.scalar

        # Finally, we want to have an array that we shouldn't be able to operate on,
        # because the precision is wrong, and one where the length is wrong.
        self.bad = Array([1,1,1],dtype = self.other_precision[self.odtype])
        self.bad2 = Array([1,1,1,1], dtype = self.dtype)

    def test_set(self):
        c = self.a1 * 1
        with self.context:
            # First we will check that get works properly for all
            # the different python syntaxes
            self.assertTrue(self.a1[:][0] == self.alist[0:3][0])
            self.assertTrue(self.a1[:][1] == self.alist[0:3][1])
            self.assertTrue(self.a1[:][2] == self.alist[0:3][2])
            self.assertRaises(IndexError,self.a1[:].__getitem__,3)
            self.assertTrue(self.a1[-1] ==self.alist[2])
            self.assertTrue(self.a1[-2] == self.alist[1])
            self.assertTrue(self.a1[1:2][0] == self.alist[1])
            self.assertRaises(IndexError,self.a1[1:2].__getitem__,1)
            self.assertTrue(self.a1[:-1][0] == self.alist[0:2][0])
            self.assertTrue(self.a1[:-1][1] == self.alist[0:2][1])
            self.assertTrue(self.a1[-1:][0] == self.alist[2])
            self.assertRaises(IndexError, self.a1.__getitem__, 3)
            self.assertRaises(IndexError, self.a1.__getitem__, -4)
                            
        if not (self.kind == 'real' and self.okind == 'complex'):   
            with self.context:
                # We will check setting from arrays on multiple contexts
                self.b1 *= 1
                c[0] = Array(self.b1[0])
                c[1] = Array(self.b2[1])
                c[2] = Array(self.b1[2])
                self.checkCurrentState((self.b1, self.b2, c),(self.blist,self.blist,self.blist), self.places)
                c = self.a1 * 1
            # And also going back to the CPU from Other
            c[0] = Array(self.b1[0])
            c[1] = Array(self.b2[1])
            c[2] = Array(self.b1[2])
            self.checkCurrentState((self.b1, self.b2, c),(self.blist,self.blist,self.blist), self.places)
                
        else:
            with self.context:
                self.assertRaises(ValueError, self.a1.__setitem__, 0, Array(self.b1[0]))

    def test_numpy_init(self):
        with self.context:                                
            in1 = numpy.array([5,3,1],dtype=self.odtype)
            in2 = numpy.array([5,3,1],dtype=self.other_precision[self.odtype])
            
            #We don't want to cast complex as real
            if not (self.kind == 'real' and self.okind == 'complex'):
                #First we must check that the dtype is correct when specified
                out1 = Array(in1, dtype=self.dtype)
                out2 = Array(in2, dtype=self.dtype)
                #to be sure that it is copied
                in1 += 1
                in2 += 1
                self.assertTrue(type(out1._scheme) == self.scheme)
                self.assertTrue(type(out1._data) is SchemeArray)
                self.assertEqual(out1[0],5)
                self.assertEqual(out1[1],3)
                self.assertEqual(out1[2],1)
                self.assertTrue(out1.dtype==self.dtype)
                
                
                self.assertTrue(type(out2._scheme) == self.scheme)
                self.assertTrue(type(out2._data) is SchemeArray)
                self.assertEqual(out2[0],5)
                self.assertEqual(out2[1],3)
                self.assertEqual(out2[2],1)
                self.assertTrue(out2.dtype==self.dtype)
                
                in1-=1
                in2-=1
            
            # Also, when it is unspecified
            out3 = Array(in1)
            in1 += 1
            self.assertTrue(type(out3._scheme) == self.scheme)
            self.assertTrue(type(out3._data) is SchemeArray)
            self.assertEqual(out3[0],5)
            self.assertEqual(out3[1],3)
            self.assertEqual(out3[2],1)
            self.assertTrue(out3.dtype==self.odtype)
                        
            # Check for copy=false
            # On the CPU, this should be possible
            in3 = numpy.array([5,3,1],dtype=self.dtype)
            if _options['scheme'] == 'cpu':
                out4 = Array(in3,copy=False)
                in3 += 1
                
                self.assertTrue(out4.dtype==self.dtype)
                self.assertTrue(type(out4._scheme) == self.scheme)
                self.assertEqual(out4[0],6)
                self.assertEqual(out4[1],4)
                self.assertEqual(out4[2],2)
                
            # If we're in different scheme, this should raise an error
            else:
                self.assertRaises(TypeError, Array, in3, copy=False)

            # We also need to check initialization using GPU arrays
            if _options['scheme'] == 'cuda':
                in4 = pycuda.gpuarray.zeros(3,self.dtype)
            elif _options['scheme'] == 'opencl':
                in4 = pyopencl.array.zeros(pycbc.scheme.mgr.state.queue,3, self.dtype)
            if _options['scheme'] != 'cpu':
                out4 = Array(in4, copy=False)
                in4 += 1
                self.assertTrue(type(out4._scheme) == self.scheme)
                self.assertTrue(type(out4._data) is SchemeArray)
                self.assertEqual(out4[0],1)
                self.assertEqual(out4[1],1)
                self.assertEqual(out4[2],1)
                self.assertTrue(out4.dtype==self.dtype)
                            
            # We should be able to create an array from the wrong dtype, and
            # it should be cast as float64
            in5 = numpy.array([1,2,3],dtype=numpy.int32)
            out5 = Array(in5)
            in5 += 1
            self.assertTrue(type(out5._scheme) == self.scheme)
            self.assertTrue(type(out5._data) is SchemeArray)
            self.assertEqual(out5[0],1)
            self.assertEqual(out5[1],2)
            self.assertEqual(out5[2],3)
            self.assertTrue(out5.dtype==numpy.float64)
            
            # We shouldn't be able to copy it though
            self.assertRaises(TypeError,Array,in5, copy=False)
            
            # Just checking that we can make an empty array correctly
            empty = numpy.array([])
            out6 = Array(empty)
            self.assertTrue(out6.dtype==numpy.float64)
            self.assertRaises(IndexError, out6.__getitem__,0)
            
        if _options['scheme'] != 'cpu':
            self.assertRaises(TypeError, Array, in4, copy=False)
        self.assertRaises(TypeError, Array, in5, copy=False)
                    

    def test_array_init(self):
        # this array is made outside the context so we can check that an error is raised when copy = false in a GPU scheme
        cpuarray = Array([1,2,3])
        with self.context:      
            in1 = Array([5,3,1],dtype=self.odtype)
            in2 = Array([5,3,1],dtype=self.other_precision[self.odtype])
            self.assertTrue(type(in1._scheme) == self.scheme)
            self.assertTrue(type(in1._data) is SchemeArray)
            self.assertTrue(type(in2._scheme) == self.scheme)
            self.assertTrue(type(in2._data) is SchemeArray)
            # We don't want to cast complex as real
            if not (self.kind=='real' and self.okind == 'complex'):
                # First we must check that the dtype is correct when specified
                out1 = Array(in1, dtype=self.dtype)
                out2 = Array(in2, dtype=self.dtype)
                # to be sure that it is copied
                in1 += 1
                in2 += 1
                
                self.assertTrue(type(out1._scheme) == self.scheme)
                self.assertTrue(type(out1._data) is SchemeArray)
                self.assertEqual(out1[0],5)
                self.assertEqual(out1[1],3)
                self.assertEqual(out1[2],1)
                self.assertTrue(out1.dtype==self.dtype)
                                
                if out1.dtype == numpy.float32:
                    self.assertTrue(out1.precision == 'single')
                    #self.assertTrue(out1.kind == 'real')
                if out1.dtype == numpy.float64:
                    self.assertTrue(out1.precision == 'double')
                    #self.assertTrue(out1.kind == 'real')
                if out1.dtype == numpy.complex64:
                    self.assertTrue(out1.precision == 'single')
                    #self.assertTrue(out1.kind == 'complex')
                if out1.dtype == numpy.complex128:
                    self.assertTrue(out1.precision == 'double')
                    #self.assertTrue(out1.kind == 'complex')                
                
                self.assertTrue(type(out2._scheme) == self.scheme)
                self.assertTrue(type(out2._data) is SchemeArray)
                self.assertEqual(out2[0],5)
                self.assertEqual(out2[1],3)
                self.assertEqual(out2[2],1)
                self.assertTrue(out2.dtype==self.dtype)
                                
                in1-=1
                in2-=1
            # Giving complex input and specifying a real dtype should raise an error
            else:
                self.assertRaises(TypeError, Array, in1, dtype = self.dtype)
                self.assertRaises(TypeError, Array, in2, dtype = self.dtype)
            
            # Also, when it is unspecified
            out3 = Array(in1)
            in1 += 1
            
            self.assertTrue(type(out3._scheme) == self.scheme)
            self.assertTrue(type(out3._data) is SchemeArray)
            self.assertEqual(out3[0],5)
            self.assertEqual(out3[1],3)
            self.assertEqual(out3[2],1)
            self.assertTrue(out3.dtype==self.odtype)
                        
            # We should also be able to create from a CPU Array
            out4 = Array(cpuarray, dtype=self.dtype)
            
            self.assertTrue(type(out4._scheme) == self.scheme)
            self.assertTrue(type(out4._data) is SchemeArray)
            self.assertEqual(out4[0],1)
            self.assertEqual(out4[1],2)
            self.assertEqual(out4[2],3)
            self.assertTrue(out4.dtype==self.dtype)
                        
            self.assertRaises(TypeError, Array,in1, dtype=numpy.int32)
            
            # Check for copy=false
            in3 = Array([5,3,1],dtype=self.dtype)
            out5 = Array(in3,copy=False)
            in3 += 1
            
            self.assertTrue(type(out5._scheme) == self.scheme)
            self.assertTrue(type(out5._data) is SchemeArray)
            self.assertEqual(out5[0],6)
            self.assertEqual(out5[1],4)
            self.assertEqual(out5[2],2)
            self.assertTrue(out5.dtype==self.dtype)
                        
            if _options['scheme'] != 'cpu':
                self.assertRaises(TypeError,Array,cpuarray,copy=False)
                
            # Just checking that we can make an empty array correctly
            empty = Array(numpy.array([]))
            out7 = Array(empty)
            self.assertTrue(out7.dtype==numpy.float64)
            self.assertRaises(IndexError, out7.__getitem__,0)
            
        # Also checking that a cpu array can't be made out of another scheme without copying
        if _options['scheme'] != 'cpu':
            self.assertRaises(TypeError, Array, out4, copy=False)
            out6 = Array(out4, dtype=self.dtype)
            self.assertTrue(type(out6._scheme) == type(None))
            self.assertTrue(type(out6._data) is numpy.ndarray)
            self.assertEqual(out6[0],1)
            self.assertEqual(out6[1],2)
            self.assertEqual(out6[2],3)
            self.assertTrue(out6.dtype==self.dtype)
            
            
    def test_list_init(self):
        with self.context:
            # When specified
            out1 = Array([5,3,1], dtype=self.dtype)
            
            self.assertTrue(type(out1._scheme) == self.scheme)
            self.assertTrue(type(out1._data) is SchemeArray)
            self.assertEqual(out1[0],5)
            self.assertEqual(out1[1],3)
            self.assertEqual(out1[2],1)
            self.assertTrue(out1.dtype==self.dtype)
                        
            if out1.dtype == numpy.float32:
                self.assertTrue(out1.precision == 'single')
                #self.assertTrue(out1.kind == 'real')
            if out1.dtype == numpy.float64:
                self.assertTrue(out1.precision == 'double')
                #self.assertTrue(out1.kind == 'real')
            if out1.dtype == numpy.complex64:
                self.assertTrue(out1.precision == 'single')
                #self.assertTrue(out1.kind == 'complex')
            if out1.dtype == numpy.complex128:
                self.assertTrue(out1.precision == 'double')
                #self.assertTrue(out1.kind == 'complex')  
            
            if self.kind == 'complex':
                out2 = Array([5+0j,3+0j,1+0j], dtype=self.dtype)
            
                self.assertTrue(type(out2._scheme) == self.scheme)
                self.assertTrue(type(out2._data) is SchemeArray)
                self.assertEqual(out2[0],5)
                self.assertEqual(out2[1],3)
                self.assertEqual(out2[2],1)
                self.assertTrue(out2.dtype==self.dtype)
                                
            else:
                self.assertRaises(TypeError, Array,[5+0j, 3+0j, 1+0j],dtype=self.dtype)
            self.assertRaises(TypeError, Array,[1,2,3], dtype=numpy.int32)
            
            #Also, when it is unspecified
            out3 = Array([5,3,1])
            
            self.assertTrue(type(out3._scheme) == self.scheme)
            self.assertTrue(type(out3._data) is SchemeArray)
            self.assertEqual(out3[0],5)
            self.assertEqual(out3[1],3)
            self.assertEqual(out3[2],1)
            self.assertTrue(out3.dtype==numpy.float64)
                        
            out4 = Array([5+0j,3+0j,1+0j])
            
            self.assertTrue(type(out4._scheme) == self.scheme)
            self.assertTrue(type(out4._data) is SchemeArray)
            self.assertEqual(out4[0],5)
            self.assertEqual(out4[1],3)
            self.assertEqual(out4[2],1)
            self.assertTrue(out4.dtype==numpy.complex128)
            
            # Just checking that we can make an empty array correctly
            out7 = Array([])
            self.assertTrue(out7.dtype==numpy.float64)
            self.assertRaises(IndexError, out7.__getitem__,0)
                        
            #We also need to check the zero function
            out5 = zeros(3,dtype=self.dtype)
            out6 = zeros(3)
            
            self.assertTrue(type(out5._scheme) == self.scheme)
            self.assertTrue(type(out5._data) is SchemeArray)
            self.assertEqual(out5[0],0)
            self.assertEqual(out5[1],0)
            self.assertEqual(out5[2],0)
            self.assertTrue(out5.dtype == self.dtype)
                        
            self.assertTrue(type(out6._scheme) == self.scheme)
            self.assertTrue(type(out6._data) is SchemeArray)
            self.assertEqual(out6[0],0)               
            self.assertEqual(out6[1],0)
            self.assertEqual(out6[2],0)
            self.assertTrue(out6.dtype == numpy.float64)
                        
            self.assertRaises(TypeError,Array,[1,2,3],copy=False)
            
def array_test_maker(context,dtype,odtype):
    class tests(ArrayTestBase,unittest.TestCase):
        def __init__(self,*args):
            self.context=context
            self.dtype=dtype
            self.odtype=odtype
            if _options['scheme'] == 'cpu':
                self.scheme = type(None)
            elif _options['scheme'] == 'cuda':
                self.scheme = pycbc.scheme.CUDAScheme
            else:
                self.scheme = pycbc.scheme.OpenCLScheme            
            unittest.TestCase.__init__(self,*args)
    tests.__name__ = _options['scheme'] + " " + dtype.__name__ + " with " + odtype.__name__
    return tests

types = [ (float32,[float32,complex64]), (float64,[float64,complex128]),
        (complex64,[complex64,float32]), (complex128,[float64,complex128]) ]

suite = unittest.TestSuite()

scs =[]
if _options['scheme'] == 'cpu':
    scs.append(DefaultScheme())
if _options['scheme'] == 'cuda':
    scs.append(CUDAScheme(device_num=_options['devicenum']))
if _options['scheme'] == 'opencl':
    scs.append(OpenCLScheme(device_num=_options['devicenum']))

ind = 0
for sc in scs:
    for ty,oktype in types:
        for ot in oktype:
            na = 'test' + str(ind)
            vars()[na] = array_test_maker(sc,ty,ot)
            suite.addTest(unittest.TestLoader().loadTestsFromTestCase(vars()[na]))
            ind += 1



# TODO More specific array tests (instatiation, failure modes, type conversion, etc)
        
        
if __name__ == '__main__':
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
#    elif len(results.errors) == 0:
#        sys.exit(2)
#    elif len(results.errors)==NotImpErrors:
#        sys.exit(3)
#    elif len(results.failures)==0:
#        sys.exit(4)
#    else:
#        sys.exit(5)
        
