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

class tests_base(object):
    def setUp(self):
    
        #We need to check for correct creation from all dtypes, 
        #and errors from incorrect operations so the other precision of 
        #odtype needs to be available as well
        self.other_precision = {numpy.complex64 : numpy.complex128,
                            numpy.complex128 : numpy.complex64,
                            numpy.float32 : numpy.float64,
                            numpy.float64 : numpy.float32}
        
        self.a = Array([5,3,1], dtype = self.dtype)
        self.v = Array([10,8,6],dtype = self.odtype)
        self.bad = Array([1,1,1],dtype = self.other_precision[self.odtype])
        # Number of decimal places to compare for single precision
        if self.a.precision == 'single':
            self.places = 5
        # Number of decimal places to compare for double precision
        if self.a.precision == 'double':
            self.places = 13
        if self.dtype == numpy.float32 or self.dtype == numpy.float64:
            self.kind = 'real'
        else:
            self.kind = 'complex'
            
        if self.kind == 'complex':
            self.comp = Array([5+1j,3+3j,1+5j],dtype=self.dtype)
        pass

    def test_numpy_init(self):
        with self.context:                                
            in1 = numpy.array([5,3,1],dtype=self.odtype)
            in2 = numpy.array([5,3,1],dtype=self.other_precision[self.odtype])
            
            #We don't want to cast complex as real
            if not (self.kind=='real' and (self.odtype == numpy.complex64 or self.odtype==numpy.complex128)):
                #First we must check that the dtype is correct when specified
                out1 = Array(in1, dtype=self.dtype)
                out2 = Array(in2, dtype=self.dtype)
                #to be sure that it is copied
                in1 += 1
                in2 += 1
                self.assertEqual(out1[0],5)
                self.assertEqual(out1[1],3)
                self.assertEqual(out1[2],1)
                self.assertTrue(out1.dtype==self.dtype)
                
                self.assertEqual(out2[0],5)
                self.assertEqual(out2[1],3)
                self.assertEqual(out2[2],1)
                self.assertTrue(out2.dtype==self.dtype)
                in1-=1
                in2-=1
            
            #Also, when it is unspecified
            out3 = Array(in1)
            in1 += 1
            self.assertEqual(out3[0],5)
            self.assertEqual(out3[1],3)
            self.assertEqual(out3[2],1)
            self.assertTrue(out3.dtype==self.odtype)
            
            #Check for copy=false
            if _options['scheme'] == 'cpu':
                in3 = numpy.array([5,3,1],dtype=self.dtype)
                out4 = Array(in3,copy=False)
                in3 += 1
                
                self.assertEqual(out4[0],6)
                self.assertEqual(out4[1],4)
                self.assertEqual(out4[2],2)
                self.assertTrue(out4.dtype==self.dtype)
                    

    def test_array_init(self):
        #this array is made outside the context so we can check that an error is raised when copy = false in a GPU scheme
        cpuarray = Array([1,2,3])
        with self.context:      
            in1 = Array([5,3,1],dtype=self.odtype)
            in2 = Array([5,3,1],dtype=self.other_precision[self.odtype])
            #We don't want to cast complex as real
            if not (self.kind=='real' and (self.odtype == numpy.complex64 or self.odtype==numpy.complex128)):
                #First we must check that the dtype is correct when specified
                out1 = Array(in1, dtype=self.dtype)
                out2 = Array(in2, dtype=self.dtype)
                #to be sure that it is copied
                in1 += 1
                in2 += 1
                
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
                
                self.assertEqual(out2[0],5)
                self.assertEqual(out2[1],3)
                self.assertEqual(out2[2],1)
                self.assertTrue(out2.dtype==self.dtype)
                
                in1-=1
                in2-=1
            
            #Also, when it is unspecified
            out3 = Array(in1)
            in1 += 1
            self.assertEqual(out3[0],5)
            self.assertEqual(out3[1],3)
            self.assertEqual(out3[2],1)
            self.assertTrue(out3.dtype==self.odtype)
            
            #Check for copy=false
            in3 = Array([5,3,1],dtype=self.dtype)
            out4 = Array(in3,copy=False)
            in3 += 1
            
            self.assertEqual(out4[0],6)
            self.assertEqual(out4[1],4)
            self.assertEqual(out4[2],2)
            self.assertTrue(out4.dtype==self.dtype)
            
            if _options['scheme'] != 'cpu':
                self.assertRaises(TypeError,Array,cpuarray,copy=False)
            
    def test_list_init(self):
        with self.context:
            #When specified
            out1 = Array([5,3,1], dtype=self.dtype)
            
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
            
                self.assertEqual(out2[0],5)
                self.assertEqual(out2[1],3)
                self.assertEqual(out2[2],1)
                self.assertTrue(out2.dtype==self.dtype)
            
            #Also, when it is unspecified
            out3 = Array([5,3,1])
            
            self.assertEqual(out3[0],5)
            self.assertEqual(out3[1],3)
            self.assertEqual(out3[2],1)
            self.assertTrue(out3.dtype==numpy.float64)
            
            out4 = Array([5+0j,3+0j,1+0j])
            
            self.assertEqual(out4[0],5)
            self.assertEqual(out4[1],3)
            self.assertEqual(out4[2],1)
            self.assertTrue(out4.dtype==numpy.complex128)
            
            #We also need to check the zero function
            out5 = zeros(3,dtype=self.dtype)
            out6 = zeros(3)
            
            self.assertEqual(out5[0],0)
            self.assertEqual(out5[1],0)
            self.assertEqual(out5[2],0)
            self.assertTrue(out5.dtype == self.dtype)
            
            self.assertEqual(out6[0],0)               
            self.assertEqual(out6[1],0)
            self.assertEqual(out6[2],0)
            self.assertTrue(out6.dtype == numpy.float64)
            
            self.assertRaises(TypeError,Array,[1,2,3],copy=False)
            
    def test_mul(self):
        with self.context:     
            self.assertRaises(TypeError, self.a.__mul__, self.bad)       
            b = self.a * 2
            c = self.a * self.v
            
            self.assertEqual(c[0],50)
            self.assertEqual(c[1],24)
            self.assertEqual(c[2],6)
            
            self.assertEqual(b[0],10)
            self.assertEqual(b[1],6)
            self.assertEqual(b[2],2)
            
            self.assertEqual(self.a[0],5)
            self.assertEqual(self.a[1],3)
            self.assertEqual(self.a[2],1)
            
            self.assertEqual(self.v[0],10)
            self.assertEqual(self.v[1],8)
            self.assertEqual(self.v[2],6)
            
            if self.kind == 'complex':
                d = self.comp * 2
                e = self.comp * self.v
                
                self.assertEqual(e[0],50+10j)
                self.assertEqual(e[1],24+24j)
                self.assertEqual(e[2],6+30j)
                
                self.assertEqual(d[0],10+2j)
                self.assertEqual(d[1],6+6j)
                self.assertEqual(d[2],2+10j)
                
                self.assertEqual(self.comp[0],5+1j)
                self.assertEqual(self.comp[1],3+3j)
                self.assertEqual(self.comp[2],1+5j)
                
                self.assertEqual(self.v[0],10)
                self.assertEqual(self.v[1],8)
                self.assertEqual(self.v[2],6)
                
                
        
    def test_rmul(self):
        with self.context:  
            self.assertRaises(TypeError, self.a.__rmul__, self.bad)       
            b = 2 * self.a 
            c = self.v * self.a
            
            self.assertEqual(c[0],50)
            self.assertEqual(c[1],24)
            self.assertEqual(c[2],6)
            
            self.assertEqual(b[0],10)
            self.assertEqual(b[1],6)
            self.assertEqual(b[2],2)
            
            self.assertEqual(self.a[0],5)
            self.assertEqual(self.a[1],3)
            self.assertEqual(self.a[2],1)
            
            self.assertEqual(self.v[0],10)
            self.assertEqual(self.v[1],8)
            self.assertEqual(self.v[2],6)
            
            if self.kind == 'complex':
                d = 2 * self.comp
                e = self.v * self.comp
                
                self.assertEqual(e[0],50+10j)
                self.assertEqual(e[1],24+24j)
                self.assertEqual(e[2],6+30j)
                
                self.assertEqual(d[0],10+2j)
                self.assertEqual(d[1],6+6j)
                self.assertEqual(d[2],2+10j)
                
                self.assertEqual(self.comp[0],5+1j)
                self.assertEqual(self.comp[1],3+3j)
                self.assertEqual(self.comp[2],1+5j)
                
                self.assertEqual(self.v[0],10)
                self.assertEqual(self.v[1],8)
                self.assertEqual(self.v[2],6)
        
    def test_imul(self):
        with self.context:
            if not (self.kind=='real' and (self.odtype == numpy.complex64 or self.odtype==numpy.complex128)):
                self.assertRaises(TypeError, self.a.__imul__, self.bad)          
                t = self.a * 1 
                t2 = self.a * 1
                t *= 2
                t2 *= self.v
                
                self.assertEqual(t2[0],50)
                self.assertEqual(t2[1],24)
                self.assertEqual(t2[2],6)
                
                self.assertEqual(t[0],10)
                self.assertEqual(t[1],6)
                self.assertEqual(t[2],2)
            
                self.assertEqual(self.a[0],5)
                self.assertEqual(self.a[1],3)
                self.assertEqual(self.a[2],1)
                
                self.assertEqual(self.v[0],10)
                self.assertEqual(self.v[1],8)
                self.assertEqual(self.v[2],6)
            
            if self.kind == 'complex':
                d = self.comp * 1 
                e = self.comp * 1
                d *= 2
                e *= self.v
                
                self.assertEqual(e[0],50+10j)
                self.assertEqual(e[1],24+24j)
                self.assertEqual(e[2],6+30j)
                
                self.assertEqual(d[0],10+2j)
                self.assertEqual(d[1],6+6j)
                self.assertEqual(d[2],2+10j)
                
                self.assertEqual(self.comp[0],5+1j)
                self.assertEqual(self.comp[1],3+3j)
                self.assertEqual(self.comp[2],1+5j)
                
                self.assertEqual(self.v[0],10)
                self.assertEqual(self.v[1],8)
                self.assertEqual(self.v[2],6)
        
    def test_add(self):
        with self.context:
            self.assertRaises(TypeError, self.a.__add__, self.bad)            
            b = self.a + 5
            c = self.a + self.v
            
            self.assertEqual(c[0],15)
            self.assertEqual(c[1],11)
            self.assertEqual(c[2],7)
            
            self.assertEqual(b[0],10)
            self.assertEqual(b[1],8)
            self.assertEqual(b[2],6)
            
            self.assertEqual(self.a[0],5)
            self.assertEqual(self.a[1],3)
            self.assertEqual(self.a[2],1)
            
            self.assertEqual(self.v[0],10)
            self.assertEqual(self.v[1],8)
            self.assertEqual(self.v[2],6)
            
            if self.kind == 'complex':
                d = self.comp + 5
                e = self.comp + self.v
                
                self.assertEqual(e[0],15+1j)
                self.assertEqual(e[1],11+3j)
                self.assertEqual(e[2],7+5j)
                
                self.assertEqual(d[0],10+1j)
                self.assertEqual(d[1],8+3j)
                self.assertEqual(d[2],6+5j)
                
                self.assertEqual(self.comp[0],5+1j)
                self.assertEqual(self.comp[1],3+3j)
                self.assertEqual(self.comp[2],1+5j)
                
                self.assertEqual(self.v[0],10)
                self.assertEqual(self.v[1],8)
                self.assertEqual(self.v[2],6)
        
    def test_radd(self):
        with self.context:    
            self.assertRaises(TypeError, self.a.__radd__, self.bad)            
            b = 5 + self.a
            c = self.v + self.a
            
            self.assertEqual(c[0],15)
            self.assertEqual(c[1],11)
            self.assertEqual(c[2],7)
            
            self.assertEqual(b[0],10)
            self.assertEqual(b[1],8)
            self.assertEqual(b[2],6)
            
            self.assertEqual(self.a[0],5)
            self.assertEqual(self.a[1],3)
            self.assertEqual(self.a[2],1)
            
            self.assertEqual(self.v[0],10)
            self.assertEqual(self.v[1],8)
            self.assertEqual(self.v[2],6)
            
            if self.kind == 'complex':
                d = 5 + self.comp
                e = self.v + self.comp
                
                self.assertEqual(e[0],15+1j)
                self.assertEqual(e[1],11+3j)
                self.assertEqual(e[2],7+5j)
                
                self.assertEqual(d[0],10+1j)
                self.assertEqual(d[1],8+3j)
                self.assertEqual(d[2],6+5j)
                
                self.assertEqual(self.comp[0],5+1j)
                self.assertEqual(self.comp[1],3+3j)
                self.assertEqual(self.comp[2],1+5j)
                
                self.assertEqual(self.v[0],10)
                self.assertEqual(self.v[1],8)
                self.assertEqual(self.v[2],6)
        
    def test_iadd(self):
        with self.context:
            if not (self.kind=='real' and (self.odtype == numpy.complex64 or self.odtype==numpy.complex128)):
                self.assertRaises(TypeError, self.a.__iadd__, self.bad)         
                t = self.a * 1
                t2 = self.a * 1         
                t += 5
                t2 += self.v
                
                self.assertEqual(t2[0],15)
                self.assertEqual(t2[1],11)
                self.assertEqual(t2[2],7)
                
                self.assertEqual(t[0],10)
                self.assertEqual(t[1],8)
                self.assertEqual(t[2],6)
            
                self.assertEqual(self.a[0],5)
                self.assertEqual(self.a[1],3)
                self.assertEqual(self.a[2],1)
                
                self.assertEqual(self.v[0],10)
                self.assertEqual(self.v[1],8)
                self.assertEqual(self.v[2],6)
            
            if self.kind == 'complex':
                d = self.comp * 1 
                e = self.comp * 1
                d += 5
                e += self.v
                
                self.assertEqual(e[0],15+1j)
                self.assertEqual(e[1],11+3j)
                self.assertEqual(e[2],7+5j)
                
                self.assertEqual(d[0],10+1j)
                self.assertEqual(d[1],8+3j)
                self.assertEqual(d[2],6+5j)
                
                self.assertEqual(self.comp[0],5+1j)
                self.assertEqual(self.comp[1],3+3j)
                self.assertEqual(self.comp[2],1+5j)
                
                self.assertEqual(self.v[0],10)
                self.assertEqual(self.v[1],8)
                self.assertEqual(self.v[2],6)

    def test_div(self):
        with self.context:  
            self.assertRaises(TypeError, self.a.__div__, self.bad)               
            b = self.a / 5
            c = self.a /self.v
            
            self.assertAlmostEqual(c[0],1.0/2.0, places=self.places)
            self.assertAlmostEqual(c[1],3.0/8.0, places=self.places)
            self.assertAlmostEqual(c[2],1.0/6.0, places=self.places)
            
            self.assertAlmostEqual(b[0],1.0, places=self.places)
            self.assertAlmostEqual(b[1],3.0/5.0, places=self.places)
            self.assertAlmostEqual(b[2],1.0/5.0, places=self.places)
            
            self.assertEqual(self.a[0],5)
            self.assertEqual(self.a[1],3)
            self.assertEqual(self.a[2],1)
            
            self.assertEqual(self.v[0],10)
            self.assertEqual(self.v[1],8)
            self.assertEqual(self.v[2],6)
            
            if self.kind == 'complex':
                d = self.comp / 5
                e = self.comp / self.v
                
                self.assertAlmostEqual(e[0],1.0/2.0+1.0j/10.0, places=self.places)
                self.assertAlmostEqual(e[1],3.0/8.0+3.0j/8.0, places=self.places)
                self.assertAlmostEqual(e[2],1.0/6.0+5.0j/6.0, places=self.places)
                
                self.assertAlmostEqual(d[0],1.0+1.0j/5.0, places=self.places)
                self.assertAlmostEqual(d[1],3.0/5.0+3.0j/5.0, places=self.places)
                self.assertAlmostEqual(d[2],1.0/5.0+1.0j, places=self.places)
                
                self.assertEqual(self.comp[0],5+1j)
                self.assertEqual(self.comp[1],3+3j)
                self.assertEqual(self.comp[2],1+5j)
                
                self.assertEqual(self.v[0],10)
                self.assertEqual(self.v[1],8)
                self.assertEqual(self.v[2],6)
        
    def test_rdiv(self):
        with self.context:
            self.assertRaises(TypeError, self.a.__rdiv__, self.bad)
            b = 5 / self.a
            c = self.v /self.a
            
            self.assertAlmostEqual(c[0],2.0, places=self.places)
            self.assertAlmostEqual(c[1],8.0/3.0, places=self.places)
            self.assertAlmostEqual(c[2],6.0, places=self.places)
            
            self.assertAlmostEqual(b[0],1.0, places=self.places)
            self.assertAlmostEqual(b[1],5.0/3.0, places=self.places)
            self.assertAlmostEqual(b[2],5.0, places=self.places)
            
            self.assertEqual(self.a[0],5)
            self.assertEqual(self.a[1],3)
            self.assertEqual(self.a[2],1)
            
            self.assertEqual(self.v[0],10)
            self.assertEqual(self.v[1],8)
            self.assertEqual(self.v[2],6)
            
            if self.kind == 'complex':
                d = 5 / self.comp
                e = self.v / self.comp
                
                self.assertAlmostEqual(e[0],25.0/13.0-5.0j/13.0, places=self.places)
                self.assertAlmostEqual(e[1],4.0/3.0-4.0j/3.0, places=self.places)
                self.assertAlmostEqual(e[2],3.0/13.0-15.0j/13.0, places=self.places)
                
                self.assertAlmostEqual(d[0],25.0/26.0-5.0j/26.0, places=self.places)
                self.assertAlmostEqual(d[1],5.0/6.0-5.0j/6.0, places=self.places)
                self.assertAlmostEqual(d[2],5.0/26.0-25.0j/26.0, places=self.places)
                
                self.assertEqual(self.comp[0],5+1j)
                self.assertEqual(self.comp[1],3+3j)
                self.assertEqual(self.comp[2],1+5j)
                
                self.assertEqual(self.v[0],10)
                self.assertEqual(self.v[1],8)
                self.assertEqual(self.v[2],6)
        
    def test_idiv(self):
        with self.context:  
            if not (self.kind=='real' and (self.odtype == numpy.complex64 or self.odtype==numpy.complex128)):
                self.assertRaises(TypeError, self.a.__idiv__, self.bad)        
                t = self.a * 1
                t2 = self.a * 1       
                t /= 5
                t2 /= self.v
                
                self.assertAlmostEqual(t2[0],1.0/2.0, places=self.places)
                self.assertAlmostEqual(t2[1],3.0/8.0, places=self.places)
                self.assertAlmostEqual(t2[2],1.0/6.0, places=self.places)
                
                self.assertAlmostEqual(t[0],1.0, places=self.places)
                self.assertAlmostEqual(t[1],3.0/5.0, places=self.places)
                self.assertAlmostEqual(t[2],1.0/5.0, places=self.places)
            
                self.assertEqual(self.a[0],5)
                self.assertEqual(self.a[1],3)
                self.assertEqual(self.a[2],1)
                
                self.assertEqual(self.v[0],10)
                self.assertEqual(self.v[1],8)
                self.assertEqual(self.v[2],6)
            
            if self.kind == 'complex':
                d = self.comp * 1 
                e = self.comp * 1
                d /= 5
                e /= self.v
                
                self.assertAlmostEqual(e[0],1.0/2.0+1.0j/10.0, places=self.places)
                self.assertAlmostEqual(e[1],3.0/8.0+3.0j/8.0, places=self.places)
                self.assertAlmostEqual(e[2],1.0/6.0+5.0j/6.0, places=self.places)
                
                self.assertAlmostEqual(d[0],1.0+1.0j/5.0, places=self.places)
                self.assertAlmostEqual(d[1],3.0/5.0+3.0j/5.0, places=self.places)
                self.assertAlmostEqual(d[2],1.0/5.0+1.0j, places=self.places)
                
                self.assertEqual(self.comp[0],5+1j)
                self.assertEqual(self.comp[1],3+3j)
                self.assertEqual(self.comp[2],1+5j)
                
                self.assertEqual(self.v[0],10)
                self.assertEqual(self.v[1],8)
                self.assertEqual(self.v[2],6)
            
    def test_sub(self):
        with self.context: 
            self.assertRaises(TypeError, self.a.__sub__, self.bad)              
            b = self.a - 5
            c = self.a - self.v
            
            self.assertEqual(c[0],-5)
            self.assertEqual(c[1],-5)
            self.assertEqual(c[2],-5)
            
            self.assertEqual(b[0],0)
            self.assertEqual(b[1],-2)
            self.assertEqual(b[2],-4)
            
            self.assertEqual(self.a[0],5)
            self.assertEqual(self.a[1],3)
            self.assertEqual(self.a[2],1)
            
            self.assertEqual(self.v[0],10)
            self.assertEqual(self.v[1],8)
            self.assertEqual(self.v[2],6)
            
            if self.kind == 'complex':
                d = self.comp - 5
                e = self.comp - self.v
                
                self.assertEqual(e[0],-5+1j)
                self.assertEqual(e[1],-5+3j)
                self.assertEqual(e[2],-5+5j)
                
                self.assertEqual(d[0],0+1j)
                self.assertEqual(d[1],-2+3j)
                self.assertEqual(d[2],-4+5j)
                
                self.assertEqual(self.comp[0],5+1j)
                self.assertEqual(self.comp[1],3+3j)
                self.assertEqual(self.comp[2],1+5j)
                
                self.assertEqual(self.v[0],10)
                self.assertEqual(self.v[1],8)
                self.assertEqual(self.v[2],6)
        
    def test_rsub(self):
        with self.context:        
            self.assertRaises(TypeError, self.a.__rsub__, self.bad)     
            b = 5 - self.a
            c = self.v - self.a
            
            self.assertEqual(c[0],5)
            self.assertEqual(c[1],5)
            self.assertEqual(c[2],5)
            
            self.assertEqual(b[0],0)
            self.assertEqual(b[1],2)
            self.assertEqual(b[2],4)
            
            self.assertEqual(self.a[0],5)
            self.assertEqual(self.a[1],3)
            self.assertEqual(self.a[2],1)
            
            self.assertEqual(self.v[0],10)
            self.assertEqual(self.v[1],8)
            self.assertEqual(self.v[2],6)
            
            if self.kind == 'complex':
                d = 5 - self.comp
                e = self.v - self.comp
                
                self.assertEqual(e[0],5-1j)
                self.assertEqual(e[1],5-3j)
                self.assertEqual(e[2],5-5j)
                
                self.assertEqual(d[0],0-1j)
                self.assertEqual(d[1],2-3j)
                self.assertEqual(d[2],4-5j)
                
                self.assertEqual(self.comp[0],5+1j)
                self.assertEqual(self.comp[1],3+3j)
                self.assertEqual(self.comp[2],1+5j)
                
                self.assertEqual(self.v[0],10)
                self.assertEqual(self.v[1],8)
                self.assertEqual(self.v[2],6)
        
    def test_isub(self):
        with self.context:
            if not (self.kind=='real' and (self.odtype == numpy.complex64 or self.odtype==numpy.complex128)):
                self.assertRaises(TypeError, self.a.__isub__, self.bad)    
                t = self.a * 1 
                t2 = self.a * 1        
                t -= 5
                t2 -= self.v
                
                self.assertEqual(t2[0],-5)
                self.assertEqual(t2[1],-5)
                self.assertEqual(t2[2],-5)
                
                self.assertEqual(t[0],0)
                self.assertEqual(t[1],-2)
                self.assertEqual(t[2],-4)
            
                self.assertEqual(self.a[0],5)
                self.assertEqual(self.a[1],3)
                self.assertEqual(self.a[2],1)
                
                self.assertEqual(self.v[0],10)
                self.assertEqual(self.v[1],8)
                self.assertEqual(self.v[2],6)      
                
            if self.kind == 'complex':
                d = self.comp * 1 
                e = self.comp * 1
                d -= 5
                e -= self.v
                
                self.assertEqual(e[0],-5+1j)
                self.assertEqual(e[1],-5+3j)
                self.assertEqual(e[2],-5+5j)
                
                self.assertEqual(d[0],0+1j)
                self.assertEqual(d[1],-2+3j)
                self.assertEqual(d[2],-4+5j)
                
                self.assertEqual(self.comp[0],5+1j)
                self.assertEqual(self.comp[1],3+3j)
                self.assertEqual(self.comp[2],1+5j)
                
                self.assertEqual(self.v[0],10)
                self.assertEqual(self.v[1],8)
                self.assertEqual(self.v[2],6)
        
    def test_pow(self):
        with self.context:
            b = self.a ** 2
           # c = 2 ** self.v
           # self.assertTrue(c[0],1024)
            self.assertAlmostEqual(b[0],25, places=self.places)
            self.assertAlmostEqual(b[1],9, places=self.places)
            self.assertAlmostEqual(b[2],1, places=self.places)
            
            self.assertEqual(self.a[0],5)
            self.assertEqual(self.a[1],3)
            self.assertEqual(self.a[2],1)
            
            if self.kind == 'complex':
                c = self.comp ** 2 
                
                self.assertAlmostEqual(c[0],24+10j, places=self.places)
                self.assertAlmostEqual(c[1],18j, places=self.places)
                self.assertAlmostEqual(c[2],-24+10j, places=self.places)
                
                self.assertEqual(self.comp[0],5+1j)
                self.assertEqual(self.comp[1],3+3j)
                self.assertEqual(self.comp[2],1+5j)
                
                self.assertEqual(self.v[0],10)
                self.assertEqual(self.v[1],8)
                self.assertEqual(self.v[2],6)
        
    def test_abs(self):
        with self.context:
            b = abs(self.a)
            c = abs(self.a * -1)
            
            self.assertEqual(b[0],5)
            self.assertEqual(b[1],3)
            self.assertEqual(b[2],1)
            
            self.assertEqual(c[0],5)
            self.assertEqual(c[1],3)
            self.assertEqual(c[2],1)
            
            self.assertEqual(self.a[0],5)
            self.assertEqual(self.a[1],3)
            self.assertEqual(self.a[2],1)
            
            if self.kind == 'complex':
                d = abs(self.comp)
                e = abs(self.comp * -1)

                self.assertAlmostEqual(d[0],pow(26,.5), places=self.places)
                self.assertAlmostEqual(d[1],3*pow(2,.5), places=self.places)
                self.assertAlmostEqual(d[2],pow(26,.5), places=self.places)
                
                self.assertAlmostEqual(e[0],pow(26,.5), places=self.places)
                self.assertAlmostEqual(e[1],3*pow(2,.5), places=self.places)
                self.assertAlmostEqual(e[2],pow(26,.5), places=self.places)
                
                self.assertEqual(self.comp[0],5+1j)
                self.assertEqual(self.comp[1],3+3j)
                self.assertEqual(self.comp[2],1+5j)
                
                self.assertEqual(self.v[0],10)
                self.assertEqual(self.v[1],8)
                self.assertEqual(self.v[2],6)
        
    def test_real(self):
        with self.context:        
            b = self.a.real()
            
            self.assertEqual(b[0],5)
            self.assertEqual(b[1],3)
            self.assertEqual(b[2],1)
            
            self.assertEqual(self.a[0],5)
            self.assertEqual(self.a[1],3)
            self.assertEqual(self.a[2],1)
            
            if self.kind == 'complex':
                c = self.comp.real() 
                
                self.assertEqual(c[0],5)
                self.assertEqual(c[1],3)
                self.assertEqual(c[2],1)
                
                self.assertEqual(self.comp[0],5+1j)
                self.assertEqual(self.comp[1],3+3j)
                self.assertEqual(self.comp[2],1+5j)
                
                self.assertEqual(self.v[0],10)
                self.assertEqual(self.v[1],8)
                self.assertEqual(self.v[2],6)
        
    def test_imag(self):
        with self.context:        
            b = self.a.imag()
            
            self.assertEqual(b[0],0)
            self.assertEqual(b[1],0)
            self.assertEqual(b[2],0)
            
            self.assertEqual(self.a[0],5)
            self.assertEqual(self.a[1],3)
            self.assertEqual(self.a[2],1)
            
            if self.kind == 'complex':
                c = self.comp.imag() 
                
                self.assertEqual(c[0],1)
                self.assertEqual(c[1],3)
                self.assertEqual(c[2],5)
                
                self.assertEqual(self.comp[0],5+1j)
                self.assertEqual(self.comp[1],3+3j)
                self.assertEqual(self.comp[2],1+5j)
                
                self.assertEqual(self.v[0],10)
                self.assertEqual(self.v[1],8)
                self.assertEqual(self.v[2],6)
        
    def test_conj(self):
        with self.context:       
            b = self.a.conj()
            
            self.assertEqual(b[0],5)
            self.assertEqual(b[1],3)
            self.assertEqual(b[2],1)
            
            self.assertEqual(self.a[0],5)
            self.assertEqual(self.a[1],3)
            self.assertEqual(self.a[2],1)
            
            if self.kind == 'complex':
                c = self.comp.conj() 
                
                self.assertEqual(c[0],5-1j)
                self.assertEqual(c[1],3-3j)
                self.assertEqual(c[2],1-5j)
                
                self.assertEqual(self.comp[0],5+1j)
                self.assertEqual(self.comp[1],3+3j)
                self.assertEqual(self.comp[2],1+5j)
                
                self.assertEqual(self.v[0],10)
                self.assertEqual(self.v[1],8)
                self.assertEqual(self.v[2],6)
            
    def test_sum(self):
        with self.context:         
            b = self.a.sum()
            
            self.assertEqual(b,(9))
            
            self.assertEqual(self.a[0],5)
            self.assertEqual(self.a[1],3)
            self.assertEqual(self.a[2],1)
            
            if self.kind == 'complex':
                c = self.comp.sum()
                
                self.assertEqual(c,(9+9j))
                self.assertEqual(self.comp[0],5+1j)
                self.assertEqual(self.comp[1],3+3j)
                self.assertEqual(self.comp[2],1+5j)
                
            
    def test_dot(self):
        with self.context:        
            b = self.a.dot(self.a)
            c = self.a.dot(self.v)
            
            self.assertEqual(b,(35))
            self.assertEqual(c,(80))
            
            self.assertEqual(self.v[0],10)
            self.assertEqual(self.v[1],8)
            self.assertEqual(self.v[2],6)
            
            self.assertEqual(self.a[0],5)
            self.assertEqual(self.a[1],3)
            self.assertEqual(self.a[2],1)
            
            if self.kind == 'complex':
                d = self.comp.dot(self.comp)
                e = self.comp.dot(self.v)
                
                self.assertEqual(d,(38j))
                self.assertEqual(e,(80+64j))
                
                self.assertEqual(self.comp[0],5+1j)
                self.assertEqual(self.comp[1],3+3j)
                self.assertEqual(self.comp[2],1+5j)
                
                self.assertEqual(self.v[0],10)
                self.assertEqual(self.v[1],8)
                self.assertEqual(self.v[2],6)
    
    def test_max(self):
        with self.context:
            #When updated, this should call self.a.kind
            if self.kind == 'real':
                self.assertEqual(self.a.max(),5)
                
                self.assertEqual(self.a[0],5)
                self.assertEqual(self.a[1],3)
                self.assertEqual(self.a[2],1)
                
    def test_min(self):
        with self.context:
            #When updated, this should call self.a.kind
            if self.kind == 'real':
                self.assertEqual(self.a.min(),1)
                
                self.assertEqual(self.a[0],5)
                self.assertEqual(self.a[1],3)
                self.assertEqual(self.a[2],1)
                
    

def array_test_maker(context,dtype,odtype):
    class tests(tests_base,unittest.TestCase):
        def __init__(self,*args):
            self.context=context
            self.dtype=dtype
            self.odtype=odtype
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
    unittest.TextTestRunner(verbosity=2).run(suite)
