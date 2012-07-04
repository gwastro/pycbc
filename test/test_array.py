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
            
        # Here two test arrays are created. Their content depends on their dtype,
        # so that we can make use of non-zero imaginary parts.
        # These arrays are just arbitrarily chosen, but in such a way that they are not
        # oversimplified, e.g. we don't want to multiply by one, or add zero, or have
        # the answer be one or zero.
        # We will make 2 of each, so that we can always check that the operation has not
        # altered anything it shouldn't.
        if self.kind == 'real':
            self.a = Array([5,3,1], dtype=self.dtype)
            self.a2 = Array([5,3,1], dtype=self.dtype)
        else:
            self.a = Array([5+1j,3+3j,1+5j], dtype=self.dtype)
            self.a2 = Array([5+1j,3+3j,1+5j], dtype=self.dtype)
            
        if self.okind == 'real':
            self.b = Array([10,8,6], dtype=self.odtype)
            self.b2 = Array([10,8,6], dtype=self.odtype)
        else:
            self.b = Array([10+6j,8+4j,6+2j], dtype=self.odtype)
            self.b2 = Array([10+6j,8+4j,6+2j], dtype=self.odtype)
        
        # We will need to also test a non-zero imaginary part scalar.
        # For simplicity, we will test with a real scalar when odtype is real, 
        # and a complex scalar when odtype is complex. This will prevent redundancy,
        # and make it easier to spot problems, just based off of the test names.
        if self.okind =='real':
            self.s = 5
            self.s2 = 5
        else:
            self.s = 5+2j
            self.s2 = 5+2j

        # Finally, we want to have an array that we shouldn't be able to operate on,
        # because the precision is wrong.
        self.bad = Array([1,1,1],dtype = self.other_precision[self.odtype])
        
        # All the answers are stored here to make it easier to read in the actual tests.
        # Again, it makes a difference whether they are complex or real valued, so there
        # are four sets of possible answers, depending on the dtypes.
        if self.kind == 'real' and self.okind == 'real':
        
            self.mul = [50, 24, 6]
            self.mul_s = [25, 15, 5]
                        
            self.add = [15, 11, 7]
            self.add_s = [10, 8, 6]
                        
            #self.div = [.5, 3./8., 1./6.]
            self.div = [.5, 0.375, .16666666666666666667]
            #self.div_s = [1., 3./5., 1./5.]
            self.div_s = [1., 0.6, 0.2]
                        
            #self.rdiv = [2., 8./3., 6.]
            self.rdiv = [2., 2.66666666666666666667, 6.]
            #self.rdiv_s = [1., 5./3., 5.]
            self.rdiv_s = [1., 1.66666666666666666667, 5.]
            
            self.sub = [-5, -5, -5]
            self.sub_s = [0, -2, -4]
            
            self.rsub = [5, 5, 5]
            self.rsub_s = [0, 2, 4]
            
            self.pow1 = [25., 9., 1.]
            #self.pow2 = [pow(5,-1.5), pow(3,-1.5), pow(1,-1.5)]
            self.pow2 = [0.08944271909999158786, 0.19245008972987525484, 1.]
            
            self.abs = [5, 3, 1]
            
            self.real = [5,3,1]
            self.imag = [0, 0, 0]
            self.conj = [5, 3, 1]
            
            self.sum = 9
            
            self.dot = 80
                        
        if self.kind =='real' and self.okind == 'complex':
            
            self.mul = [50+30j, 24+12j, 6+2j]
            self.mul_s = [25+10j, 15+6j, 5+2j]
            
            self.add = [15+6j, 11+4j, 7+2j]
            self.add_s = [10+2j, 8+2j, 6+2j]
            
            #self.div = [25./68.-15.j/68., 3./10.-3.j/20., 3./20.-1.j/20.] 
            self.div = [0.36764705882352941176-0.22058823529411764706j, 0.3-0.15j, 0.15-0.05j] 
            #self.div_s = [25./29.-10.j/29., 15./29.-6.j/29., 5./29.-2.j/29.]
            self.div_s = [0.86206896551724137931-0.34482758620689655172j,
                          0.51724137931034482759-0.20689655172413793103j,
                          0.17241379310344827586-0.06896551724137931034j]
            
            #self.rdiv = [2.+6.j/5., 8./3.+4.j/3, 6.+2.j]
            self.rdiv = [2.+1.2j, 2.66666666666666666667+1.33333333333333333333j, 6.+2.j]
            #self.rdiv_s = [1.+2.j/5., 5./3.+2.j/3., 5.+2.j]
            self.rdiv_s = [1.+0.4j, 1.66666666666666666667+0.666666666666666666667j, 5.+2.j]
            
            self.sub = [-5-6j, -5-4j, -5-2j]
            self.sub_s = [0-2j, -2-2j, -4-2j]
            
            self.rsub = [5+6j, 5+4j, 5+2j]
            self.rsub_s = [0+2j, 2+2j, 4+2j]
            
            self.pow1 = [25., 9., 1.]
            #self.pow2 = [pow(5,-1.5), pow(3,-1.5), pow(1,-1.5)]
            self.pow2 = [0.08944271909999158786, 0.19245008972987525484, 1.]
            
            self.abs = [5, 3, 1]
            
            self.real = [5,3,1]
            self.imag = [0, 0, 0]
            self.conj = [5, 3, 1]
            
            self.sum = 9
            
            self.dot = 80+44j
            
        if self.kind == 'complex' and self.okind == 'real':
            
            self.mul = [50+10j, 24+24j, 6+30j]
            self.mul_s = [25+5j, 15+15j, 5+25j]
            
            self.add = [15+1j, 11+3j, 7+5j]
            self.add_s = [10+1j, 8+3j, 6+5j]
            
            #self.div = [1./2.+1.j/10., 3./8.+3.j/8., 1./6.+5.j/6.]
            self.div = [0.5+0.1j, 0.375+0.375j, 0.16666666666666666667+0.83333333333333333333j]
            #self.div_s = [1.+1.j/5., 3./5.+3.j/5., 1./5.+1.j]
            self.div_s = [1.+0.2j, 0.6+0.6j, 0.2+1.j]
            
            #self.rdiv = [25./13.-5.j/13., 4./3.-4.j/3., 3./13.-15.j/13.]
            self.rdiv = [1.92307692307692307692-0.38461538461538461538j,
                         1.33333333333333333333-1.33333333333333333333j,
                         0.23076923076923076923-1.15384615384615384615j]
            #self.rdiv_s = [25./26.-5.j/26., 5./6.-5.j/6., 5./26.-25.j/26.]
            self.rdiv_s = [0.96153846153846153846-0.19230769230769230769j,
                           0.83333333333333333333-0.83333333333333333333j,
                           0.19230769230769230769-0.96153846153846153846j]
            
            self.sub = [-5+1j, -5+3j, -5+5j]
            self.sub_s = [0+1j, -2+3j, -4+5j]
            
            self.rsub = [5-1j, 5-3j, 5-5j]
            self.rsub_s = [0-1j, 2-3j, 4-5j]
            
            self.pow1 = [24.+10.j, 0.+18.j, -24.+10.j]
            #self.pow2 = [pow(5+1j,-1.5), pow(3+3j,-1.5), pow(1+5j,-1.5)]
            self.pow2 = [0.08307064054041229214-0.0253416052125975132j,
                         0.04379104225017853491-0.1057209281108342370j,
                        -0.04082059235165559671-0.0766590341356157206j]
            
            #self.abs = [pow(26,.5), 3*pow(2,.5), pow(26,.5)]
            self.abs = [5.09901951359278483003,
                        4.24264068711928514641,
                        5.09901951359278483003]
            
            self.real = [5,3,1]
            self.imag = [1, 3, 5]
            self.conj = [5-1j, 3-3j, 1-5j]
            
            self.sum = 9+9j
            
            self.dot = 80+64j
            
        if self.kind =='complex' and self.okind =='complex':
            
            self.mul = [44+40j, 12+36j, -4+32j]
            self.mul_s = [23+15j, 9+21j, -5+27j]
            
            self.add = [15+7j, 11+7j, 7+7j]
            self.add_s = [10+3j, 8+5j, 6+7j]
            
            #self.div = [7./17.-5.j/34., 9./20.+3.j/20., 2./5.+7.j/10.]
            self.div = [0.41176470588235294118-0.14705882352941176471j, 0.45+0.15j, 0.4+0.7j]
            #self.div_s = [27./29.-5.j/29., 21./29.+9.j/29., 15./29.+23.j/29.]
            self.div_s = [0.93103448275862068966-0.17241379310344827586j,
                          0.72413793103448275862+0.31034482758620689655j,
                          0.51724137931034482759+0.79310344827586206897j]
            
            #self.rdiv = [28./13.+10.j/13., 2.-2.j/3., 8./13.-14.j/13.]
            self.rdiv = [2.15384615384615384615+0.76923076923076923077j,
                         2.                    -0.66666666666666666667j,
                         0.61538461538461538462-1.07692307692307692308j]
            #self.rdiv_s = [27./26.+5.j/26., 7./6.-1.j/2., 15./26.-23.j/26]             
            self.rdiv_s = [1.03846153846153846154+0.19230769230769230769j,
                           1.16666666666666666667-0.5j,
                           0.57692307692307692308-0.88461538461538461538j]
            
            self.sub = [-5-5j, -5-1j, -5+3j]
            self.sub_s = [0-1j, -2+1j, -4+3j]
            
            self.rsub = [5+5j, 5+1j, 5-3j]
            self.rsub_s = [0+1j, 2-1j, 4-3j]
            
            self.pow1 = [24.+10.j, 0.+18.j, -24.+10.j]
            #self.pow2 = [pow(5+1j,-1.5), pow(3+3j,-1.5), pow(1+5j,-1.5)]
            self.pow2 = [0.08307064054041229214-0.0253416052125975132j,
                         0.04379104225017853491-0.1057209281108342370j,
                        -0.04082059235165559671-0.0766590341356157206j]
            
            #self.abs = [pow(26,.5), 3*pow(2,.5), pow(26,.5)]
            self.abs = [5.09901951359278483003,
                        4.24264068711928514641,
                        5.09901951359278483003]
            
            self.real = [5,3,1]
            self.imag = [1, 3, 5]
            self.conj = [5-1j, 3-3j, 1-5j]
            
            self.sum = 9+9j
            
            self.dot = 52+108j
            

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
                self.assertEqual(out1[0],5)
                self.assertEqual(out1[1],3)
                self.assertEqual(out1[2],1)
                self.assertTrue(out1.dtype==self.dtype)
                self.assertTrue(out1._scheme == self.scheme)
                
                self.assertEqual(out2[0],5)
                self.assertEqual(out2[1],3)
                self.assertEqual(out2[2],1)
                self.assertTrue(out2.dtype==self.dtype)
                self.assertTrue(out2._scheme == self.scheme)
                in1-=1
                in2-=1
            
            #Also, when it is unspecified
            out3 = Array(in1)
            in1 += 1
            self.assertEqual(out3[0],5)
            self.assertEqual(out3[1],3)
            self.assertEqual(out3[2],1)
            self.assertTrue(out3.dtype==self.odtype)
            self.assertTrue(out3._scheme == self.scheme)
            
            #Check for copy=false
            if _options['scheme'] == 'cpu':
                in3 = numpy.array([5,3,1],dtype=self.dtype)
                out4 = Array(in3,copy=False)
                in3 += 1
                
                self.assertEqual(out4[0],6)
                self.assertEqual(out4[1],4)
                self.assertEqual(out4[2],2)
                self.assertTrue(out4.dtype==self.dtype)
                self.assertTrue(out4._scheme == self.scheme)
                    

    def test_array_init(self):
        #this array is made outside the context so we can check that an error is raised when copy = false in a GPU scheme
        cpuarray = Array([1,2,3])
        with self.context:      
            in1 = Array([5,3,1],dtype=self.odtype)
            in2 = Array([5,3,1],dtype=self.other_precision[self.odtype])
            self.assertTrue(in1._scheme == self.scheme)
            self.assertTrue(in2._scheme == self.scheme)
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
                self.assertTrue(out1._scheme == self.scheme)
                
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
                self.assertTrue(out2._scheme == self.scheme)
                
                in1-=1
                in2-=1
            
            #Also, when it is unspecified
            out3 = Array(in1)
            in1 += 1
            self.assertEqual(out3[0],5)
            self.assertEqual(out3[1],3)
            self.assertEqual(out3[2],1)
            self.assertTrue(out3.dtype==self.odtype)
            self.assertTrue(out3._scheme == self.scheme)
            
            #Check for copy=false
            in3 = Array([5,3,1],dtype=self.dtype)
            out4 = Array(in3,copy=False)
            in3 += 1
            
            self.assertEqual(out4[0],6)
            self.assertEqual(out4[1],4)
            self.assertEqual(out4[2],2)
            self.assertTrue(out4.dtype==self.dtype)
            self.assertTrue(out4._scheme == self.scheme)
            
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
            self.assertTrue(out1._scheme == self.scheme)
            
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
                self.assertTrue(out2._scheme == self.scheme)
            
            #Also, when it is unspecified
            out3 = Array([5,3,1])
            
            self.assertEqual(out3[0],5)
            self.assertEqual(out3[1],3)
            self.assertEqual(out3[2],1)
            self.assertTrue(out3.dtype==numpy.float64)
            self.assertTrue(out3._scheme == self.scheme)
            
            out4 = Array([5+0j,3+0j,1+0j])
            
            self.assertEqual(out4[0],5)
            self.assertEqual(out4[1],3)
            self.assertEqual(out4[2],1)
            self.assertTrue(out4.dtype==numpy.complex128)
            self.assertTrue(out4._scheme == self.scheme)
            
            #We also need to check the zero function
            out5 = zeros(3,dtype=self.dtype)
            out6 = zeros(3)
            
            self.assertEqual(out5[0],0)
            self.assertEqual(out5[1],0)
            self.assertEqual(out5[2],0)
            self.assertTrue(out5.dtype == self.dtype)
            self.assertTrue(out5._scheme == self.scheme)
            
            self.assertEqual(out6[0],0)               
            self.assertEqual(out6[1],0)
            self.assertEqual(out6[2],0)
            self.assertTrue(out6.dtype == numpy.float64)
            self.assertTrue(out6._scheme == self.scheme)
            
            self.assertRaises(TypeError,Array,[1,2,3],copy=False)
            
    def test_set(self):
        with self.context:
            t1 = self.a * 1
            
            if not (self.kind == 'real' and self.okind == 'complex'):                
                t1[0] = Array(self.b[0])
                t1[2] = Array(self.b[2])
                
                self.assertEqual(t1[0],self.b2[0])
                self.assertEqual(t1[1],self.a2[1])
                self.assertEqual(t1[2],self.b2[2])
                
            else:
                self.assertRaises(TypeError, t1.__setitem__, 0, Array(self.b[0]))
                
            self.assertEqual(self.a[0],self.a2[0])
            self.assertEqual(self.a[1],self.a2[1])
            self.assertEqual(self.a[2],self.a2[2])
            
            self.assertEqual(self.b[0],self.b2[0])
            self.assertEqual(self.b[1],self.b2[1])
            self.assertEqual(self.b[2],self.b2[2])
            
    def test_mul(self):
        with self.context:  
           
            self.assertRaises(TypeError, self.a.__mul__, self.bad)  
                 
            t1 = self.a * self.s
            t2 = self.a * self.b
            
            self.assertEqual(t1[0],self.mul_s[0])
            self.assertEqual(t1[1],self.mul_s[1])
            self.assertEqual(t1[2],self.mul_s[2])
            
            self.assertEqual(t2[0],self.mul[0])
            self.assertEqual(t2[1],self.mul[1])
            self.assertEqual(t2[2],self.mul[2])

            self.assertEqual(self.a[0],self.a2[0])
            self.assertEqual(self.a[1],self.a2[1])
            self.assertEqual(self.a[2],self.a2[2])
            
            self.assertEqual(self.b[0],self.b2[0])
            self.assertEqual(self.b[1],self.b2[1])
            self.assertEqual(self.b[2],self.b2[2])
            
            self.assertEqual(self.s,self.s2)
        
    def test_rmul(self):
        with self.context:  
        
            self.assertRaises(TypeError, self.a.__rmul__, self.bad)  
                 
            t1 = self.s * self.a
            t2 = self.a.__rmul__(self.b)
            
            self.assertEqual(t1[0],self.mul_s[0])
            self.assertEqual(t1[1],self.mul_s[1])
            self.assertEqual(t1[2],self.mul_s[2])
            
            self.assertEqual(t2[0],self.mul[0])
            self.assertEqual(t2[1],self.mul[1])
            self.assertEqual(t2[2],self.mul[2])

            self.assertEqual(self.a[0],self.a2[0])
            self.assertEqual(self.a[1],self.a2[1])
            self.assertEqual(self.a[2],self.a2[2])
            
            self.assertEqual(self.b[0],self.b2[0])
            self.assertEqual(self.b[1],self.b2[1])
            self.assertEqual(self.b[2],self.b2[2])
            
            self.assertEqual(self.s,self.s2)
                
    def test_imul(self):
        with self.context:  
            self.assertRaises(TypeError, self.a.__imul__, self.bad)
            
            t1 = self.a * 1
            t2 = self.a * 1
            
            if not (self.kind == 'real' and self.okind == 'complex'):
                t1 *= self.s
                t2 *= self.b
                
                self.assertEqual(t1[0],self.mul_s[0])
                self.assertEqual(t1[1],self.mul_s[1])
                self.assertEqual(t1[2],self.mul_s[2])
                
                self.assertEqual(t2[0],self.mul[0])
                self.assertEqual(t2[1],self.mul[1])
                self.assertEqual(t2[2],self.mul[2])

            else:
                self.assertRaises(TypeError, t1.__imul__,self.s)
                self.assertRaises(TypeError, t1.__imul__,self.b)
                
                self.assertEqual(t1[0],self.a2[0])
                self.assertEqual(t1[1],self.a2[1])
                self.assertEqual(t1[2],self.a2[2])
                
            self.assertEqual(self.a[0],self.a2[0])
            self.assertEqual(self.a[1],self.a2[1])
            self.assertEqual(self.a[2],self.a2[2])
            
            self.assertEqual(self.b[0],self.b2[0])
            self.assertEqual(self.b[1],self.b2[1])
            self.assertEqual(self.b[2],self.b2[2])
            
            self.assertEqual(self.s,self.s2)
        
    def test_add(self):
        with self.context:  
           
            self.assertRaises(TypeError, self.a.__add__, self.bad)  
                 
            t1 = self.a + self.s
            t2 = self.a + self.b
            
            self.assertEqual(t1[0],self.add_s[0])
            self.assertEqual(t1[1],self.add_s[1])
            self.assertEqual(t1[2],self.add_s[2])
            
            self.assertEqual(t2[0],self.add[0])
            self.assertEqual(t2[1],self.add[1])
            self.assertEqual(t2[2],self.add[2])

            self.assertEqual(self.a[0],self.a2[0])
            self.assertEqual(self.a[1],self.a2[1])
            self.assertEqual(self.a[2],self.a2[2])
            
            self.assertEqual(self.b[0],self.b2[0])
            self.assertEqual(self.b[1],self.b2[1])
            self.assertEqual(self.b[2],self.b2[2])
            
            self.assertEqual(self.s,self.s2)
        
    def test_radd(self):
        with self.context:  
           
            self.assertRaises(TypeError, self.a.__radd__, self.bad)  
                 
            t1 = self.s + self.a
            t2 = self.a.__radd__(self.b)
            
            self.assertEqual(t1[0],self.add_s[0])
            self.assertEqual(t1[1],self.add_s[1])
            self.assertEqual(t1[2],self.add_s[2])
            
            self.assertEqual(t2[0],self.add[0])
            self.assertEqual(t2[1],self.add[1])
            self.assertEqual(t2[2],self.add[2])

            self.assertEqual(self.a[0],self.a2[0])
            self.assertEqual(self.a[1],self.a2[1])
            self.assertEqual(self.a[2],self.a2[2])
            
            self.assertEqual(self.b[0],self.b2[0])
            self.assertEqual(self.b[1],self.b2[1])
            self.assertEqual(self.b[2],self.b2[2])
            
            self.assertEqual(self.s,self.s2)
        
    def test_iadd(self):
        with self.context:  
            self.assertRaises(TypeError, self.a.__iadd__, self.bad)
            
            t1 = self.a * 1
            t2 = self.a * 1
            
            if not (self.kind == 'real' and self.okind == 'complex'):
                t1 += self.s
                t2 += self.b
                
                self.assertEqual(t1[0],self.add_s[0])
                self.assertEqual(t1[1],self.add_s[1])
                self.assertEqual(t1[2],self.add_s[2])
                
                self.assertEqual(t2[0],self.add[0])
                self.assertEqual(t2[1],self.add[1])
                self.assertEqual(t2[2],self.add[2])

            else:
                self.assertRaises(TypeError, t1.__iadd__,self.s)
                self.assertRaises(TypeError, t1.__iadd__,self.b)
                
                self.assertEqual(t1[0],self.a2[0])
                self.assertEqual(t1[1],self.a2[1])
                self.assertEqual(t1[2],self.a2[2])
                
            self.assertEqual(self.a[0],self.a2[0])
            self.assertEqual(self.a[1],self.a2[1])
            self.assertEqual(self.a[2],self.a2[2])
            
            self.assertEqual(self.b[0],self.b2[0])
            self.assertEqual(self.b[1],self.b2[1])
            self.assertEqual(self.b[2],self.b2[2])
            
            self.assertEqual(self.s,self.s2)

    def test_div(self):
        with self.context:  
           
            self.assertRaises(TypeError, self.a.__div__, self.bad)  
                 
            t1 = self.a / self.s
            t2 = self.a / self.b
            
            self.assertAlmostEqual(t1[0],self.div_s[0], places=self.places)
            self.assertAlmostEqual(t1[1],self.div_s[1], places=self.places)
            self.assertAlmostEqual(t1[2],self.div_s[2], places=self.places)
            
            self.assertAlmostEqual(t2[0],self.div[0], places=self.places)
            self.assertAlmostEqual(t2[1],self.div[1], places=self.places)
            self.assertAlmostEqual(t2[2],self.div[2], places=self.places)

            self.assertEqual(self.a[0],self.a2[0])
            self.assertEqual(self.a[1],self.a2[1])
            self.assertEqual(self.a[2],self.a2[2])
            
            self.assertEqual(self.b[0],self.b2[0])
            self.assertEqual(self.b[1],self.b2[1])
            self.assertEqual(self.b[2],self.b2[2])
            
            self.assertEqual(self.s,self.s2)
        
    def test_rdiv(self):
        with self.context:  
           
            self.assertRaises(TypeError, self.a.__rdiv__, self.bad)  
                 
            t1 = self.s / self.a
            t2 = self.a.__rdiv__(self.b)
            
            self.assertAlmostEqual(t1[0],self.rdiv_s[0], places=self.places)
            self.assertAlmostEqual(t1[1],self.rdiv_s[1], places=self.places)
            self.assertAlmostEqual(t1[2],self.rdiv_s[2], places=self.places)
            
            self.assertAlmostEqual(t2[0],self.rdiv[0], places=self.places)
            self.assertAlmostEqual(t2[1],self.rdiv[1], places=self.places)
            self.assertAlmostEqual(t2[2],self.rdiv[2], places=self.places)

            self.assertEqual(self.a[0],self.a2[0])
            self.assertEqual(self.a[1],self.a2[1])
            self.assertEqual(self.a[2],self.a2[2])
            
            self.assertEqual(self.b[0],self.b2[0])
            self.assertEqual(self.b[1],self.b2[1])
            self.assertEqual(self.b[2],self.b2[2])
            
            self.assertEqual(self.s,self.s2)
        
    def test_idiv(self):
        with self.context:  
            self.assertRaises(TypeError, self.a.__idiv__, self.bad)
            
            t1 = self.a * 1
            t2 = self.a * 1
            
            if not (self.kind == 'real' and self.okind == 'complex'):
                t1 /= self.s
                t2 /= self.b
                
                self.assertAlmostEqual(t1[0],self.div_s[0], places=self.places)
                self.assertAlmostEqual(t1[1],self.div_s[1], places=self.places)
                self.assertAlmostEqual(t1[2],self.div_s[2], places=self.places)
                
                self.assertAlmostEqual(t2[0],self.div[0], places=self.places)
                self.assertAlmostEqual(t2[1],self.div[1], places=self.places)
                self.assertAlmostEqual(t2[2],self.div[2], places=self.places)

            else:
                self.assertRaises(TypeError, t1.__idiv__,self.s)
                self.assertRaises(TypeError, t1.__idiv__,self.b)
                
                self.assertEqual(t1[0],self.a2[0])
                self.assertEqual(t1[1],self.a2[1])
                self.assertEqual(t1[2],self.a2[2])
                
            self.assertEqual(self.a[0],self.a2[0])
            self.assertEqual(self.a[1],self.a2[1])
            self.assertEqual(self.a[2],self.a2[2])
            
            self.assertEqual(self.b[0],self.b2[0])
            self.assertEqual(self.b[1],self.b2[1])
            self.assertEqual(self.b[2],self.b2[2])
            
            self.assertEqual(self.s,self.s2)
            
    def test_sub(self):
        with self.context:  
           
            self.assertRaises(TypeError, self.a.__sub__, self.bad)  
                 
            t1 = self.a - self.s
            t2 = self.a - self.b
            
            self.assertEqual(t1[0],self.sub_s[0])
            self.assertEqual(t1[1],self.sub_s[1])
            self.assertEqual(t1[2],self.sub_s[2])
            
            self.assertEqual(t2[0],self.sub[0])
            self.assertEqual(t2[1],self.sub[1])
            self.assertEqual(t2[2],self.sub[2])

            self.assertEqual(self.a[0],self.a2[0])
            self.assertEqual(self.a[1],self.a2[1])
            self.assertEqual(self.a[2],self.a2[2])
            
            self.assertEqual(self.b[0],self.b2[0])
            self.assertEqual(self.b[1],self.b2[1])
            self.assertEqual(self.b[2],self.b2[2])
            
            self.assertEqual(self.s,self.s2)
        
    def test_rsub(self):
        with self.context:  
           
            self.assertRaises(TypeError, self.a.__rsub__, self.bad)  
                 
            t1 = self.s - self.a
            t2 = self.a.__rsub__(self.b)
            
            self.assertEqual(t1[0],self.rsub_s[0])
            self.assertEqual(t1[1],self.rsub_s[1])
            self.assertEqual(t1[2],self.rsub_s[2])
            
            self.assertEqual(t2[0],self.rsub[0])
            self.assertEqual(t2[1],self.rsub[1])
            self.assertEqual(t2[2],self.rsub[2])

            self.assertEqual(self.a[0],self.a2[0])
            self.assertEqual(self.a[1],self.a2[1])
            self.assertEqual(self.a[2],self.a2[2])
            
            self.assertEqual(self.b[0],self.b2[0])
            self.assertEqual(self.b[1],self.b2[1])
            self.assertEqual(self.b[2],self.b2[2])
            
            self.assertEqual(self.s,self.s2)
        
    def test_isub(self):
        with self.context:  
            self.assertRaises(TypeError, self.a.__isub__, self.bad)
            
            t1 = self.a * 1
            t2 = self.a * 1
            
            if not (self.kind == 'real' and self.okind == 'complex'):
                t1 -= self.s
                t2 -= self.b
                
                self.assertEqual(t1[0],self.sub_s[0])
                self.assertEqual(t1[1],self.sub_s[1])
                self.assertEqual(t1[2],self.sub_s[2])
                
                self.assertEqual(t2[0],self.sub[0])
                self.assertEqual(t2[1],self.sub[1])
                self.assertEqual(t2[2],self.sub[2])

            else:
                self.assertRaises(TypeError, t1.__isub__,self.s)
                self.assertRaises(TypeError, t1.__isub__,self.b)
                
                self.assertEqual(t1[0],self.a2[0])
                self.assertEqual(t1[1],self.a2[1])
                self.assertEqual(t1[2],self.a2[2])
                
            self.assertEqual(self.a[0],self.a2[0])
            self.assertEqual(self.a[1],self.a2[1])
            self.assertEqual(self.a[2],self.a2[2])
            
            self.assertEqual(self.b[0],self.b2[0])
            self.assertEqual(self.b[1],self.b2[1])
            self.assertEqual(self.b[2],self.b2[2])
            
            self.assertEqual(self.s,self.s2)
        
    def test_pow(self):
        with self.context:
            t1 = self.a ** 2
            t2 = self.a ** -1.5
            
            self.assertAlmostEqual(t1[0],self.pow1[0], places=self.places)
            self.assertAlmostEqual(t1[1],self.pow1[1], places=self.places)
            self.assertAlmostEqual(t1[2],self.pow1[2], places=self.places)
            
            self.assertAlmostEqual(t2[0],self.pow2[0], places=self.places)
            self.assertAlmostEqual(t2[1],self.pow2[1], places=self.places)
            self.assertAlmostEqual(t2[2],self.pow2[2], places=self.places)
            
            self.assertEqual(self.a[0],self.a2[0])
            self.assertEqual(self.a[1],self.a2[1])
            self.assertEqual(self.a[2],self.a2[2])
        
    def test_abs(self):
        with self.context:
            
            # We want to check that absolute value behaves correctly no matter
            # what quadran it's in
            
            t1 = abs(self.a * 1)
            t2 = abs(self.a * -1)
            t3 = abs(self.a * 1j)
            t4 = abs(self.a * -1j)
            
            self.assertAlmostEqual(t1[0],self.abs[0], places=self.places)
            self.assertAlmostEqual(t1[1],self.abs[1], places=self.places)
            self.assertAlmostEqual(t1[2],self.abs[2], places=self.places)
            
            self.assertAlmostEqual(t2[0],self.abs[0], places=self.places)
            self.assertAlmostEqual(t2[1],self.abs[1], places=self.places)
            self.assertAlmostEqual(t2[2],self.abs[2], places=self.places)
            
            self.assertAlmostEqual(t3[0],self.abs[0], places=self.places)
            self.assertAlmostEqual(t3[1],self.abs[1], places=self.places)
            self.assertAlmostEqual(t3[2],self.abs[2], places=self.places)
            
            self.assertAlmostEqual(t4[0],self.abs[0], places=self.places)
            self.assertAlmostEqual(t4[1],self.abs[1], places=self.places)
            self.assertAlmostEqual(t4[2],self.abs[2], places=self.places)
            
            self.assertEqual(self.a[0],self.a2[0])
            self.assertEqual(self.a[1],self.a2[1])
            self.assertEqual(self.a[2],self.a2[2])
        
    def test_real(self):
        with self.context:        
            t1 = self.a.real()
            
            self.assertEqual(t1[0],self.real[0])
            self.assertEqual(t1[1],self.real[1])
            self.assertEqual(t1[2],self.real[2])
            
            self.assertEqual(self.a[0],self.a2[0])
            self.assertEqual(self.a[1],self.a2[1])
            self.assertEqual(self.a[2],self.a2[2])
        
    def test_imag(self):
        with self.context:        
            t1 = self.a.imag()
            
            self.assertEqual(t1[0],self.imag[0])
            self.assertEqual(t1[1],self.imag[1])
            self.assertEqual(t1[2],self.imag[2])
            
            self.assertEqual(self.a[0],self.a2[0])
            self.assertEqual(self.a[1],self.a2[1])
            self.assertEqual(self.a[2],self.a2[2])
        
    def test_conj(self):
        with self.context:       
            t1 = self.a.conj()
            
            self.assertEqual(t1[0],self.conj[0])
            self.assertEqual(t1[1],self.conj[1])
            self.assertEqual(t1[2],self.conj[2])
            
            self.assertEqual(self.a[0],self.a2[0])
            self.assertEqual(self.a[1],self.a2[1])
            self.assertEqual(self.a[2],self.a2[2])
            
    def test_sum(self):
        with self.context:         
            t1 = self.a.sum()
            
            self.assertEqual(t1,self.sum)
            
            self.assertEqual(self.a[0],self.a2[0])
            self.assertEqual(self.a[1],self.a2[1])
            self.assertEqual(self.a[2],self.a2[2])
            
    def test_dot(self):
        with self.context:
            t1 = self.a.dot(self.b)
            
            self.assertEqual(t1,self.dot)
            
            self.assertEqual(self.a[0],self.a2[0])
            self.assertEqual(self.a[1],self.a2[1])
            self.assertEqual(self.a[2],self.a2[2])
            
            self.assertEqual(self.b[0],self.b2[0])
            self.assertEqual(self.b[1],self.b2[1])
            self.assertEqual(self.b[2],self.b2[2])
    
    def test_max(self):
        with self.context:
            #When updated, this should call self.a.kind
            if self.kind == 'real':
                self.assertEqual(self.a.max(),5)
                
            self.assertEqual(self.a[0],self.a2[0])
            self.assertEqual(self.a[1],self.a2[1])
            self.assertEqual(self.a[2],self.a2[2])
                
    def test_min(self):
        with self.context:
            #When updated, this should call self.a.kind
            if self.kind == 'real':
                self.assertEqual(self.a.min(),1)
                
            self.assertEqual(self.a[0],self.a2[0])
            self.assertEqual(self.a[1],self.a2[1])
            self.assertEqual(self.a[2],self.a2[2])
                
    

def array_test_maker(context,dtype,odtype):
    class tests(tests_base,unittest.TestCase):
        def __init__(self,*args):
            self.context=context
            self.dtype=dtype
            self.odtype=odtype
            if _options['scheme'] == 'cpu':
                self.scheme = None
            elif _options['scheme'] =='cuda':
                self.scheme == pycbc.scheme.CUDAScheme
            else:
                self.scheme == pycbc.scheme.OpenCLScheme            
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
