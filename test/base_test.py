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
This module tests basic operations including when the starting arrays are on different schemes.
The checkCurrentState and checkCurrentState functions must be declared in the child class, because they require
different checks for Arrays, TimeSeries or FrequencySeries.
'''
import pycbc
if pycbc.HAVE_CUDA:
    import pycuda
    import pycuda.gpuarray
if pycbc.HAVE_OPENCL:
    import pyopencl
    import pyopencl.array
import numpy

class checks(object):
    def checkCurrentState(self, inputs, results, places):
        # First importing the correct dtype the context
        if type(pycbc.scheme.mgr.state) is pycbc.scheme.CUDAScheme:
            SchemeArray = pycuda.gpuarray.GPUArray
        elif type(pycbc.scheme.mgr.state) is pycbc.scheme.OpenCLScheme:
            SchemeArray = pyopencl.array.Array
        else:
            SchemeArray = numpy.ndarray
        # Check will specify all the arguments that should be tested, and what their values should be.
        for k in range(len(inputs)):
            val = inputs[k]
            trueval = results[k]
            if isinstance(val,pycbc.types.Array):
                # We will need to check the scheme, and data type
                self.assertTrue(type(val.data)==SchemeArray)
                if type(pycbc.scheme.mgr.state) != pycbc.scheme.DefaultScheme:
                    self.assertTrue(type(val._scheme)==type(pycbc.scheme.mgr.state))
                else:
                    self.assertTrue(val._scheme==None)
                
                # These can be almostequal tests, because the child class should perform
                # the numerical tests, this is just to be sure that it isn't garbled
                self.assertEqual(len(val), len(trueval))
                for i in range(len(trueval)):
                    self.assertAlmostEqual(val[i], trueval[i], places=places)
            else:
                self.assertTrue(val==trueval)

class function_base(checks):
    
    def scheme_test(self, function, inputs, results, places, *args, **kwargs):                    
        # We will need to make lists of inputs, involving the CPU and current scheme, covering all possible permutations
        finalargs = []
        for i in range(pow(2,len(inputs))):
            finalargs.append([])
            for j in range(len(inputs)):
                finalargs[i].append(inputs[j] * 1)
                # Now that it has been added to the list (on the CPU) we will move it to the Current Scheme if necessary
                with self.context:
                    if (i/(pow(2,j))%2 == 1):
                        finalargs[i][j]*=1
            # After these inputs, we will tack on all the other arguments
            for arg in args:
                finalargs[i].append(arg)
        # Now we'll enter the context, run and check all these
        with self.context:
            for i in range(pow(2,len(inputs))):
                function(*finalargs[i],**kwargs)
                self.checkCurrentState(finalargs[i][0:len(inputs)],results,places)
                
    # This function is here so that when different kwargs need to be specified for running on the cpu,
    # the tests can call them correctly.
    def cpu_test(self, function, inputs, results, places, *args, **kwargs):        
        finalargs = []
        for i in range(pow(2,len(inputs))):
            finalargs.append([])
            for j in range(len(inputs)):
                finalargs[i].append(inputs[j] * 1)
                # Now that it has been added to the list (on the CPU) we will move it to the Current Scheme if necessary
                with self.context:
                    if (i/(pow(2,j))%2 == 1):
                        finalargs[i][j]*=1
            # After these inputs, we will tack on all the other arguments
            for arg in args:
                finalargs[i].append(arg)
        
        for i in range(pow(2,len(inputs))):
            function(*finalargs[i],**kwargs)
            self.checkCurrentState(finalargs[i][0:len(inputs)],results,places)

class array_base(checks):
    def setNumbers(self):
        # These are the values that should be used to initialize the test arrays.
        if self.kind == 'real':
            self.alist = [5,3,1]
        else:
            self.alist = [5+1j,3+3j,1+5j]
        if self.okind == 'real':
            self.blist = [10,8,6]
        else:
            self.blist = [10+6j,8+4j,6+2j]
        # And the scalar to test on
        if self.okind == 'real':
            self.scalar = 5
        else:
            self.scalar = 5+2j

        # All the answers are stored here to make it easier to read in the actual tests.
        # Again, it makes a difference whether they are complex or real valued, so there
        # are four sets of possible answers, depending on the dtypes. The fourth value in
        # the list (True or False) lets the checker know whether to use Equal or AlmostEqual
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
        self.min = 1
        self.max = 5
            
    def test_mul(self):
        with self.context:
            # CPU with CPU
            c = self.a1 * self.b1
            self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.mul), self.places)
            # CPU with Other
            c = self.a2 * self.b1
            self.checkCurrentState((self.a2, self.b1, c), (self.alist, self.blist, self.mul), self.places)
            # Current Scheme with Current Scheme
            c = self.a1 * self.b1
            self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.mul), self.places)
            # Current Scheme with CPU
            c = self.a1 * self.b2
            self.checkCurrentState((self.a1, self.b2, c), (self.alist, self.blist, self.mul), self.places)
            # CPU with scalar
            c = self.a3 * self.s
            self.checkCurrentState((self.a3, self.s, c), (self.alist, self.scalar, self.mul_s), self.places)
            # Current Scheme with scalar
            c = self.a1 * self.s
            self.checkCurrentState((self.a1, self.s, c), (self.alist, self.scalar, self.mul_s), self.places)
            
            self.assertRaises(TypeError, self.a1.__mul__, self.bad)
            self.assertRaises(ValueError, self.a1.__mul__, self.bad2)
            
        # Now taking Current Scheme Array and going back to the CPU
        # Current Scheme with Current Scheme
        c = self.a1 * self.b1
        self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.mul), self.places)
        # CPU with Current Scheme
        c = self.a1 * self.b2
        self.checkCurrentState((self.a1, self.b2, c), (self.alist, self.blist, self.mul), self.places)
        # CPU with CPU
        c = self.a1 * self.b1
        self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.mul), self.places)
        # Current Scheme with CPU
        c = self.a2 * self.b1
        self.checkCurrentState((self.a2, self.b1, c), (self.alist, self.blist, self.mul), self.places)
        # Current Scheme with scalar
        c = self.a3 * self.s
        self.checkCurrentState((self.a3, self.s, c), (self.alist, self.scalar, self.mul_s), self.places)
        # CPU with scalar
        c = self.a1 * self.s
        self.checkCurrentState((self.a1, self.s, c), (self.alist, self.scalar, self.mul_s), self.places)
        
    def test_rmul(self):
        with self.context:
            # CPU with CPU
            c = self.a1.__rmul__(self.b1)
            self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.mul), self.places)
            # CPU with Current Scheme
            c = self.a2.__rmul__(self.b1)
            self.checkCurrentState((self.a2, self.b1, c), (self.alist, self.blist, self.mul), self.places)
            # Current Scheme with Current Scheme
            c = self.a1.__rmul__(self.b1)
            self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.mul), self.places)
            # Current Scheme with CPU
            c = self.a1.__rmul__(self.b2)
            self.checkCurrentState((self.a1, self.b2, c), (self.alist, self.blist, self.mul), self.places)
            # CPU with scalar
            c = self.s * self.a3
            self.checkCurrentState((self.a3, self.s, c), (self.alist, self.scalar, self.mul_s), self.places)
            # Current Scheme with scalar
            c = self.s * self.a1
            self.checkCurrentState((self.a1, self.s, c), (self.alist, self.scalar, self.mul_s), self.places)
            
            self.assertRaises(TypeError, self.a1.__rmul__, self.bad)
            self.assertRaises(ValueError, self.a1.__rmul__, self.bad2)
            
        # Now taking Current Scheme Array and going back to the CPU
        # Current Scheme with Current Scheme
        c = self.a1.__rmul__(self.b1)
        self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.mul), self.places)
        # CPU with Current Scheme
        c = self.a1.__rmul__(self.b2)
        self.checkCurrentState((self.a1, self.b2, c), (self.alist, self.blist, self.mul), self.places)
        # CPU with CPU
        c = self.a1.__rmul__(self.b1)
        self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.mul), self.places)
        # Current Scheme with CPU
        c = self.a2.__rmul__(self.b1)
        self.checkCurrentState((self.a2, self.b1, c), (self.alist, self.blist, self.mul), self.places)
        # Current Scheme with scalar
        c = self.s * self.a3
        self.checkCurrentState((self.a3, self.s, c), (self.alist, self.scalar, self.mul_s), self.places)
        # CPU with scalar
        c = self.s * self.a1
        self.checkCurrentState((self.a1, self.s, c), (self.alist, self.scalar, self.mul_s), self.places)
                
    def test_imul(self):
        if not (self.kind == 'real' and self.okind == 'complex'):
            # We need three cs on the cpu
            c1 = self.a1 * 1
            c2 = self.a1 * 1
            c3 = self.a1 * 1
            with self.context:
                # and three on the current scheme
                c4 = self.a1 * 1
                c5 = self.a1 * 1
                c6 = self.a1 * 1
                
                # CPU with CPU
                c1 *= self.b1
                self.checkCurrentState((self.b1, c1), (self.blist, self.mul), self.places)
                # CPU with Current Scheme
                c2 *= self.b1
                self.checkCurrentState((self.b1, c2), (self.blist, self.mul), self.places)
                # Current Scheme with Current Scheme
                c4 *= self.b1
                self.checkCurrentState((self.b1, c4), (self.blist, self.mul), self.places)
                # Current Scheme with CPU
                c5 *= self.b2
                self.checkCurrentState((self.b2, c5), (self.blist, self.mul), self.places)
                # CPU with scalar
                c3 *= self.s
                self.checkCurrentState((self.s, c3), (self.scalar, self.mul_s), self.places)
                # Current Scheme with scalar
                c6 *= self.s
                self.checkCurrentState((self.s, c6), (self.scalar, self.mul_s), self.places)
                
                self.assertRaises(TypeError, self.a1.__imul__, self.bad)
                self.assertRaises(ValueError, self.a1.__imul__, self.bad2)
                
                # We now need to set cs back to the correct values and locations
                c1 = self.a1 * 1
                c2 = self.a1 * 1
                c3 = self.a1 * 1
            c4 = self.a1 * 1
            c5 = self.a1 * 1
            c6 = self.a1 * 1
            # Now taking Current Scheme Array and going back to the CPU
            # Current Scheme with Current Scheme
            c1 *= self.b1
            self.checkCurrentState((self.b1, c1), (self.blist, self.mul), self.places)
            # CPU with Current Scheme
            c4 *= self.b2
            self.checkCurrentState((self.b2, c4), (self.blist, self.mul), self.places)
            # CPU with CPU
            c5 *= self.b1
            self.checkCurrentState((self.b1, c5), (self.blist, self.mul), self.places)
            # Current Scheme with CPU
            c2 *= self.b1
            self.checkCurrentState((self.b1, c2), (self.blist, self.mul), self.places)
            # Current Scheme with scalar
            c3 *= self.s
            self.checkCurrentState((self.s, c3), (self.scalar, self.mul_s), self.places)
            # CPU with scalar
            c6 *= self.s
            self.checkCurrentState((self.s, c6), (self.scalar, self.mul_s), self.places)
            
        else:
            with self.context:
                self.assertRaises(TypeError, self.a1.__imul__,self.s)
                self.assertRaises(TypeError, self.a1.__imul__,self.b1)
            
    def test_add(self):
        with self.context:
            # CPU with CPU
            c = self.a1 + self.b1
            self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.add), self.places)
            # CPU with Other
            c = self.a2 + self.b1
            self.checkCurrentState((self.a2, self.b1, c), (self.alist, self.blist, self.add), self.places)
            # Current Scheme with Current Scheme
            c = self.a1 + self.b1
            self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.add), self.places)
            # Current Scheme with CPU
            c = self.a1 + self.b2
            self.checkCurrentState((self.a1, self.b2, c), (self.alist, self.blist, self.add), self.places)
            # CPU with scalar
            c = self.a3 + self.s
            self.checkCurrentState((self.a3, self.s, c), (self.alist, self.scalar, self.add_s), self.places)
            # Current Scheme with scalar
            c = self.a1 + self.s
            self.checkCurrentState((self.a1, self.s, c), (self.alist, self.scalar, self.add_s), self.places)
            
            self.assertRaises(TypeError, self.a1.__add__, self.bad)
            self.assertRaises(ValueError, self.a1.__add__, self.bad2)
            
        # Now taking Current Scheme Array and going back to the CPU
        # Current Scheme with Current Scheme
        c = self.a1 + self.b1
        self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.add), self.places)
        # CPU with Current Scheme
        c = self.a1 + self.b2
        self.checkCurrentState((self.a1, self.b2, c), (self.alist, self.blist, self.add), self.places)
        # CPU with CPU
        c = self.a1 + self.b1
        self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.add), self.places)
        # Current Scheme with CPU
        c = self.a2 + self.b1
        self.checkCurrentState((self.a2, self.b1, c), (self.alist, self.blist, self.add), self.places)
        # Current Scheme with scalar
        c = self.a3 + self.s
        self.checkCurrentState((self.a3, self.s, c), (self.alist, self.scalar, self.add_s), self.places)
        # CPU with scalar
        c = self.a1 + self.s
        self.checkCurrentState((self.a1, self.s, c), (self.alist, self.scalar, self.add_s), self.places)
        
    def test_radd(self):
        with self.context:
            # CPU with CPU
            c = self.a1.__radd__(self.b1)
            self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.add), self.places)
            # CPU with Current Scheme
            c = self.a2.__radd__(self.b1)
            self.checkCurrentState((self.a2, self.b1, c), (self.alist, self.blist, self.add), self.places)
            # Current Scheme with Current Scheme
            c = self.a1.__radd__(self.b1)
            self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.add), self.places)
            # Current Scheme with CPU
            c = self.a1.__radd__(self.b2)
            self.checkCurrentState((self.a1, self.b2, c), (self.alist, self.blist, self.add), self.places)
            # CPU with scalar
            c = self.s + self.a3
            self.checkCurrentState((self.a3, self.s, c), (self.alist, self.scalar, self.add_s), self.places)
            # Current Scheme with scalar
            c = self.s + self.a1
            self.checkCurrentState((self.a1, self.s, c), (self.alist, self.scalar, self.add_s), self.places)
            
            self.assertRaises(TypeError, self.a1.__radd__, self.bad)
            self.assertRaises(ValueError, self.a1.__radd__, self.bad2)
            
        # Now taking Current Scheme Array and going back to the CPU
        # Current Scheme with Current Scheme
        c = self.a1.__radd__(self.b1)
        self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.add), self.places)
        # CPU with Current Scheme
        c = self.a1.__radd__(self.b2)
        self.checkCurrentState((self.a1, self.b2, c), (self.alist, self.blist, self.add), self.places)
        # CPU with CPU
        c = self.a1.__radd__(self.b1)
        self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.add), self.places)
        # Current Scheme with CPU
        c = self.a2.__radd__(self.b1)
        self.checkCurrentState((self.a2, self.b1, c), (self.alist, self.blist, self.add), self.places)
        # Current Scheme with scalar
        c = self.s + self.a3
        self.checkCurrentState((self.a3, self.s, c), (self.alist, self.scalar, self.add_s), self.places)
        # CPU with scalar
        c = self.s + self.a1
        self.checkCurrentState((self.a1, self.s, c), (self.alist, self.scalar, self.add_s), self.places)
        
    def test_iadd(self):
        if not (self.kind == 'real' and self.okind == 'complex'):
            # We need three cs on the cpu
            c1 = self.a1 * 1
            c2 = self.a1 * 1
            c3 = self.a1 * 1
            with self.context:
                # and three on the current scheme
                c4 = self.a1 * 1
                c5 = self.a1 * 1
                c6 = self.a1 * 1
                
                # CPU with CPU
                c1 += self.b1
                self.checkCurrentState((self.b1, c1), (self.blist, self.add), self.places)
                # CPU with Current Scheme
                c2 += self.b1
                self.checkCurrentState((self.b1, c2), (self.blist, self.add), self.places)
                # Current Scheme with Current Scheme
                c4 += self.b1
                self.checkCurrentState((self.b1, c4), (self.blist, self.add), self.places)
                # Current Scheme with CPU
                c5 += self.b2
                self.checkCurrentState((self.b2, c5), (self.blist, self.add), self.places)
                # CPU with scalar
                c3 += self.s
                self.checkCurrentState((self.s, c3), (self.scalar, self.add_s), self.places)
                # Current Scheme with scalar
                c6 += self.s
                self.checkCurrentState((self.s, c6), (self.scalar, self.add_s), self.places)
                
                self.assertRaises(TypeError, self.a1.__iadd__, self.bad)
                self.assertRaises(ValueError, self.a1.__iadd__, self.bad2)
                
                # We now need to set cs back to the correct values and locations
                c1 = self.a1 * 1
                c2 = self.a1 * 1
                c3 = self.a1 * 1
            c4 = self.a1 * 1
            c5 = self.a1 * 1
            c6 = self.a1 * 1
            # Now taking Current Scheme Array and going back to the CPU
            # Current Scheme with Current Scheme
            c1 += self.b1
            self.checkCurrentState((self.b1, c1), (self.blist, self.add), self.places)
            # CPU with Current Scheme
            c4 += self.b2
            self.checkCurrentState((self.b2, c4), (self.blist, self.add), self.places)
            # CPU with CPU
            c5 += self.b1
            self.checkCurrentState((self.b1, c5), (self.blist, self.add), self.places)
            # Current Scheme with CPU
            c2 += self.b1
            self.checkCurrentState((self.b1, c2), (self.blist, self.add), self.places)
            # Current Scheme with scalar
            c3 += self.s
            self.checkCurrentState((self.s, c3), (self.scalar, self.add_s), self.places)
            # CPU with scalar
            c6 += self.s
            self.checkCurrentState((self.s, c6), (self.scalar, self.add_s), self.places)
            
        else:
            with self.context:
                self.assertRaises(TypeError, self.a1.__iadd__,self.s)
                self.assertRaises(TypeError, self.a1.__iadd__,self.b1)
    
    def test_div(self):
        with self.context:
            # CPU with CPU
            c = self.a1 / self.b1
            self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.div), self.places)
            # CPU with Current Scheme
            c = self.a2 / self.b1
            self.checkCurrentState((self.a2, self.b1, c), (self.alist, self.blist, self.div), self.places)
            # Current Scheme with Current Scheme
            c = self.a1 / self.b1
            self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.div), self.places)
            # Current Scheme with CPU
            c = self.a1 / self.b2
            self.checkCurrentState((self.a1, self.b2, c), (self.alist, self.blist, self.div), self.places)
            # CPU with scalar
            c = self.a3 / self.s
            self.checkCurrentState((self.a3, self.s, c), (self.alist, self.scalar, self.div_s), self.places)
            # Current Scheme with scalar
            c = self.a1 / self.s
            self.checkCurrentState((self.a1, self.s, c), (self.alist, self.scalar, self.div_s), self.places)
            
            self.assertRaises(TypeError, self.a1.__div__, self.bad)
            self.assertRaises(ValueError, self.a1.__div__, self.bad2)
            
        # Now taking Current Scheme Array and going back to the CPU
        # Current Scheme with Current Scheme
        c = self.a1 / self.b1
        self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.div), self.places)
        # CPU with Current Scheme
        c = self.a1 / self.b2
        self.checkCurrentState((self.a1, self.b2, c), (self.alist, self.blist, self.div), self.places)
        # CPU with CPU
        c = self.a1 / self.b1
        self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.div), self.places)
        # Current Scheme with CPU
        c = self.a2 / self.b1
        self.checkCurrentState((self.a2, self.b1, c), (self.alist, self.blist, self.div), self.places)
        # Current Scheme with scalar
        c = self.a3 / self.s
        self.checkCurrentState((self.a3, self.s, c), (self.alist, self.scalar, self.div_s), self.places)
        # CPU with scalar
        c = self.a1 / self.s
        self.checkCurrentState((self.a1, self.s, c), (self.alist, self.scalar, self.div_s), self.places)
        
    def test_rdiv(self):
        with self.context:
            # CPU with CPU
            c = self.a1.__rdiv__(self.b1)
            self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.rdiv), self.places)
            # CPU with Current Scheme
            c = self.a2.__rdiv__(self.b1)
            self.checkCurrentState((self.a2, self.b1, c), (self.alist, self.blist, self.rdiv), self.places)
            # Current Scheme with Current Scheme
            c = self.a1.__rdiv__(self.b1)
            self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.rdiv), self.places)
            # Current Scheme with CPU
            c = self.a1.__rdiv__(self.b2)
            self.checkCurrentState((self.a1, self.b2, c), (self.alist, self.blist, self.rdiv), self.places)
            # CPU with scalar
            c = self.s / self.a3
            self.checkCurrentState((self.a3, self.s, c), (self.alist, self.scalar, self.rdiv_s), self.places)
            # Current Scheme with scalar
            c = self.s / self.a1
            self.checkCurrentState((self.a1, self.s, c), (self.alist, self.scalar, self.rdiv_s), self.places)
            
            self.assertRaises(TypeError, self.a1.__rdiv__, self.bad)
            self.assertRaises(ValueError, self.a1.__rdiv__, self.bad2)
            
        # Now taking Current Scheme Array and going back to the CPU
        # Current Scheme with Current Scheme
        c = self.a1.__rdiv__(self.b1)
        self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.rdiv), self.places)
        # CPU with Current Scheme
        c = self.a1.__rdiv__(self.b2)
        self.checkCurrentState((self.a1, self.b2, c), (self.alist, self.blist, self.rdiv), self.places)
        # CPU with CPU
        c = self.a1.__rdiv__(self.b1)
        self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.rdiv), self.places)
        # Current Scheme with CPU
        c = self.a2.__rdiv__(self.b1)
        self.checkCurrentState((self.a2, self.b1, c), (self.alist, self.blist, self.rdiv), self.places)
        # Current Scheme with scalar
        c = self.s / self.a3
        self.checkCurrentState((self.a3, self.s, c), (self.alist, self.scalar, self.rdiv_s), self.places)
        # CPU with scalar
        c = self.s / self.a1
        self.checkCurrentState((self.a1, self.s, c), (self.alist, self.scalar, self.rdiv_s), self.places)
        
    def test_idiv(self):
        if not (self.kind == 'real' and self.okind == 'complex'):
            # We need three cs on the cpu
            c1 = self.a1 * 1
            c2 = self.a1 * 1
            c3 = self.a1 * 1
            with self.context:
                # and three on the current scheme
                c4 = self.a1 * 1
                c5 = self.a1 * 1
                c6 = self.a1 * 1
                
                # CPU with CPU
                c1 /= self.b1
                self.checkCurrentState((self.b1, c1), (self.blist, self.div), self.places)
                # CPU with Current Scheme
                c2 /= self.b1
                self.checkCurrentState((self.b1, c2), (self.blist, self.div), self.places)
                # Current Scheme with Current Scheme
                c4 /= self.b1
                self.checkCurrentState((self.b1, c4), (self.blist, self.div), self.places)
                # Current Scheme with CPU
                c5 /= self.b2
                self.checkCurrentState((self.b2, c5), (self.blist, self.div), self.places)
                # CPU with scalar
                c3 /= self.s
                self.checkCurrentState((self.s, c3), (self.scalar, self.div_s), self.places)
                # Current Scheme with scalar
                c6 /= self.s
                self.checkCurrentState((self.s, c6), (self.scalar, self.div_s), self.places)
                
                self.assertRaises(TypeError, self.a1.__idiv__, self.bad)
                self.assertRaises(ValueError, self.a1.__idiv__, self.bad2)
                
                # We now need to set cs back to the correct values and locations
                c1 = self.a1 * 1
                c2 = self.a1 * 1
                c3 = self.a1 * 1
            c4 = self.a1 * 1
            c5 = self.a1 * 1
            c6 = self.a1 * 1
            # Now taking Current Scheme Array and going back to the CPU
            # Current Scheme with Current Scheme
            c1 /= self.b1
            self.checkCurrentState((self.b1, c1), (self.blist, self.div), self.places)
            # CPU with Current Scheme
            c4 /= self.b2
            self.checkCurrentState((self.b2, c4), (self.blist, self.div), self.places)
            # CPU with CPU
            c5 /= self.b1
            self.checkCurrentState((self.b1, c5), (self.blist, self.div), self.places)
            # Current Scheme with CPU
            c2 /= self.b1
            self.checkCurrentState((self.b1, c2), (self.blist, self.div), self.places)
            # Current Scheme with scalar
            c3 /= self.s
            self.checkCurrentState((self.s, c3), (self.scalar, self.div_s), self.places)
            # CPU with scalar
            c6 /= self.s
            self.checkCurrentState((self.s, c6), (self.scalar, self.div_s), self.places)
            
        else:
            with self.context:
                self.assertRaises(TypeError, self.a1.__idiv__,self.s)
                self.assertRaises(TypeError, self.a1.__idiv__,self.b1)
            
    def test_sub(self):
        with self.context:
            # CPU with CPU
            c = self.a1 - self.b1
            self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.sub), self.places)
            # CPU with Current Scheme
            c = self.a2 - self.b1
            self.checkCurrentState((self.a2, self.b1, c), (self.alist, self.blist, self.sub), self.places)
            # Current Scheme with Current Scheme
            c = self.a1 - self.b1
            self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.sub), self.places)
            # Current Scheme with CPU
            c = self.a1 - self.b2
            self.checkCurrentState((self.a1, self.b2, c), (self.alist, self.blist, self.sub), self.places)
            # CPU with scalar
            c = self.a3 - self.s
            self.checkCurrentState((self.a3, self.s, c), (self.alist, self.scalar, self.sub_s), self.places)
            # Current Scheme with scalar
            c = self.a1 - self.s
            self.checkCurrentState((self.a1, self.s, c), (self.alist, self.scalar, self.sub_s), self.places)
            
            self.assertRaises(TypeError, self.a1.__sub__, self.bad)
            self.assertRaises(ValueError, self.a1.__sub__, self.bad2)
            
        # Now taking Current Scheme Array and going back to the CPU
        # Current Scheme with Current Scheme
        c = self.a1 - self.b1
        self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.sub), self.places)
        # CPU with Current Scheme
        c = self.a1 - self.b2
        self.checkCurrentState((self.a1, self.b2, c), (self.alist, self.blist, self.sub), self.places)
        # CPU with CPU
        c = self.a1 - self.b1
        self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.sub), self.places)
        # Current Scheme with CPU
        c = self.a2 - self.b1
        self.checkCurrentState((self.a2, self.b1, c), (self.alist, self.blist, self.sub), self.places)
        # Current Scheme with scalar
        c = self.a3 - self.s
        self.checkCurrentState((self.a3, self.s, c), (self.alist, self.scalar, self.sub_s), self.places)
        # CPU with scalar
        c = self.a1 - self.s
        self.checkCurrentState((self.a1, self.s, c), (self.alist, self.scalar, self.sub_s), self.places)
        
    def test_rsub(self):
        with self.context:
            # CPU with CPU
            c = self.a1.__rsub__(self.b1)
            self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.rsub), self.places)
            # CPU with Current Scheme
            c = self.a2.__rsub__(self.b1)
            self.checkCurrentState((self.a2, self.b1, c), (self.alist, self.blist, self.rsub), self.places)
            # Current Scheme with Current Scheme
            c = self.a1.__rsub__(self.b1)
            self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.rsub), self.places)
            # Current Scheme with CPU
            c = self.a1.__rsub__(self.b2)
            self.checkCurrentState((self.a1, self.b2, c), (self.alist, self.blist, self.rsub), self.places)
            # CPU with scalar
            c = self.s - self.a3
            self.checkCurrentState((self.a3, self.s, c), (self.alist, self.scalar, self.rsub_s), self.places)
            # Current Scheme with scalar
            c = self.s - self.a1
            self.checkCurrentState((self.a1, self.s, c), (self.alist, self.scalar, self.rsub_s), self.places)
            
            self.assertRaises(TypeError, self.a1.__rsub__, self.bad)
            self.assertRaises(ValueError, self.a1.__rsub__, self.bad2)
            
        # Now taking Current Scheme Array and going back to the CPU
        # Current Scheme with Current Scheme
        c = self.a1.__rsub__(self.b1)
        self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.rsub), self.places)
        # CPU with Current Scheme
        c = self.a1.__rsub__(self.b2)
        self.checkCurrentState((self.a1, self.b2, c), (self.alist, self.blist, self.rsub), self.places)
        # CPU with CPU
        c = self.a1.__rsub__(self.b1)
        self.checkCurrentState((self.a1, self.b1, c), (self.alist, self.blist, self.rsub), self.places)
        # Current Scheme with CPU
        c = self.a2.__rsub__(self.b1)
        self.checkCurrentState((self.a2, self.b1, c), (self.alist, self.blist, self.rsub), self.places)
        # Current Scheme with scalar
        c = self.s - self.a3
        self.checkCurrentState((self.a3, self.s, c), (self.alist, self.scalar, self.rsub_s), self.places)
        # CPU with scalar
        c = self.s - self.a1
        self.checkCurrentState((self.a1, self.s, c), (self.alist, self.scalar, self.rsub_s), self.places)
        
    def test_isub(self):
        if not (self.kind == 'real' and self.okind == 'complex'):
            # We need three cs on the cpu
            c1 = self.a1 * 1
            c2 = self.a1 * 1
            c3 = self.a1 * 1
            with self.context:
                # and three on the current scheme
                c4 = self.a1 * 1
                c5 = self.a1 * 1
                c6 = self.a1 * 1
                
                # CPU with CPU
                c1 -= self.b1
                self.checkCurrentState((self.b1, c1), (self.blist, self.sub), self.places)
                # CPU with Current Scheme
                c2 -= self.b1
                self.checkCurrentState((self.b1, c2), (self.blist, self.sub), self.places)
                # Current Scheme with Current Scheme
                c4 -= self.b1
                self.checkCurrentState((self.b1, c4), (self.blist, self.sub), self.places)
                # Current Scheme with CPU
                c5 -= self.b2
                self.checkCurrentState((self.b2, c5), (self.blist, self.sub), self.places)
                # CPU with scalar
                c3 -= self.s
                self.checkCurrentState((self.s, c3), (self.scalar, self.sub_s), self.places)
                # Current Scheme with scalar
                c6 -= self.s
                self.checkCurrentState((self.s, c6), (self.scalar, self.sub_s), self.places)
                
                self.assertRaises(TypeError, self.a1.__isub__, self.bad)
                self.assertRaises(ValueError, self.a1.__isub__, self.bad2)
                
                # We now need to set cs back to the correct values and locations
                c1 = self.a1 * 1
                c2 = self.a1 * 1
                c3 = self.a1 * 1
            c4 = self.a1 * 1
            c5 = self.a1 * 1
            c6 = self.a1 * 1
            # Now taking Current Scheme Array and going back to the CPU
            # Current Scheme with Current Scheme
            c1 -= self.b1
            self.checkCurrentState((self.b1, c1), (self.blist, self.sub), self.places)
            # CPU with Current Scheme
            c4 -= self.b2
            self.checkCurrentState((self.b2, c4), (self.blist, self.sub), self.places)
            # CPU with CPU
            c5 -= self.b1
            self.checkCurrentState((self.b1, c5), (self.blist, self.sub), self.places)
            # Current Scheme with CPU
            c2 -= self.b1
            self.checkCurrentState((self.b1, c2), (self.blist, self.sub), self.places)
            # Current Scheme with scalar
            c3 -= self.s
            self.checkCurrentState((self.s, c3), (self.scalar, self.sub_s), self.places)
            # CPU with scalar
            c6 -= self.s
            self.checkCurrentState((self.s, c6), (self.scalar, self.sub_s), self.places)
            
        else:
            with self.context:
                self.assertRaises(TypeError, self.a1.__isub__,self.s)
                self.assertRaises(TypeError, self.a1.__isub__,self.b1)
        
    def test_pow(self):
        with self.context:
            self.b1 *= 1
            # From CPU
            c1 = self.a1 ** 2
            c2 = self.a2 ** -1.5
            
            self.checkCurrentState((self.a1, c1), (self.alist, self.pow1), self.places)
            self.checkCurrentState((self.a2, c2), (self.alist, self.pow2), self.places)
        # From Current Scheme
        c1 = self.a1 ** 2
        c2 = self.a2 ** -1.5
        
        self.checkCurrentState((self.a1, c1), (self.alist, self.pow1), self.places)
        self.checkCurrentState((self.a2, c2), (self.alist, self.pow2), self.places)
        
    def test_abs(self):
        # We want to check that absolute value behaves correctly no matter
        # what quadrant it's in. First we will check with cpu arrays
        t1 = self.a1 * 1
        t2 = self.a1 * -1
        t3 = self.a1 * 1j
        t4 = self.a1 * -1j
        with self.context:
            self.b1 *= 1
            c1 = abs(t1)
            c2 = abs(t2)
            c3 = abs(t3)
            c4 = abs(t4)
            
            self.checkCurrentState((t1, c1), (self.alist,self.abs), self.places)
            self.checkCurrentState((t1, c2), (self.alist,self.abs), self.places)
            self.checkCurrentState((t1, c3), (self.alist,self.abs), self.places)
            self.checkCurrentState((t1, c4), (self.alist,self.abs), self.places)
            # Now coming from the current scheme
            c1 = abs(t1)
            c2 = abs(t2)
            c3 = abs(t3)
            c4 = abs(t4)
            
            self.checkCurrentState((t1, c1), (self.alist,self.abs), self.places)
            self.checkCurrentState((t1, c2), (self.alist,self.abs), self.places)
            self.checkCurrentState((t1, c3), (self.alist,self.abs), self.places)
            self.checkCurrentState((t1, c4), (self.alist,self.abs), self.places)
        # Now taking abs of the current scheme on the CPU
        c1 = abs(t1)
        c2 = abs(t2)
        c3 = abs(t3)
        c4 = abs(t4)
        
        self.checkCurrentState((t1, c1), (self.alist,self.abs), self.places)
        self.checkCurrentState((t1, c2), (self.alist,self.abs), self.places)
        self.checkCurrentState((t1, c3), (self.alist,self.abs), self.places)
        self.checkCurrentState((t1, c4), (self.alist,self.abs), self.places)
        #And finally, from the CPU to the CPU
        c1 = abs(t1)
        c2 = abs(t2)
        c3 = abs(t3)
        c4 = abs(t4)
        
        self.checkCurrentState((t1, c1), (self.alist,self.abs), self.places)
        self.checkCurrentState((t1, c2), (self.alist,self.abs), self.places)
        self.checkCurrentState((t1, c3), (self.alist,self.abs), self.places)
        self.checkCurrentState((t1, c4), (self.alist,self.abs), self.places)
            
        
    def test_real(self):
        with self.context:
            # From CPU
            c = self.a1.real()
            self.checkCurrentState((self.a1, c),(self.alist,self.real), self.places)
            
            # From Current Scheme
            c = self.a1.real()
            self.checkCurrentState((self.a1, c),(self.alist,self.real), self.places)
        # Now on the CPU, from Current Scheme
        c = self.a1.real()
        self.checkCurrentState((self.a1, c),(self.alist,self.real), self.places)
        # And finally CPU on the CPU
        c = self.a1.real()
        self.checkCurrentState((self.a1, c),(self.alist,self.real), self.places)
            

        
    def test_imag(self):
        with self.context:
            # From CPU
            c = self.a1.imag()
            self.checkCurrentState((self.a1, c),(self.alist,self.imag), self.places)
            
            # From Current Scheme
            c = self.a1.imag()
            self.checkCurrentState((self.a1, c),(self.alist,self.imag), self.places)
        # Now on the CPU, from Current Scheme
        c = self.a1.imag()
        self.checkCurrentState((self.a1, c),(self.alist,self.imag), self.places)
        # And finally CPU on the CPU
        c = self.a1.imag()
        self.checkCurrentState((self.a1, c),(self.alist,self.imag), self.places)
        
    def test_conj(self):
        with self.context:
            # From CPU
            c = self.a1.conj()
            self.checkCurrentState((self.a1, c),(self.alist,self.conj), self.places)
            
            # From Current Scheme
            c = self.a1.conj()
            self.checkCurrentState((self.a1, c),(self.alist,self.conj), self.places)
        # Now on the CPU, from Current Scheme
        c = self.a1.conj()
        self.checkCurrentState((self.a1, c),(self.alist,self.conj), self.places)
        # And finally CPU on the CPU
        c = self.a1.conj()
        self.checkCurrentState((self.a1, c),(self.alist,self.conj), self.places)
            
    def test_sum(self):
        with self.context:
            # From CPU
            c = self.a1.sum()
            self.checkCurrentState((self.a1, c),(self.alist,self.sum), self.places)
            
            # From Current Scheme
            c = self.a1.sum()
            self.checkCurrentState((self.a1, c),(self.alist,self.sum), self.places)
        # Now on the CPU, from Current Scheme
        c = self.a1.sum()
        self.checkCurrentState((self.a1, c),(self.alist,self.sum), self.places)
        # And finally CPU on the CPU
        c = self.a1.sum()
        self.checkCurrentState((self.a1, c),(self.alist,self.sum), self.places)
            
    def test_dot(self):
        with self.context:
            # CPU with CPU
            c = self.a1.dot(self.b1)
            self.checkCurrentState((self.a1, self.b1, c),(self.alist, self.blist, self.dot), self.places)
            # CPU with Current Scheme
            c = self.a2.dot(self.b1)
            self.checkCurrentState((self.a2, self.b1, c),(self.alist, self.blist, self.dot), self.places)
            # Current Scheme with Current Scheme
            c = self.a1.dot(self.b1)
            self.checkCurrentState((self.a1, self.b1, c),(self.alist, self.blist, self.dot), self.places)
            # Current Scheme with CPU
            c = self.a1.dot(self.b2)
            self.checkCurrentState((self.a1, self.b2, c),(self.alist, self.blist, self.dot), self.places)
            
            self.assertRaises(TypeError, self.a1.dot, self.bad)
            self.assertRaises(ValueError, self.a1.dot, self.bad2)
            
        # Now taking Current Scheme Array and going back to the CPU
        # Current Scheme with Current Scheme
        c = self.a1.dot(self.b1)
        self.checkCurrentState((self.a1, self.b1, c),(self.alist, self.blist, self.dot), self.places)
        # CPU with Current Scheme
        c = self.a1.dot(self.b2)
        self.checkCurrentState((self.a1, self.b2, c),(self.alist, self.blist, self.dot), self.places)
        # CPU with CPU
        c = self.a1.dot(self.b1)
        self.checkCurrentState((self.a1, self.b1, c),(self.alist, self.blist, self.dot), self.places)
        # Current Scheme with CPU
        c = self.a2.dot(self.b1)
        self.checkCurrentState((self.a2, self.b1, c),(self.alist, self.blist, self.dot), self.places)
    
    def test_max(self):
        if self.kind == 'real':
            with self.context:
                # From CPU
                c = self.a1.max()
                self.checkCurrentState((self.a1, c),(self.alist,self.max), self.places)
                # From Current Scheme
                c = self.a1.max()
                self.checkCurrentState((self.a1, c),(self.alist,self.max), self.places)
            # From Current Scheme
            c = self.a1.max()
            self.checkCurrentState((self.a1, c),(self.alist,self.max), self.places)
            # From CPU
            c = self.a1.max()
            self.checkCurrentState((self.a1, c),(self.alist,self.max), self.places)

    def test_min(self):
        if self.kind == 'real':
            with self.context:
                self.b1 *= 1
                # From CPU
                c = self.a1.min()
                self.checkCurrentState((self.a1, c),(self.alist,self.min), self.places)
                # From Current Scheme
                c = self.a1.min()
                self.checkCurrentState((self.a1, c),(self.alist,self.min), self.places)
            # From Current Scheme
            c = self.a1.min()
            self.checkCurrentState((self.a1, c),(self.alist,self.min), self.places)
            # From CPU
            c = self.a1.min()
            self.checkCurrentState((self.a1, c),(self.alist,self.min), self.places)
                
