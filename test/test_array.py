# Copyright (C) 2012  Alex Nitz, Andrew Miller, Tito Dal Canton, Josh Willis
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
from utils import array_base, parse_args_all_schemes, simple_exit
import sys
import os
import tempfile

_scheme, _context = parse_args_all_schemes("Array")

# By importing the current schemes array type, it will make it
# easier to check the  array types later
if _scheme == 'cuda':
    import pycuda
    import pycuda.gpuarray
    from pycuda.gpuarray import GPUArray as SchemeArray
elif _scheme == 'cpu':
    from pycbc.types.aligned import ArrayWithAligned as SchemeArray

from pycbc.types.aligned import ArrayWithAligned as CPUArray

# ********************GENERIC ARRAY TESTS ***********************

class ArrayTestBase(array_base,unittest.TestCase):
    def setUp(self):
        self.scheme = _scheme
        self.context = _context
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
            self.tol = 1e-5
        # Number of decimal places to compare for double precision
        else:
            self.places = 13
            self.tol = 1e-13

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
        # We need to tell the arithmetic test functions what our type is:
        self.type = Array
        # and for Array, there are no additional keyword arguments needed:
        self.kwds = {}

        # Now that the kinds are set, we need to call our parent method to set up all the
        # inputs and answers for our functions
        self.setNumbers()

        # The above call created instances for all of our inputs and various correct
        # outputs.  But we make a copy of the scalar to check later.
        self.s = self.scalar

        # Finally, we want to have an array that we shouldn't be able to operate on,
        # first because the precision is wrong, and seconds because the length is wrong.
        self.bad = Array([1,1,1],dtype = self.other_precision[self.odtype])
        self.bad2 = Array([1,1,1,1], dtype = self.dtype)

    def test_set(self):
        c = self.a * 1
        with self.context:
            # First we will check that get works properly for all
            # the different python syntaxes
            self.assertTrue(self.a[:][0] == self.alist[0:3][0])
            self.assertTrue(self.a[:][1] == self.alist[0:3][1])
            self.assertTrue(self.a[:][2] == self.alist[0:3][2])
            self.assertRaises(IndexError,self.a[:].__getitem__,3)
            self.assertTrue(self.a[-1] ==self.alist[2])
            self.assertTrue(self.a[-2] == self.alist[1])
            self.assertTrue(self.a[1:2][0] == self.alist[1])
            self.assertRaises(IndexError,self.a[1:2].__getitem__,1)
            self.assertTrue(self.a[:-1][0] == self.alist[0:2][0])
            self.assertTrue(self.a[:-1][1] == self.alist[0:2][1])
            self.assertTrue(self.a[-1:][0] == self.alist[2])
            self.assertRaises(IndexError, self.a.__getitem__, 3)
            self.assertRaises(IndexError, self.a.__getitem__, -4)

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
                self.assertTrue(type(out1._scheme) == type(self.context))
                self.assertTrue(type(out1._data) is SchemeArray)
                self.assertEqual(out1[0],5)
                self.assertEqual(out1[1],3)
                self.assertEqual(out1[2],1)
                self.assertTrue(out1.dtype==self.dtype)

                self.assertTrue(type(out2._scheme) == type(self.context))
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
            self.assertTrue(type(out3._scheme) == type(self.context))
            self.assertTrue(type(out3._data) is SchemeArray)
            self.assertEqual(out3[0],5)
            self.assertEqual(out3[1],3)
            self.assertEqual(out3[2],1)
            self.assertTrue(out3.dtype==self.odtype)

            # Check for copy=false
            # On the CPU, this should be possible
            in3 = numpy.array([5,3,1],dtype=self.dtype)
            if self.scheme == 'cpu':
                out4 = Array(in3,copy=False)
                in3 += 1

                self.assertTrue(out4.dtype==self.dtype)
                self.assertTrue(type(out4._scheme) == type(self.context))
                self.assertEqual(out4[0],6)
                self.assertEqual(out4[1],4)
                self.assertEqual(out4[2],2)

            # If we're in different scheme, this should raise an error
            else:
                self.assertRaises(TypeError, Array, in3, copy=False)

            # We also need to check initialization using GPU arrays
            if self.scheme == 'cuda':
                in4 = pycuda.gpuarray.zeros(3,self.dtype)
            if self.scheme != 'cpu':
                out4 = Array(in4, copy=False)
                in4 += 1
                self.assertTrue(type(out4._scheme) == type(self.context))
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
            self.assertTrue(type(out5._scheme) == type(self.context))
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

        if self.scheme != 'cpu':
            self.assertRaises(TypeError, Array, in4, copy=False)
        self.assertRaises(TypeError, Array, in5, copy=False)


    def test_array_init(self):
        # this array is made outside the context so we can check that an error is raised when copy = false in a GPU scheme
        cpuarray = Array([1,2,3])
        with self.context:
            in1 = Array([5,3,1],dtype=self.odtype)
            in2 = Array([5,3,1],dtype=self.other_precision[self.odtype])
            self.assertTrue(type(in1._scheme) == type(self.context))
            self.assertTrue(type(in1._data) is SchemeArray)
            self.assertTrue(type(in2._scheme) == type(self.context))
            self.assertTrue(type(in2._data) is SchemeArray)
            # We don't want to cast complex as real
            if not (self.kind=='real' and self.okind == 'complex'):
                # First we must check that the dtype is correct when specified
                out1 = Array(in1, dtype=self.dtype)
                out2 = Array(in2, dtype=self.dtype)
                # to be sure that it is copied
                in1 += 1
                in2 += 1

                self.assertTrue(type(out1._scheme) == type(self.context))
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

                self.assertTrue(type(out2._scheme) == type(self.context))
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

            self.assertTrue(type(out3._scheme) == type(self.context))
            self.assertTrue(type(out3._data) is SchemeArray)
            self.assertEqual(out3[0],5)
            self.assertEqual(out3[1],3)
            self.assertEqual(out3[2],1)
            self.assertTrue(out3.dtype==self.odtype)

            # We should also be able to create from a CPU Array
            out4 = Array(cpuarray, dtype=self.dtype)

            self.assertTrue(type(out4._scheme) == type(self.context))
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

            self.assertTrue(type(out5._scheme) == type(self.context))
            self.assertTrue(type(out5._data) is SchemeArray)
            self.assertEqual(out5[0],6)
            self.assertEqual(out5[1],4)
            self.assertEqual(out5[2],2)
            self.assertTrue(out5.dtype==self.dtype)

            if self.scheme != 'cpu':
                self.assertRaises(TypeError,Array,cpuarray,copy=False)

            # Just checking that we can make an empty array correctly
            empty = Array(numpy.array([]))
            out7 = Array(empty)
            self.assertTrue(out7.dtype==numpy.float64)
            self.assertRaises(IndexError, out7.__getitem__,0)

        # Also checking that a cpu array can't be made out of another scheme without copying
        if self.scheme != 'cpu':
            self.assertRaises(TypeError, Array, out4, copy=False)
            out6 = Array(out4, dtype=self.dtype)
            self.assertTrue(type(out6._scheme) == CPUScheme)
            self.assertTrue(type(out6._data) is CPUArray)
            self.assertEqual(out6[0],1)
            self.assertEqual(out6[1],2)
            self.assertEqual(out6[2],3)
            self.assertTrue(out6.dtype==self.dtype)

    def test_take(self):
        with self.context:
            if self.kind == 'real':
                a = Array([1,2,3,4,5,6], dtype=self.dtype)

            if self.kind == 'complex':
                a = Array([1+2j, 2+3j, 3+4j, 4+5j, 5+6j, 6+7j], dtype=self.dtype)

            i = numpy.array([0,4,2], dtype=numpy.int64)
            b = a.take(i)
            self.assertEqual(b[0], a[0])
            self.assertEqual(b[1], a[4])
            self.assertEqual(b[2], a[2])

    def test_abs_max_loc(self):
        with self.context:
            if self.kind == 'real':
                a = Array([-1,2,3,4,5,-6], dtype=self.dtype)
                v = abs(-6)

            if self.kind == 'complex':
                a = Array([1+2j, 2+3j, -3+4j, 4+5j, 5+6j, -6+7j], dtype=self.dtype)
                v = abs(-6+7j)

            m, l = a.abs_max_loc()
            self.assertAlmostEqual(m, v, places=5)
            self.assertEqual(l, 5)

    def test_clear(self):
        with self.context:
            if self.kind == 'real':
                a = Array([1,2,3,4,5,6], dtype=self.dtype)

            if self.kind == 'complex':
                a = Array([1+2j, 2+3j, 3+4j, 4+5j, 5+6j, 6+7j], dtype=self.dtype)

            a.clear()
            for i in range(len(a)):
                self.assertEqual(a[i], 0)


    def test_max_loc(self):
        with self.context:
            if self.kind == 'real':
                a = Array([1,2,3,4,5,6], dtype=self.dtype)

                m, l = a.max_loc()
                self.assertEqual(m, 6)
                self.assertEqual(l, 5)

    def test_list_init(self):
        with self.context:
            # When specified
            out1 = Array([5,3,1], dtype=self.dtype)

            self.assertTrue(type(out1._scheme) == type(self.context))
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

                self.assertTrue(type(out2._scheme) == type(self.context))
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

            self.assertTrue(type(out3._scheme) == type(self.context))
            self.assertTrue(type(out3._data) is SchemeArray)
            self.assertEqual(out3[0],5)
            self.assertEqual(out3[1],3)
            self.assertEqual(out3[2],1)
            self.assertTrue(out3.dtype==numpy.float64)

            out4 = Array([5+0j,3+0j,1+0j])

            self.assertTrue(type(out4._scheme) == type(self.context))
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

            self.assertTrue(type(out5._scheme) == type(self.context))
            self.assertTrue(type(out5._data) is SchemeArray)
            self.assertEqual(out5[0],0)
            self.assertEqual(out5[1],0)
            self.assertEqual(out5[2],0)
            self.assertTrue(out5.dtype == self.dtype)

            self.assertTrue(type(out6._scheme) == type(self.context))
            self.assertTrue(type(out6._data) is SchemeArray)
            self.assertEqual(out6[0],0)
            self.assertEqual(out6[1],0)
            self.assertEqual(out6[2],0)
            self.assertTrue(out6.dtype == numpy.float64)

            self.assertRaises(TypeError,Array,[1,2,3],copy=False)

    def test_save(self):
        with self.context:
            # make temporary file paths
            temp_file = tempfile.NamedTemporaryFile()
            temp_path_npy = temp_file.name + '.npy'
            temp_path_txt = temp_file.name + '.txt'
            # make a test array
            a_numpy = numpy.arange(100, dtype=self.dtype)
            a = Array(a_numpy)
            # test saving to Numpy array
            a.save(temp_path_npy)
            b = numpy.load(temp_path_npy)
            self.assertEqual(b.shape, a_numpy.shape)
            self.assertEqual(numpy.abs(b - a_numpy).max(), 0)
            os.remove(temp_path_npy)
            # test saving to text file
            a.save(temp_path_txt)
            b = numpy.loadtxt(temp_path_txt)
            if a.kind == 'complex':
                self.assertEqual(b.shape, (a_numpy.shape[0], 2))
                b = b[:,0] + 1j * b[:,1]
            elif a.kind == 'real':
                self.assertEqual(b.shape, a_numpy.shape)
            self.assertEqual(numpy.abs(b - a_numpy).max(), 0)
            os.remove(temp_path_txt)

def array_test_maker(dtype,odtype):
    class tests(ArrayTestBase):
        def __init__(self,*args):
            self.dtype = dtype
            self.odtype = odtype
            unittest.TestCase.__init__(self,*args)
    tests.__name__ = _scheme + " " + dtype.__name__ + " with " + odtype.__name__
    return tests

types = [ (float32,[float32,complex64]), (float64,[float64,complex128]),
        (complex64,[complex64,float32]), (complex128,[float64,complex128]) ]

suite = unittest.TestSuite()

ind = 0
for ty,oktype in types:
    for ot in oktype:
        na = 'test' + str(ind)
        vars()[na] = array_test_maker(ty,ot)
        suite.addTest(unittest.TestLoader().loadTestsFromTestCase(vars()[na]))
        ind += 1



# TODO More specific array tests (instatiation, failure modes, type conversion, etc)


if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
