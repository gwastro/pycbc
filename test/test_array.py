



import pycbc

if pycbc.have_cuda:
    import pycuda
    import pycuda.gpuarray

if pycbc.have_opencl:
    import pyopencl
    import pyopencl.array

import unittest
from pycbc.array import Array
from pycbc.processingcontext import *
import numpy 

# ********************GENERIC ARRAY TESTS ***********************

class tests_base(object):

    def setUp(self):
        pass

    def test_init(self):
        a = Array([5],context=self.context,dtype=self.dtype)
        b = numpy.array([5])
        c = Array(b,context=self.context,dtype=self.dtype)
        self.assertEqual(a[0],b[0])
        self.assertEqual(c[0],b[0])
        
    def test_mul(self):
        a = Array([5],context=self.context,dtype=self.dtype)
        b = a * 2
        self.assertEqual(b[0],10)
        
    def test_rmul(self):
        a = Array([5],context=self.context,dtype=self.dtype)
        b = 2 * a
        self.assertEqual(b[0],10)
        
    def test_imul(self):
        a = Array([5],context=self.context,dtype=self.dtype)
        a *= 2
        self.assertEqual(a[0],10)
        
    def test_add(self):
        a = Array([5],context=self.context,dtype=self.dtype)
        b = a + 5
        self.assertEqual(b[0],10)
        
    def test_radd(self):
        a = Array([5],context=self.context,dtype=self.dtype)
        b = 5 + a
        self.assertEqual(b[0],10)
        
    def test_iadd(self):
        a = Array([5],context=self.context,dtype=self.dtype)
        a += 5
        self.assertEqual(a[0],10)

    def test_div(self):
        a = Array([5],context=self.context,dtype=self.dtype)
        b = a / 5
        self.assertEqual(b[0],1)
        
    def test_rdiv(self):
        a = Array([5],context=self.context,dtype=self.dtype)
        b = 5 / a
        self.assertEqual(b[0],1)
        
    def test_idiv(self):
        a = Array([5],context=self.context,dtype=self.dtype)
        a /= 5
        self.assertEqual(a[0],1)       
    
    def test_sub(self):
        a = Array([5],context=self.context,dtype=self.dtype)
        b = a - 5
        self.assertEqual(b[0],0)
        
    def test_rsub(self):
        a = Array([5],context=self.context,dtype=self.dtype)
        b = 5 - a
        self.assertEqual(b[0],0)
        
    def test_isub(self):
        a = Array([5],context=self.context,dtype=self.dtype)
        a -= 5
        self.assertEqual(a[0],0)       
        
    def test_pow(self):
        a = Array([2],context=self.context,dtype=self.dtype)
        b = a ** 2
        self.assertEqual(b[0],4)
        
    def test_abs(self):
        a = Array([-5],context=self.context,dtype=self.dtype)
        b = abs(a)
        self.assertEqual(b[0],5)
        
    def test_real(self):
        a = Array([5],context=self.context,dtype=self.dtype)
        b = a.real()
        self.assertEqual(b[0],5)
        
    def test_imag(self):
        a = Array([5],context=self.context,dtype=self.dtype)
        b = a.imag()
        self.assertEqual(b[0],0)
        
    def test_conj(self):
        a = Array([5],context=self.context,dtype=self.dtype)
        b = a.conj()
        self.assertEqual(b[0],(5))
        
    

def array_test_maker(context,dtype):
    class tests(tests_base,unittest.TestCase):
        def __init__(self,*args):
            self.context=context
            self.dtype=dtype
            unittest.TestCase.__init__(self,*args)
    tests.__name__ = str(type(context)) + " " + dtype.__name__
    return tests

cpu_test = array_test_maker(CPUContext(),numpy.float64)
cpu_test2 = array_test_maker(CPUContext(),numpy.float32)
cpu_test3 = array_test_maker(CPUContext(),numpy.complex128)
cpu_test4 = array_test_maker(CPUContext(),numpy.complex64)

if pycbc.have_cuda:
    cuda_test = array_test_maker(CUDAContext(),numpy.float64)
    cuda_test2 = array_test_maker(CUDAContext(),numpy.float32)
    cuda_test3 = array_test_maker(CUDAContext(),numpy.complex128)
    cuda_test4 = array_test_maker(CUDAContext(),numpy.complex64)

if pycbc.have_opencl:
    cl_test = array_test_maker(OpenCLContext(),numpy.float64)
    cl_test2 = array_test_maker(OpenCLContext(),numpy.float32)
    cl_test3 = array_test_maker(OpenCLContext(),numpy.complex128)
    cl_test4 = array_test_maker(OpenCLContext(),numpy.complex64)


# TODO More specific array tests (instatiation, failure modes, type conversion, etc)
        
        
if __name__ == '__main__':
    unittest.main()
