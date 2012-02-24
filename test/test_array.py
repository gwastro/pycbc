import pycbc

if pycbc.HAVE_CUDA:
    import pycuda
    import pycuda.gpuarray

if pycbc.HAVE_OPENCL:
    import pyopencl
    import pyopencl.array

import unittest
from pycbc.array import Array
from pycbc.scheme import *
import numpy 

# ********************GENERIC ARRAY TESTS ***********************

class tests_base(object):

    def setUp(self):
        pass

    def test_init(self):
        with self.context:
            a = Array([5]  ,dtype=self.dtype)
            b = numpy.array([5])
            c = Array(b  ,dtype=self.dtype)
            self.assertEqual(a[0],b[0])
            self.assertEqual(c[0],b[0])
        
    def test_mul(self):
        with self.context:
            a = Array([5]  ,dtype=self.dtype)
            b = a * 2
            self.assertEqual(b[0],10)
        
    def test_rmul(self):
        with self.context:
            a = Array([5]  ,dtype=self.dtype)
            b = 2 * a
            self.assertEqual(b[0],10)
        
    def test_imul(self):
        with self.context:
            a = Array([5]  ,dtype=self.dtype)
            a *= 2
            self.assertEqual(a[0],10)
        
    def test_add(self):
        with self.context:
            a = Array([5]  ,dtype=self.dtype)
            b = a + 5
            self.assertEqual(b[0],10)
        
    def test_radd(self):
        with self.context:
            a = Array([5]  ,dtype=self.dtype)
            b = 5 + a
            self.assertEqual(b[0],10)
        
    def test_iadd(self):
        with self.context:
            a = Array([5]  ,dtype=self.dtype)
            a += 5
            self.assertEqual(a[0],10)

    def test_div(self):
        with self.context:
            a = Array([5]  ,dtype=self.dtype)
            b = a / 5
            self.assertEqual(b[0],1)
        
    def test_rdiv(self):
        with self.context:
            a = Array([5]  ,dtype=self.dtype)
            b = 5 / a
            self.assertEqual(b[0],1)
        
    def test_idiv(self):
        with self.context:
            a = Array([5]  ,dtype=self.dtype)
            a /= 5
            self.assertEqual(a[0],1)       
    
    def test_sub(self):
        with self.context:
            a = Array([5]  ,dtype=self.dtype)
            b = a - 5
            self.assertEqual(b[0],0)
        
    def test_rsub(self):
        with self.context:
            a = Array([5]  ,dtype=self.dtype)
            b = 5 - a
            self.assertEqual(b[0],0)
        
    def test_isub(self):
        with self.context:
            a = Array([5]  ,dtype=self.dtype)
            a -= 5
            self.assertEqual(a[0],0)       
        
    def test_pow(self):
        with self.context:
            a = Array([2]  ,dtype=self.dtype)
            b = a ** 2
            self.assertEqual(b[0],4)
        
    def test_abs(self):
        with self.context:
            a = Array([-5]  ,dtype=self.dtype)
            b = abs(a)
            self.assertEqual(b[0],5)
        
    def test_real(self):
        with self.context:
            a = Array([5]  ,dtype=self.dtype)
            b = a.real()
            self.assertEqual(b[0],5)
        
    def test_imag(self):
        with self.context:
            a = Array([5]  ,dtype=self.dtype)
            b = a.imag()
            self.assertEqual(b[0],0)
        
    def test_conj(self):
        with self.context:
            a = Array([5]  ,dtype=self.dtype)
            b = a.conj()
            self.assertEqual(b[0],(5))
            
    def test_sum(self):
        with self.context:
            a = Array([5]  ,dtype=self.dtype)
            b = a.sum()
            self.assertEqual(b,(5))
            
    def test_dot(self):
        with self.context:
            a = Array([5]  ,dtype=self.dtype)
            b = a.dot(a)
            self.assertEqual(b,(25))
        
    

def array_test_maker(context,dtype):
    class tests(tests_base,unittest.TestCase):
        def __init__(self,*args):
            self.context=context
            self.dtype=dtype
            unittest.TestCase.__init__(self,*args)
    tests.__name__ = str(type(context)) + " " + dtype.__name__
    return tests

cpu_test = array_test_maker(CPUScheme(),numpy.float64)
cpu_test2 = array_test_maker(CPUScheme(),numpy.float32)
cpu_test3 = array_test_maker(CPUScheme(),numpy.complex128)
cpu_test4 = array_test_maker(CPUScheme(),numpy.complex64)

if pycbc.HAVE_CUDA:
    cuda_test = array_test_maker(CUDAScheme(),numpy.float64)
    cuda_test2 = array_test_maker(CUDAScheme(),numpy.float32)
    cuda_test3 = array_test_maker(CUDAScheme(),numpy.complex128)
    cuda_test4 = array_test_maker(CUDAScheme(),numpy.complex64)

if pycbc.HAVE_OPENCL:
    cl_test = array_test_maker(OpenCLScheme(),numpy.float64)
    cl_test2 = array_test_maker(OpenCLScheme(),numpy.float32)
    cl_test3 = array_test_maker(OpenCLScheme(),numpy.complex128)
    cl_test4 = array_test_maker(OpenCLScheme(),numpy.complex64)


# TODO More specific array tests (instatiation, failure modes, type conversion, etc)
        
        
if __name__ == '__main__':
    unittest.main()
