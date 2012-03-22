import pycbc

if pycbc.HAVE_CUDA:
    import pycuda
    import pycuda.gpuarray

if pycbc.HAVE_OPENCL:
    import pyopencl
    import pyopencl.array

import unittest
from pycbc.array import *
from pycbc.scheme import *
import numpy 

# ********************GENERIC ARRAY TESTS ***********************

class tests_base(object):
    def setUp(self):
        self.a = Array([5], dtype = self.dtype)
        self.v = Array([10],dtype = self.odtype)
        self.n = 5
        pass

    def test_init(self):
        with self.context:      
            b = numpy.array([5])
            c = Array(b  ,dtype=self.dtype)
            self.assertEqual(self.a[0],b[0])
            self.assertEqual(c[0],b[0])
        
    def test_mul(self):
        with self.context:   
            b = self.a * 2
            c = self.a * self.v
            self.assertEqual(c[0],50)
            self.assertEqual(b[0],10)
        
    def test_rmul(self):
        with self.context:       
            b = 2 * self.a 
            c = self.v * self.a
            self.assertEqual(c[0],50) 
            self.assertEqual(b[0],10)
        
    def test_imul(self):
        with self.context:         
            t = self.a * 1 
            t2 = self.a * 1
            t *= 2
            t2 *= self.v
            self.assertEqual(t2[0],50) 
            self.assertEqual(t[0],10)
        
    def test_add(self):
        with self.context:         
            b = self.a + 5
            c = self.a + self.v
            self.assertEqual(c[0],15)
            self.assertEqual(b[0],10)
        
    def test_radd(self):
        with self.context:        
            b = 5 + self.a
            c = self.v + self.a
            self.assertEqual(c[0],15)
            self.assertEqual(b[0],10)
        
    def test_iadd(self):
        with self.context:
            t = self.a * 1
            t2 = self.a * 1         
            t += 5
            t2 += self.v
            self.assertEqual(t2[0],15)
            self.assertEqual(t[0],10)

    def test_div(self):
        with self.context:         
            b = self.a / 5
            c = self.a /self.v
            self.assertEqual(c[0],1.0/2)
            self.assertEqual(b[0],1)
        
    def test_rdiv(self):
        with self.context:        
            b = 5 / self.a
            c = self.v /self.a
            self.assertEqual(c[0],2)
            self.assertEqual(b[0],1)
        
    def test_idiv(self):
        with self.context:  
            t = self.a * 1
            t2 = self.a * 1       
            t /= 5
            t2 /= self.v
            self.assertEqual(t2[0],1./2)
            self.assertEqual(t[0],1)       
    
    def test_sub(self):
        with self.context:       
            b = self.a - 5
            c = self.a - self.v
            self.assertEqual(c[0],-5)
            self.assertEqual(b[0],0)
        
    def test_rsub(self):
        with self.context:         
            b = 5 - self.a
            c = self.v - self.a
            self.assertEqual(c[0],5)
            self.assertEqual(b[0],0)
        
    def test_isub(self):
        with self.context:
            t = self.a * 1 
            t2 = self.a * 1        
            t -= 5
            t2 -= self.v
            self.assertEqual(t2[0],-5)
            self.assertEqual(t[0],0)       
        
    def test_pow(self):
        with self.context:
            b = self.a ** 2
            n = abs(b[0]-25)
           # c = 2 ** self.v
           # self.assertTrue(c[0],1024)
            self.assertTrue(n<1e-5)
        
    def test_abs(self):
        with self.context:
            b = abs(self.a)
            self.assertEqual(b[0],5)
        
    def test_real(self):
        with self.context:        
            b = self.a.real()
            self.assertEqual(b[0],5)
        
    def test_imag(self):
        with self.context:        
            b = self.a.imag()
            self.assertEqual(b[0],0)
        
    def test_conj(self):
        with self.context:       
            b = self.a.conj()
            self.assertEqual(b[0],(5))
            
    def test_sum(self):
        with self.context:         
            b = self.a.sum()
            self.assertEqual(b,(5))
            
    def test_dot(self):
        with self.context:        
            b = self.a.dot(self.a)
            self.assertEqual(b,(25))
        
    

def array_test_maker(context,dtype,odtype):
    class tests(tests_base,unittest.TestCase):
        def __init__(self,*args):
            self.context=context
            self.dtype=dtype
            self.odtype=odtype
            unittest.TestCase.__init__(self,*args)
    tests.__name__ = str(type(context)) + " " + dtype.__name__ + " with " + odtype.__name__
#    DefaultScheme._single = None
    return tests

scs = [DefaultScheme()]
types = [float32,float64,complex64,complex128]
if pycbc.HAVE_CUDA:
    DefaultScheme._single = None
    scs.append(CUDAScheme())
if pycbc.HAVE_OPENCL:
    DefaultScheme._single = None
    scs.append(OpenCLScheme())
tests = []

ind = 0
for sc in scs:
    for ty in types:
        for ot in types:
            na = 'test' + str(ind)
            vars()[na] = array_test_maker(sc,ty,ot)
            ind += 1



# TODO More specific array tests (instatiation, failure modes, type conversion, etc)
        
        
if __name__ == '__main__':
    unittest.main()
