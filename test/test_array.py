
# Check for optional components
have_cuda=False
have_opencl=False

try:
    import pycuda
    import pycuda.gpuarray
    import pycuda.autoinit
    have_cuda=True
except ImportError:
    have_cuda=False
    
try:
    import pyopencl
    import pyopencl.array
    have_opencl=True
except ImportError:
    have_opencl=False


import unittest
from pycbc.array import Array
from pycbc.processingcontext import *
import numpy 

class cpu(unittest.TestCase):

    def setUp(self):
        pass

    def test_init(self):
        a = Array([5],context=CPUContext())
        b = numpy.array([5])
        c = Array(b,context=CPUContext())
        self.assertEqual(a._data[0],b[0])
        self.assertEqual(c._data[0],b[0])
        
    def test_mul(self):
        a = Array([5],context=CPUContext())
        b = a * 2
        self.assertEqual(b._data[0],10)
        
    def test_rmul(self):
        a = Array([5],context=CPUContext())
        b = 2 * a
        self.assertEqual(b._data[0],10)
        
    def test_imul(self):
        a = Array([5],context=CPUContext())
        a *= 2
        self.assertEqual(a._data[0],10)
        
    def test_add(self):
        a = Array([5],context=CPUContext())
        b = a + 5
        self.assertEqual(b._data[0],10)
        
    def test_radd(self):
        a = Array([5],context=CPUContext())
        b = 5 + a
        self.assertEqual(b._data[0],10)
        
    def test_iadd(self):
        a = Array([5],context=CPUContext())
        a += 5
        self.assertEqual(a._data[0],10)

    def test_div(self):
        a = Array([5],context=CPUContext())
        b = a / 5
        self.assertEqual(b._data[0],1)
        
    def test_rdiv(self):
        a = Array([5],context=CPUContext())
        b = 5 / a
        self.assertEqual(b._data[0],1)
        
    def test_idiv(self):
        a = Array([5],context=CPUContext())
        a /= 5
        self.assertEqual(a._data[0],1)       
    
    def test_sub(self):
        a = Array([5],context=CPUContext())
        b = a - 5
        self.assertEqual(b._data[0],0)
        
    def test_rsub(self):
        a = Array([5],context=CPUContext())
        b = 5 - a
        self.assertEqual(b._data[0],0)
        
    def test_isub(self):
        a = Array([5],context=CPUContext())
        a -= 5
        self.assertEqual(a._data[0],0)       
        
    def test_pow(self):
        a = Array([5],context=CPUContext())
        b = a ** 2
        self.assertEqual(b._data[0],25)
        
    def test_abs(self):
        a = Array([-5],context=CPUContext())
        b = abs(a)
        self.assertEqual(b._data[0],5)
        
    def test_real(self):
        a = Array([1+1j],context=CPUContext())
        b = a.real()
        self.assertEqual(b._data[0],1)
        
    def test_imag(self):
        a = Array([1+2j],context=CPUContext())
        b = a.imag()
        self.assertEqual(b._data[0],2)
        
    def test_conj(self):
        a = Array([1+2j],context=CPUContext())
        b = a.conj()
        self.assertEqual(b._data[0],(1-2j))
        
    

if have_opencl: 

    class OpenCL(unittest.TestCase):

        def setUp(self):
            pass

        def test_init(self):
            a = Array([5],context=OpenCLContext())
            b = numpy.array([5])
            c = Array(b,context=OpenCLContext())
            a_cpu = a._data.get()
            c_cpu = c._data.get()
            self.assertEqual(a_cpu[0],b[0])
            self.assertEqual(c_cpu[0],b[0])
            
        def test_mul(self):
            a = Array([5],context=OpenCLContext())
            b = a * 2
            b_cpu = b._data.get()
            self.assertEqual(b_cpu[0],10)
            
        def test_rmul(self):
            a = Array([5],context=OpenCLContext())
            b = 2 * a
            b_cpu = b._data.get()
            self.assertEqual(b_cpu[0],10)
            
        def test_imul(self):
            a = Array([5],context=OpenCLContext())
            a *= 2
            a_cpu = a._data.get()
            self.assertEqual(a_cpu[0],10)
            
        def test_add(self):
            a = Array([5],context=OpenCLContext())
            b = a + 5
            b_cpu = b._data.get()
            self.assertEqual(b_cpu[0],10)
            
        def test_radd(self):
            a = Array([5],context=OpenCLContext())
            b = 5 + a
            b_cpu = b._data.get()
            self.assertEqual(b_cpu[0],10)
            
        def test_iadd(self):
            a = Array([5],context=OpenCLContext())
            a += 5
            a_cpu = a._data.get()
            self.assertEqual(a_cpu[0],10)

        def test_div(self):
            a = Array([5.0],context=OpenCLContext())
            b = a / 5.0
            b_cpu = b._data.get()
            self.assertEqual(b_cpu[0],1)
            
        def test_rdiv(self):
            a = Array([5.0],context=OpenCLContext())
            b = 5.0 / a
            b_cpu = b._data.get()
            self.assertEqual(b_cpu[0],1)
            
        def test_idiv(self):
            a = Array([5.0],context=OpenCLContext())
            a /= 5.0
            a_cpu = a._data.get()
            self.assertEqual(a_cpu[0],1)       
        
        def test_sub(self):
            a = Array([5],context=OpenCLContext())
            b = a - 5
            b_cpu = b._data.get()
            self.assertEqual(b_cpu[0],0)
            
        def test_rsub(self):
            a = Array([5],context=OpenCLContext())
            b = 5 - a
            b_cpu = b._data.get()
            self.assertEqual(b_cpu[0],0)
            
        def test_isub(self):
            a = Array([5],context=OpenCLContext())
            a -= 5
            a_cpu = a._data.get()
            self.assertEqual(a_cpu[0],0)       
            
        def test_pow(self):
            a = Array([5.0],context=OpenCLContext())
            b = a ** 2.0
            b_cpu = b._data.get()
            self.assertEqual(b_cpu[0],25)
            
        def test_abs(self):
            a = Array([-5],context=OpenCLContext())
            b = abs(a)
            b_cpu = b._data.get()
            self.assertEqual(b_cpu[0],5)
            
        def test_real(self):
            a = Array([1+1j],context=OpenCLContext())
            b = a.real()
            b_cpu = b._data.get()
            self.assertEqual(b_cpu[0],1)
            
        def test_imag(self):
            a = Array([1+2j],context=OpenCLContext())
            b = a.imag()
            b_cpu = b._data.get()
            self.assertEqual(b_cpu[0],2)
            
        def test_conj(self):
            a = Array([1+2j],context=OpenCLContext())
            b = a.conj()
            b_cpu = b._data.get()
            self.assertEqual(b_cpu[0],(1-2j))
       
if have_cuda:
 
    class CUDA(unittest.TestCase):


        def setUp(self):
            pass

        def test_init(self):
            a = Array([5],context=CUDAContext())
            b = numpy.array([5])
            c = Array(b,context=CUDAContext())
            a_cpu = a._data.get()
            c_cpu = c._data.get()
            self.assertEqual(a_cpu[0],b[0])
            self.assertEqual(c_cpu[0],b[0])
            
        def test_mul(self):
            a = Array([5],context=CUDAContext())
            b = a * 2
            b_cpu = b._data.get()
            self.assertEqual(b_cpu[0],10)
            
        def test_rmul(self):
            a = Array([5],context=CUDAContext())
            b = 2 * a
            b_cpu = b._data.get()
            self.assertEqual(b_cpu[0],10)
            
        def test_imul(self):
            a = Array([5],context=CUDAContext())
            a *= 2
            a_cpu = a._data.get()
            self.assertEqual(a_cpu[0],10)
            
        def test_add(self):
            a = Array([5],context=CUDAContext())
            b = a + 5
            b_cpu = b._data.get()
            self.assertEqual(b_cpu[0],10)
            
        def test_radd(self):
            a = Array([5],context=CUDAContext())
            b = 5 + a
            b_cpu = b._data.get()
            self.assertEqual(b_cpu[0],10)
            
        def test_iadd(self):
            a = Array([5],context=CUDAContext())
            a += 5
            a_cpu = a._data.get()
            self.assertEqual(a_cpu[0],10)

        def test_div(self):
            a = Array([5.0],context=CUDAContext())
            b = a / 5.0
            b_cpu = b._data.get()
            self.assertEqual(b_cpu[0],1)
            
        def test_rdiv(self):
            a = Array([5.0],context=CUDAContext(),dtype=numpy.float32)
            b = 5.0 / a
            b_cpu = b._data.get()
            self.assertEqual(b_cpu[0],1)
            
        def test_idiv(self):
            a = Array([5.0],context=CUDAContext())
            a /= 5.0
            a_cpu = a._data.get()
            self.assertEqual(a_cpu[0],1)       
        
        def test_sub(self):
            a = Array([5],context=CUDAContext())
            b = a - 5
            b_cpu = b._data.get()
            self.assertEqual(b_cpu[0],0)
            
        def test_rsub(self):
            a = Array([5],context=CUDAContext())
            b = 5 - a
            b_cpu = b._data.get()
            self.assertEqual(b_cpu[0],0)
            
        def test_isub(self):
            a = Array([5],context=CUDAContext())
            a -= 5
            a_cpu = a._data.get()
            self.assertEqual(a_cpu[0],0)       
            
        def test_pow(self):
            a = Array([5.0],context=CUDAContext())
            b = a ** 2.0
            b_cpu = b._data.get()
            self.assertEqual(b_cpu[0],25)
            
        def test_abs(self):
            a = Array([-5],context=CUDAContext())
            b = abs(a)
            b_cpu = b._data.get()
            self.assertEqual(b_cpu[0],5)
            
        def test_real(self):
            a = Array([1+1j],context=CUDAContext())
            b = a.real()
            b_cpu = b._data.get()
            self.assertEqual(b_cpu[0],1)
            
        def test_imag(self):
            a = Array([1+2j],context=CUDAContext())
            b = a.imag()
            b_cpu = b._data.get()
            self.assertEqual(b_cpu[0],2)
            
        def test_conj(self):
            a = Array([1+2j],context=CUDAContext())
            b = a.conj()
            b_cpu = b._data.get()
            self.assertEqual(b_cpu[0],(1-2j))
        
if __name__ == '__main__':
    unittest.main()
