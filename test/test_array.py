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
        self.a = Array([5], dtype = self.dtype)
        self.v = Array([10],dtype = self.odtype)
        self.d = Array([1,2,3],dtype = self.dtype)
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
    
    def test_max(self):
        with self.context:
            self.assertEqual(self.d.max(),3)
            
    def test_min(self):
        with self.context:
            self.assertEqual(self.d.min(),1)
                
    

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

tests = []

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
