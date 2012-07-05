import pycbc
import unittest
from pycbc.types import *
from pycbc.scheme import *
import numpy 
import swiglal
import array_base

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

class TestTimeSeriesBase(array_base.array_base):

    def checkScheme(self, a, b, s, c, c_ans):
    
            self.assertTrue(type(a._scheme) == self.scheme)
            self.assertTrue(type(a._data) is SchemeArray)
            self.assertEqual(a[0],self.alist[0])
            self.assertEqual(a[1],self.alist[1])
            self.assertEqual(a[2],self.alist[2])
            self.assertEqual(a.delta_t, self.alist[4])
            self.assertEqual(a.start_time, self.epoch)
            
            self.assertTrue(type(b._scheme) == self.scheme)
            self.assertTrue(type(b._data) is SchemeArray)
            self.assertEqual(b[0],self.blist[0])
            self.assertEqual(b[1],self.blist[1])
            self.assertEqual(b[2],self.blist[2])
            self.assertEqual(b.delta_t, self.blist[4])
            self.assertEqual(b.start_time, self.epoch)
            
            self.assertEqual(s, self.s2)
            
            if type(c_ans) == list:
                if c_ans[3]:
                    self.assertTrue(type(c._scheme) == self.scheme)
                    self.assertTrue(type(c._data) is SchemeArray)
                    self.assertEqual(c[0], c_ans[0])
                    self.assertEqual(c[1], c_ans[1])
                    self.assertEqual(c[2], c_ans[2])
                    self.assertEqual(c.delta_t, self.alist[4])
                    self.assertEqual(c.start_time, self.epoch)
                    
                else:
                    self.assertTrue(type(c._scheme) == self.scheme)
                    self.assertTrue(type(c._data) is SchemeArray)
                    self.assertAlmostEqual(c[0], c_ans[0], self.places)
                    self.assertAlmostEqual(c[1], c_ans[1], self.places)
                    self.assertAlmostEqual(c[2], c_ans[2], self.places)
                    self.assertEqual(c.delta_t, self.alist[4])
                    self.assertEqual(c.start_time, self.epoch)
                    
            else:
                self.assertEqual(c, c_ans)
                
    def checkCPU(self, a, b, s, c, c_ans):
            self.assertTrue(a._scheme == None)
            self.assertTrue(type(a._data) is numpy.ndarray)
            self.assertEqual(a[0],self.alist[0])
            self.assertEqual(a[1],self.alist[1])
            self.assertEqual(a[2],self.alist[2])
            self.assertEqual(a.delta_t, self.alist[4])
            self.assertEqual(a.start_time, self.epoch)
            
            self.assertTrue(b._scheme == None)
            self.assertTrue(type(b._data) is numpy.ndarray)
            self.assertEqual(b[0],self.blist[0])
            self.assertEqual(b[1],self.blist[1])
            self.assertEqual(b[2],self.blist[2])
            self.assertEqual(b.delta_t, self.blist[4])
            self.assertEqual(b.start_time, self.epoch)            
            
            self.assertEqual(s, self.s2)
            
            if type(c_ans) == list:
                if c_ans[3]:
                    self.assertTrue(c._scheme == None)
                    self.assertTrue(type(c._data) is numpy.ndarray)
                    self.assertEqual(c[0], c_ans[0])
                    self.assertEqual(c[1], c_ans[1])
                    self.assertEqual(c[2], c_ans[2])
                    self.assertEqual(c.delta_t, self.alist[4])
                    self.assertEqual(c.start_time, self.epoch)
                    
                else:
                    self.assertTrue(c._scheme == None)
                    self.assertTrue(type(c._data) is numpy.ndarray)
                    self.assertAlmostEqual(c[0], c_ans[0], self.places)
                    self.assertAlmostEqual(c[1], c_ans[1], self.places)
                    self.assertAlmostEqual(c[2], c_ans[2], self.places)
                    self.assertEqual(c.delta_t, self.alist[4])
                    self.assertEqual(c.start_time, self.epoch)
                    
            else:
                self.assertEqual(c, c_ans)
                
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
            
        # Here two test timeseries are created. Their content depends on their dtype,
        # so that we can make use of non-zero imaginary parts.
        # These arrays are just arbitrarily chosen, but in such a way that they are not
        # oversimplified, e.g. we don't want to multiply by one, or add zero, or have
        # the answer be one or zero.
        # We will make multiples of each to check things are on the cpu/gpu when they should be
        # and a list containing the values so we can tell it hasn't altered anything it shouldn't.
        
        if self.kind == 'real':
            self.a1 = TimeSeries([5,3,1], 0.1, epoch=self.epoch, dtype=self.dtype)
            self.a2 = TimeSeries([5,3,1], 0.1, epoch=self.epoch, dtype=self.dtype)
            self.a3 = TimeSeries([5,3,1], 0.1, epoch=self.epoch, dtype=self.dtype)
            self.alist = [5,3,1, True, 0.1]
        else:
            self.a1 = TimeSeries([5+1j,3+3j,1+5j], 0.1, epoch=self.epoch, dtype=self.dtype)
            self.a2 = TimeSeries([5+1j,3+3j,1+5j], 0.1, epoch=self.epoch, dtype=self.dtype)
            self.a3 = TimeSeries([5+1j,3+3j,1+5j], 0.1, epoch=self.epoch, dtype=self.dtype)
            self.alist = [5+1j,3+3j,1+5j, True, 0.1]
            
        if self.okind == 'real':
            self.b1 = TimeSeries([10,8,6], 0.1, epoch=self.epoch, dtype=self.odtype)
            self.b2 = TimeSeries([10,8,6], 0.1, epoch=self.epoch, dtype=self.odtype)
            self.blist = [10,8,6, True, 0.1]
        else:
            self.b1 = TimeSeries([10+6j,8+4j,6+2j], 0.1, epoch=self.epoch, dtype=self.odtype)
            self.b2 = TimeSeries([10+6j,8+4j,6+2j], 0.1, epoch=self.epoch, dtype=self.odtype)
            self.blist = [10+6j,8+4j,6+2j, True, 0.1]
        
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
        # because the precision is wrong, and one where the length is wrong.
        self.bad = TimeSeries([1,1,1], 0.1, epoch=self.epoch, dtype = self.other_precision[self.odtype])
        self.bad2 = TimeSeries([1,1,1,1], 0.1, epoch=self.epoch, dtype = self.dtype)
        
        # These are timeseries that have problems specific to timeseries
        self.bad3 = TimeSeries([1,1,1], 0.2, epoch=self.epoch, dtype = self.dtype)
        if self.epoch is None:
            self.bad4 = TimeSeries([1,1,1], 0.1, epoch = swiglal.LIGOTimeGPS(1000, 1000), dtype = self.dtype)
        else:
            self.bad4 = TimeSeries([1,1,1], 0.1, epoch=None, dtype = self.dtype)
        
        # We now need to set all the answer arrays up
        self.setAnswers()
        
    def test_mul(self):
        super(TestTimeSeriesBase,self).test_mul()
        self.assertRaises(ValueError, self.a1.__mul__,self.bad3)
        self.assertRaises(ValueError, self.a1.__mul__,self.bad4)
        
    def test_rmul(self):
        super(TestTimeSeriesBase,self).test_rmul()
        self.assertRaises(ValueError, self.a1.__rmul__,self.bad3)
        self.assertRaises(ValueError, self.a1.__rmul__,self.bad4)
        
    def test_imul(self):
        super(TestTimeSeriesBase,self).test_imul()
        self.assertRaises(ValueError, self.a1.__imul__,self.bad3)
        self.assertRaises(ValueError, self.a1.__imul__,self.bad4)
        
    def test_add(self):
        super(TestTimeSeriesBase,self).test_add()
        self.assertRaises(ValueError, self.a1.__add__,self.bad3)
        self.assertRaises(ValueError, self.a1.__add__,self.bad4)
        
    def test_radd(self):
        super(TestTimeSeriesBase,self).test_radd()
        self.assertRaises(ValueError, self.a1.__radd__,self.bad3)
        self.assertRaises(ValueError, self.a1.__radd__,self.bad4)
        
    def test_iadd(self):
        super(TestTimeSeriesBase,self).test_iadd()
        self.assertRaises(ValueError, self.a1.__iadd__,self.bad3)
        self.assertRaises(ValueError, self.a1.__iadd__,self.bad4)
        
    def test_sub(self):
        super(TestTimeSeriesBase,self).test_sub()
        self.assertRaises(ValueError, self.a1.__sub__,self.bad3)
        self.assertRaises(ValueError, self.a1.__sub__,self.bad4)
        
    def test_rsub(self):
        super(TestTimeSeriesBase,self).test_rsub()
        self.assertRaises(ValueError, self.a1.__rsub__,self.bad3)
        self.assertRaises(ValueError, self.a1.__rsub__,self.bad4)
        
    def test_isub(self):
        super(TestTimeSeriesBase,self).test_isub()
        self.assertRaises(ValueError, self.a1.__isub__,self.bad3)
        self.assertRaises(ValueError, self.a1.__isub__,self.bad4)
        
    def test_div(self):
        super(TestTimeSeriesBase,self).test_div()
        self.assertRaises(ValueError, self.a1.__div__,self.bad3)
        self.assertRaises(ValueError, self.a1.__div__,self.bad4)
        
    def test_rdiv(self):
        super(TestTimeSeriesBase,self).test_rdiv()
        self.assertRaises(ValueError, self.a1.__rdiv__,self.bad3)
        self.assertRaises(ValueError, self.a1.__rdiv__,self.bad4)
        
    def test_idiv(self):
        super(TestTimeSeriesBase,self).test_idiv()
        self.assertRaises(ValueError, self.a1.__idiv__,self.bad3)
        self.assertRaises(ValueError, self.a1.__idiv__,self.bad4)
        
    def test_duration(self):
        with self.context:
            # Moving these to the current scheme
            self.a1*=1
            self.b1*=1
            self.bad3*=1
            self.assertAlmostEqual(self.a1.duration, 0.3)
            self.assertAlmostEqual(self.b1.duration, 0.3)
            self.assertAlmostEqual(self.bad3.duration, 0.6)

    def test_sample_times(self):
        with self.context:
            # Moving these to the current scheme
            self.a1*=1
            self.b1*=1
            self.bad3*=1
            self.assertEqual(len(self.a1.sample_times), 3)
            self.assertAlmostEqual(self.a1.sample_times[-1] - self.a1.sample_times[0], 0.2)
            self.assertEqual(len(self.b1.sample_times), 3)
            self.assertAlmostEqual(self.b1.sample_times[-1] - self.b1.sample_times[0], 0.2)
            self.assertEqual(len(self.bad3.sample_times), 3)
            self.assertAlmostEqual(self.bad3.sample_times[-1] - self.bad3.sample_times[0], 0.4)

def test_maker(context, dtype, odtype, epoch):
    class TestTimeSeries(TestTimeSeriesBase, unittest.TestCase):
        def __init__(self, *args):
            self.context = context
            self.dtype = dtype
            self.odtype = odtype
            if _options['scheme'] == 'cpu':
                self.scheme = type(None)
            elif _options['scheme'] == 'cuda':
                self.scheme = pycbc.scheme.CUDAScheme
            else:
                self.scheme = pycbc.scheme.OpenCLScheme
            self.epoch = epoch
            unittest.TestCase.__init__(self, *args)
    TestTimeSeries.__name__ = _options['scheme'] + " " + dtype.__name__ + " with " + odtype.__name__
    return TestTimeSeries

types = [ (float32,[float32,complex64]), (float64,[float64,complex128]),
        (complex64,[complex64,float32]), (complex128,[float64,complex128]) ]

suite = unittest.TestSuite()

# Unlike the regular array tests, we will need to test with an epoch, and with none
epochs = [swiglal.LIGOTimeGPS(1000, 1000),None]

schemes = []

if _options['scheme'] == 'cpu':
    schemes.append(DefaultScheme())
if _options['scheme'] == 'cuda':
    schemes.append(CUDAScheme(device_num=_options['devicenum']))
if _options['scheme'] == 'opencl':
    schemes.append(OpenCLScheme(device_num=_options['devicenum']))

i = 0
for s in schemes:
    for t,otypes in types:
        for ot in otypes:
            for epoch in epochs:
                na = 'test' + str(i)
                vars()[na] = test_maker(s, t, ot, epoch)
                suite.addTest(unittest.TestLoader().loadTestsFromTestCase(vars()[na]))
                i += 1

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite)
