# Copyright (C) 2012  Alex Nitz, Andrew Miller, Josh Willis, Tito Dal Canton
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
These are the unittests for the pycbc timeseries type
'''

import pycbc
import unittest
from pycbc.types import *
from pycbc.scheme import *
import numpy
import lal
from utils import array_base, parse_args_all_schemes, simple_exit
import sys
import os
import tempfile

_scheme, _context = parse_args_all_schemes("TimeSeries")

# By importing the current schemes array type, it will make it
# easier to check the  array types later
if _scheme == 'cuda':
    import pycuda
    import pycuda.gpuarray
    from pycuda.gpuarray import GPUArray as SchemeArray
elif _scheme == 'cpu':
    from numpy import ndarray as SchemeArray

from numpy import ndarray as CPUArray

class TestTimeSeriesBase(array_base,unittest.TestCase):
    __test__ = False
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

        # Note that self.epoch is set in the factory class constructor at the end;
        # we need only set self.delta_t here.
        self.delta_t = 0.1
        # We need to tell the arithmetic test functions what our type is:
        self.type = TimeSeries
        # and the extra keyword arguments the constructors will need:
        self.kwds = {'epoch': self.epoch, 'delta_t': self.delta_t}
        # Now that the kinds are set, we need to call our parent method to set up all the
        # inputs and answers for our functions
        self.setNumbers()

        # The above call created instances for all of our inputs and various correct
        # outputs.  But we make a copy of the scalar to check later.
        self.s = self.scalar

        # Finally, we want to have an array that we shouldn't be able to operate on,
        # because the precision is wrong, and one where the length is wrong.
        self.bad = TimeSeries([1,1,1], self.delta_t, epoch=self.epoch, dtype = self.other_precision[self.odtype])
        self.bad2 = TimeSeries([1,1,1,1], self.delta_t, epoch=self.epoch, dtype = self.dtype)

        # These are timeseries that have problems specific to timeseries
        self.bad3 = TimeSeries([1,1,1], 0.2, epoch=self.epoch, dtype = self.dtype)
        if self.epoch == 0:
            self.bad4 = TimeSeries([1,1,1], self.delta_t, epoch = lal.LIGOTimeGPS(1000, 1000), dtype = self.dtype)
        else:
            self.bad4 = TimeSeries([1,1,1], self.delta_t, epoch=None, dtype = self.dtype)

    def test_numpy_init(self):
        with self.context:
            in1 = numpy.array([5,3,1],dtype=self.odtype)
            in2 = numpy.array([5,3,1],dtype=self.other_precision[self.odtype])

            #We don't want to cast complex as real
            if not (self.kind == 'real' and self.okind == 'complex'):
                #First we must check that the dtype is correct when specified
                out1 = TimeSeries(in1,0.1, dtype=self.dtype, epoch=self.epoch)
                out2 = TimeSeries(in2,0.1, dtype=self.dtype, epoch=self.epoch)
                #to be sure that it is copied
                in1 += 1
                in2 += 1
                self.assertTrue(type(out1._scheme) == type(self.context))
                self.assertTrue(type(out1._data) is SchemeArray)
                self.assertEqual(out1[0],5)
                self.assertEqual(out1[1],3)
                self.assertEqual(out1[2],1)
                self.assertTrue(out1.dtype==self.dtype)
                self.assertEqual(out1.delta_t, 0.1)
                self.assertEqual(out1.start_time, self.epoch)


                self.assertTrue(type(out2._scheme) == type(self.context))
                self.assertTrue(type(out2._data) is SchemeArray)
                self.assertEqual(out2[0],5)
                self.assertEqual(out2[1],3)
                self.assertEqual(out2[2],1)
                self.assertTrue(out2.dtype==self.dtype)
                self.assertEqual(out2.delta_t,0.1)
                self.assertEqual(out2.start_time, self.epoch)

                in1-=1
                in2-=1

            # Also, when it is unspecified
            out3 = TimeSeries(in1,0.1,epoch=self.epoch)
            in1 += 1
            self.assertTrue(type(out3._scheme) == type(self.context))
            self.assertTrue(type(out3._data) is SchemeArray)
            self.assertEqual(out3[0],5)
            self.assertEqual(out3[1],3)
            self.assertEqual(out3[2],1)
            self.assertTrue(out3.dtype==self.odtype)
            self.assertEqual(out3.delta_t,0.1)
            self.assertEqual(out3.start_time, self.epoch)

            # Check for copy=false
            # On the CPU, this should be possible
            in3 = numpy.array([5,3,1],dtype=self.dtype)
            if self.scheme == 'cpu':
                out4 = TimeSeries(in3,0.1,epoch=self.epoch,copy=False)
                in3 += 1

                self.assertTrue(out4.dtype==self.dtype)
                self.assertTrue(type(out4._scheme) == type(self.context))
                self.assertEqual(out4[0],6)
                self.assertEqual(out4[1],4)
                self.assertEqual(out4[2],2)
                self.assertEqual(out4.delta_t,0.1)
                self.assertEqual(out4.start_time, self.epoch)

            # If we're in different scheme, this should raise an error
            else:
                self.assertRaises(TypeError, TimeSeries, in3, 0.1, copy=False)

            # We also need to check initialization using GPU arrays
            if self.scheme == 'cuda':
                in4 = pycuda.gpuarray.zeros(3,self.dtype)
            if self.scheme != 'cpu':
                out4 = TimeSeries(in4,0.1, copy=False, epoch=self.epoch)
                in4 += 1
                self.assertTrue(type(out4._scheme) == type(self.context))
                self.assertTrue(type(out4._data) is SchemeArray)
                self.assertEqual(out4[0],1)
                self.assertEqual(out4[1],1)
                self.assertEqual(out4[2],1)
                self.assertTrue(out4.dtype==self.dtype)
                self.assertEqual(out4.delta_t,0.1)
                self.assertEqual(out4.start_time, self.epoch)

            # We should be able to create an array from the wrong dtype, and
            # it should be cast as float64
            in5 = numpy.array([1,2,3],dtype=numpy.int32)
            out5 = TimeSeries(in5,0.1,epoch=self.epoch)
            in5 += 1
            self.assertTrue(type(out5._scheme) == type(self.context))
            self.assertTrue(type(out5._data) is SchemeArray)
            self.assertEqual(out5[0],1)
            self.assertEqual(out5[1],2)
            self.assertEqual(out5[2],3)
            #self.assertTrue(out5.dtype==numpy.float64)
            self.assertEqual(out5.delta_t,0.1)
            self.assertEqual(out5.start_time, self.epoch)

            # We shouldn't be able to copy it though
            #self.assertRaises(TypeError,TimeSeries,in5, 0.1, copy=False)

            # Finally, just checking a few things specific to timeseries
            inbad = numpy.array([],dtype=float64)
            self.assertRaises(ValueError, TimeSeries, in1, -1)
            self.assertRaises(ValueError, TimeSeries, inbad, .1)
            self.assertRaises(TypeError, TimeSeries, in1, .1, epoch=(5,1))

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
                out1 = TimeSeries(in1, 0.1, epoch=self.epoch, dtype=self.dtype)
                out2 = TimeSeries(in2, 0.1, epoch=self.epoch, dtype=self.dtype)
                # to be sure that it is copied
                in1 += 1
                in2 += 1

                self.assertTrue(type(out1._scheme) == type(self.context))
                self.assertTrue(type(out1._data) is SchemeArray)
                self.assertEqual(out1[0],5)
                self.assertEqual(out1[1],3)
                self.assertEqual(out1[2],1)
                self.assertTrue(out1.dtype==self.dtype)
                self.assertEqual(out1.delta_t, 0.1)
                self.assertEqual(out1.start_time, self.epoch)

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
                self.assertEqual(out2.delta_t, 0.1)
                self.assertEqual(out2.start_time, self.epoch)

                in1-=1
                in2-=1
            # Giving complex input and specifying a real dtype should raise an error
            else:
                self.assertRaises(TypeError, TimeSeries, in1,0.1, dtype = self.dtype)
                self.assertRaises(TypeError, TimeSeries, in2,0.1, dtype = self.dtype)

            # Also, when it is unspecified
            out3 = TimeSeries(in1,0.1,epoch=self.epoch)
            in1 += 1

            self.assertTrue(type(out3._scheme) == type(self.context))
            self.assertTrue(type(out3._data) is SchemeArray)
            self.assertEqual(out3[0],5)
            self.assertEqual(out3[1],3)
            self.assertEqual(out3[2],1)
            self.assertTrue(out3.dtype==self.odtype)
            self.assertEqual(out3.delta_t, 0.1)
            self.assertEqual(out3.start_time, self.epoch)

            # We should also be able to create from a CPU Array
            out4 = TimeSeries(cpuarray,0.1, dtype=self.dtype, epoch=self.epoch)

            self.assertTrue(type(out4._scheme) == type(self.context))
            self.assertTrue(type(out4._data) is SchemeArray)
            self.assertEqual(out4[0],1)
            self.assertEqual(out4[1],2)
            self.assertEqual(out4[2],3)
            self.assertTrue(out4.dtype==self.dtype)
            self.assertEqual(out4.delta_t, 0.1)
            self.assertEqual(out4.start_time, self.epoch)

            # Check for copy=false
            in3 = Array([5,3,1],dtype=self.dtype)
            out5 = TimeSeries(in3,0.1,copy=False,epoch=self.epoch)
            in3 += 1

            self.assertTrue(type(out5._scheme) == type(self.context))
            self.assertTrue(type(out5._data) is SchemeArray)
            self.assertEqual(out5[0],6)
            self.assertEqual(out5[1],4)
            self.assertEqual(out5[2],2)
            self.assertTrue(out5.dtype==self.dtype)
            self.assertEqual(out5.delta_t, 0.1)
            self.assertEqual(out5.start_time, self.epoch)

            if self.scheme != 'cpu':
                self.assertRaises(TypeError,TimeSeries,0.1,cpuarray,copy=False)
            # Things specific to timeseries
            inbad = Array(numpy.array([],dtype=float64))
            self.assertRaises(ValueError, TimeSeries, in1, -1)
            self.assertRaises(ValueError, TimeSeries, inbad, .1)
            self.assertRaises(TypeError, TimeSeries, in1, .1, epoch=(5,1))

        # Also checking that a cpu array can't be made out of another scheme without copying
        if self.scheme != 'cpu':
            self.assertRaises(TypeError, TimeSeries, out4, 0.1, copy=False, epoch=self.epoch)
            out6 = TimeSeries(out4, 0.1, dtype=self.dtype)
            self.assertTrue(type(out6._scheme) == CPUScheme)
            self.assertTrue(type(out6._data) is CPUArray)
            self.assertEqual(out6[0],1)
            self.assertEqual(out6[1],2)
            self.assertEqual(out6[2],3)
            self.assertTrue(out6.dtype==self.dtype)
            self.assertEqual(out6.delta_t, 0.1)
            self.assertEqual(out6.start_time, self.epoch)

    def test_list_init(self):
        with self.context:
            # When specified
            out1 = TimeSeries([5,3,1],0.1, dtype=self.dtype, epoch=self.epoch)

            self.assertTrue(type(out1._scheme) == type(self.context))
            self.assertTrue(type(out1._data) is SchemeArray)
            self.assertEqual(out1[0],5)
            self.assertEqual(out1[1],3)
            self.assertEqual(out1[2],1)
            self.assertTrue(out1.dtype==self.dtype)
            self.assertEqual(out1.delta_t, 0.1)
            self.assertEqual(out1.start_time, self.epoch)

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
                out2 = TimeSeries([5.0+0j,3+0j,1+0j], 0.1, dtype=self.dtype, epoch=self.epoch)

                self.assertTrue(type(out2._scheme) == type(self.context))
                self.assertTrue(type(out2._data) is SchemeArray)
                self.assertEqual(out2[0],5)
                self.assertEqual(out2[1],3)
                self.assertEqual(out2[2],1)
                self.assertTrue(out2.dtype==self.dtype)
                self.assertEqual(out2.delta_t, 0.1)
                self.assertEqual(out2.start_time, self.epoch)

            else:
                self.assertRaises(TypeError, TimeSeries,[5+0j, 3+0j, 1+0j], 0.1, dtype=self.dtype)

            #Also, when it is unspecified
            out3 = TimeSeries([5.0,3,1],0.1,epoch=self.epoch)

            self.assertTrue(type(out3._scheme) == type(self.context))
            self.assertTrue(type(out3._data) is SchemeArray)
            self.assertEqual(out3[0],5)
            self.assertEqual(out3[1],3)
            self.assertEqual(out3[2],1)
            self.assertTrue(out3.dtype==numpy.float64)
            self.assertEqual(out3.delta_t, 0.1)
            self.assertEqual(out3.start_time, self.epoch)

            out4 = TimeSeries([5+0j,3+0j,1+0j],0.1,epoch = self.epoch)

            self.assertTrue(type(out4._scheme) == type(self.context))
            self.assertTrue(type(out4._data) is SchemeArray)
            self.assertEqual(out4[0],5)
            self.assertEqual(out4[1],3)
            self.assertEqual(out4[2],1)
            self.assertTrue(out4.dtype==numpy.complex128)
            self.assertEqual(out4.delta_t, 0.1)
            self.assertEqual(out4.start_time, self.epoch)

            self.assertRaises(TypeError,TimeSeries,[1,2,3],copy=False)

            # Things specific to timeseries
            self.assertRaises(ValueError, TimeSeries, [1,2,3], -1)
            self.assertRaises(ValueError, TimeSeries, [], .1)
            self.assertRaises(TypeError, TimeSeries, [1,2,3], .1, epoch=(5,1))

    def test_mul(self):
        super(TestTimeSeriesBase,self).test_mul()
        self.assertRaises(ValueError, self.a.__mul__,self.bad3)
        self.assertRaises(ValueError, self.a.__mul__,self.bad4)

    def test_rmul(self):
        super(TestTimeSeriesBase,self).test_rmul()
        self.assertRaises(ValueError, self.a.__rmul__,self.bad3)
        self.assertRaises(ValueError, self.a.__rmul__,self.bad4)

    def test_imul(self):
        super(TestTimeSeriesBase,self).test_imul()
        self.assertRaises(ValueError, self.a.__imul__,self.bad3)
        self.assertRaises(ValueError, self.a.__imul__,self.bad4)

    def test_add(self):
        super(TestTimeSeriesBase,self).test_add()
        self.assertRaises(ValueError, self.a.__add__, self.bad3)
        self.assertRaises(ValueError, self.a.__add__, self.bad4)

    def test_radd(self):
        super(TestTimeSeriesBase,self).test_radd()
        self.assertRaises(ValueError, self.a.__radd__,self.bad3)
        self.assertRaises(ValueError, self.a.__radd__,self.bad4)

    def test_iadd(self):
        super(TestTimeSeriesBase,self).test_iadd()
        self.assertRaises(ValueError, self.a.__iadd__,self.bad3)
        self.assertRaises(ValueError, self.a.__iadd__,self.bad4)

    def test_sub(self):
        super(TestTimeSeriesBase,self).test_sub()
        self.assertRaises(ValueError, self.a.__sub__,self.bad3)
        self.assertRaises(ValueError, self.a.__sub__,self.bad4)

    def test_rsub(self):
        super(TestTimeSeriesBase,self).test_rsub()
        self.assertRaises(ValueError, self.a.__rsub__,self.bad3)
        self.assertRaises(ValueError, self.a.__rsub__,self.bad4)

    def test_isub(self):
        super(TestTimeSeriesBase,self).test_isub()
        self.assertRaises(ValueError, self.a.__isub__,self.bad3)
        self.assertRaises(ValueError, self.a.__isub__,self.bad4)

    def test_div(self):
        super(TestTimeSeriesBase,self).test_div()
        self.assertRaises(ValueError, self.a.__div__,self.bad3)
        self.assertRaises(ValueError, self.a.__div__,self.bad4)

    def test_rdiv(self):
        super(TestTimeSeriesBase,self).test_rdiv()
        self.assertRaises(ValueError, self.a.__rdiv__,self.bad3)
        self.assertRaises(ValueError, self.a.__rdiv__,self.bad4)

    def test_idiv(self):
        super(TestTimeSeriesBase,self).test_idiv()
        self.assertRaises(ValueError, self.a.__idiv__,self.bad3)
        self.assertRaises(ValueError, self.a.__idiv__,self.bad4)

    def test_dot(self):
        super(TestTimeSeriesBase,self).test_dot()
        self.assertRaises(ValueError, self.a.dot,self.bad3)
        self.assertRaises(ValueError, self.a.dot,self.bad4)

    def test_duration(self):
        with self.context:
            # Moving these to the current scheme
            self.a*=1
            self.b*=1
            self.bad3*=1
            self.assertAlmostEqual(self.a.duration, 0.3)
            self.assertAlmostEqual(self.b.duration, 0.3)
            self.assertAlmostEqual(self.bad3.duration, 0.6)

    def test_inject(self):
        a = TimeSeries(numpy.zeros(2**20, dtype=numpy.float32),
                                   delta_t=1.0)
        a[2**19] = 1

        # Check that the obvious case reduces to an add operation
        r = a.inject(a)
        self.assertAlmostEqual(r.max(), 2.0, places=7)

        # Check adding an offset vector
        b = a.cyclic_time_shift(-0.21)
        r = a.inject(b)
        self.assertAlmostEqual(r.max(), 2.0, places=5)

        # check adding shoter offset vector
        c = a.time_slice(2**19-5000, 2**19+5000).cyclic_time_shift(32.12)
        r = a.inject(c)
        self.assertAlmostEqual(r.max(), 2.0, places=4)

    def test_sample_times(self):
        with self.context:
            # Moving these to the current scheme
            self.a*=1
            self.b*=1
            self.bad3*=1
            self.assertEqual(len(self.a.sample_times), 3)
            self.assertAlmostEqual(self.a.sample_times[-1] - self.a.sample_times[0], 0.2)
            self.assertEqual(len(self.b.sample_times), 3)
            self.assertAlmostEqual(self.b.sample_times[-1] - self.b.sample_times[0], 0.2)
            self.assertEqual(len(self.bad3.sample_times), 3)
            self.assertAlmostEqual(self.bad3.sample_times[-1] - self.bad3.sample_times[0], 0.4)

    def test_save(self):
        with self.context:
            # make temporary file paths
            temp_file = tempfile.NamedTemporaryFile()
            temp_path_npy = temp_file.name + '.npy'
            temp_path_txt = temp_file.name + '.txt'
            # make a test time series
            a_numpy = numpy.arange(100, dtype=self.dtype)
            a = TimeSeries(a_numpy, delta_t=0.1)
            # test saving to Numpy array
            a.save(temp_path_npy)
            b = numpy.load(temp_path_npy)
            self.assertEqual(b.shape, (a_numpy.shape[0], 2))
            self.assertEqual(numpy.abs(b[:,0] - a.sample_times.numpy()).max(), 0)
            self.assertEqual(numpy.abs(b[:,1] - a_numpy).max(), 0)
            os.remove(temp_path_npy)
            # test saving to text file
            a.save(temp_path_txt)
            b = numpy.loadtxt(temp_path_txt)
            if a.kind == 'complex':
                self.assertEqual(b.shape, (a_numpy.shape[0], 3))
                b = numpy.vstack((b[:,0], b[:,1] + 1j * b[:,2])).T
            elif a.kind == 'real':
                self.assertEqual(b.shape, (a_numpy.shape[0], 2))
            self.assertEqual(numpy.abs(b[:,0] - a.sample_times.numpy()).max(), 0)
            self.assertEqual(numpy.abs(b[:,1] - a_numpy).max(), 0)
            os.remove(temp_path_txt)

def ts_test_maker(dtype, odtype, epoch):
    class TestTimeSeries(TestTimeSeriesBase):
        __test__ = True
        def __init__(self, *args):
            self.dtype = dtype
            self.odtype = odtype
            self.epoch = epoch if epoch is not None else lal.LIGOTimeGPS(0, 0)
            unittest.TestCase.__init__(self, *args)
    TestTimeSeries.__name__ = _scheme + " " + dtype.__name__ + " with " + odtype.__name__
    return TestTimeSeries

types = [ (float32,[float32,complex64]), (float64,[float64,complex128]),
        (complex64,[complex64,float32]), (complex128,[float64,complex128]) ]

suite = unittest.TestSuite()

# Unlike the regular array tests, we will need to test with an epoch, and with none
epochs = [lal.LIGOTimeGPS(1000, 1000), None]

i = 0
for t,otypes in types:
    for ot in otypes:
        for epoch in epochs:
            na = 'test' + str(i)
            vars()[na] = ts_test_maker(t, ot, epoch)
            suite.addTest(unittest.TestLoader().loadTestsFromTestCase(vars()[na]))
            i += 1

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
