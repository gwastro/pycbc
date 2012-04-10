import pycbc

if pycbc.HAVE_CUDA:
    import pycuda
    import pycuda.gpuarray

if pycbc.HAVE_OPENCL:
    import pyopencl
    import pyopencl.array

import unittest
from pycbc.array import *
from pycbc.timeseries import TimeSeries
from pycbc.scheme import *
import numpy 
import swiglal

class TestTimeSeriesBase(object):
    def setUp(self):
        # two timeseries without epoch
        self.ts1 = TimeSeries([1, 2, 3], 0.1, dtype=self.dtype)
        self.ts2 = TimeSeries([10, 20, 30], 0.1, dtype=self.dtype)
        # two timeseries with epoch
        self.epoch = swiglal.LIGOTimeGPS(1000, 1000)
        self.ts3 = TimeSeries([1, 2, 3], 0.1, epoch=self.epoch, dtype=self.odtype)
        self.ts4 = TimeSeries([10, 20, 30], 0.1, epoch=self.epoch, dtype=self.odtype)
        # timeseries with different delta_t
        self.ts5 = TimeSeries([10, 20, 30], 0.01, epoch=self.epoch, dtype=self.odtype)

    def test_duration(self):
        with self.context:
            self.assertEqual(self.ts1.duration, 0.1 * 3)
            self.assertEqual(self.ts2.duration, 0.1 * 3)
            self.assertEqual(self.ts3.duration, 0.1 * 3)
            self.assertEqual(self.ts4.duration, 0.1 * 3)
            self.assertEqual(self.ts5.duration, 0.01 * 3)

    def test_sum(self):
        with self.context:
            s = self.ts1 + self.ts2
            self.assertEqual(s[0], 11)
            self.assertEqual(s[1], 22)
            self.assertEqual(s[2], 33)
            self.assertEqual(s.delta_t, 0.1)
            s = self.ts3 + self.ts4
            self.assertEqual(s[0], 11)
            self.assertEqual(s[1], 22)
            self.assertEqual(s[2], 33)
            self.assertEqual(s.delta_t, 0.1)
            self.assertEqual(s.start_time, self.epoch)
            try:
                s = self.ts1 + self.ts3
                fail()
            except ValueError:
                pass
            try:
                s = self.ts3 + self.ts5
                fail()
            except ValueError:
                pass

def test_maker(context, dtype, odtype):
    class TestTimeSeries(TestTimeSeriesBase, unittest.TestCase):
        def __init__(self, *args):
            self.context = context
            self.dtype = dtype
            self.odtype = odtype
            unittest.TestCase.__init__(self, *args)
    TestTimeSeries.__name__ = str(type(context)) + " " + dtype.__name__ + " with " + odtype.__name__
    return TestTimeSeries

DefaultScheme._single = None
schemes = [DefaultScheme()]
types = [float32, float64, complex64, complex128]
if pycbc.HAVE_CUDA:
    DefaultScheme._single = None
    schemes.append(CUDAScheme())
if pycbc.HAVE_OPENCL:
    DefaultScheme._single = None
    schemes.append(OpenCLScheme())
tests = []

i = 0
for s in schemes:
    for t in types:
        for ot in types:
            na = 'test' + str(i)
            vars()[na] = test_maker(s, t, ot)
            i += 1

if __name__ == '__main__':
    unittest.main()
