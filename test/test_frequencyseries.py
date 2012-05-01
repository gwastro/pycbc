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
import swiglal

class TestFrequencySeriesBase(object):
    def setUp(self):
        self.fs1 = FrequencySeries([1, 2, 3], 0.1, dtype=self.dtype)
        self.fs2 = FrequencySeries([10, 20, 30], 0.2, dtype=self.dtype)
        self.epoch = swiglal.LIGOTimeGPS(1000, 1000)
        self.fs3 = FrequencySeries([1, 2, 3], 0.1, epoch=self.epoch, dtype=self.odtype)

    def test_sample_frequencies(self):
        with self.context:
            self.assertEqual(len(self.fs1.sample_frequencies), 3)
            self.assertEqual(self.fs1.sample_frequencies[-1], 0.2)
            self.assertEqual(len(self.fs2.sample_frequencies), 3)
            self.assertEqual(self.fs2.sample_frequencies[-1], 0.4)
            self.assertEqual(len(self.fs3.sample_frequencies), 3)
            self.assertEqual(self.fs3.sample_frequencies[-1], 0.2)

    def test_sum(self):
        with self.context:
            s = self.fs1 + self.fs3
            self.assertEqual(s[0], 2)
            self.assertEqual(s[1], 4)
            self.assertEqual(s[2], 6)
            self.assertEqual(s.delta_f, 0.1)
            try:
                s = self.fs1 + self.fs2
                fail()
            except ValueError:
                pass

def test_maker(context, dtype, odtype):
    class TestFrequencySeries(TestFrequencySeriesBase, unittest.TestCase):
        def __init__(self, *args):
            self.context = context
            self.dtype = dtype
            self.odtype = odtype
            unittest.TestCase.__init__(self, *args)
    TestFrequencySeries.__name__ = str(type(context)) + " " + dtype.__name__ + " with " + odtype.__name__
    return TestFrequencySeries

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
