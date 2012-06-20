import pycbc
import unittest
from pycbc.types import *
from pycbc.scheme import *
import numpy 
import swiglal

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
    TestFrequencySeries.__name__ = _options['scheme'] + " " + dtype.__name__ + " with " + odtype.__name__
    return TestFrequencySeries

types = [ (float32,[float32,complex64]), (float64,[float64,complex128]),
        (complex64,[complex64,float32]), (complex128,[float64,complex128]) ]

schemes = []

if _options['scheme']=='cpu':
    schemes.append(DefaultScheme())
    
if _options['scheme']=='cuda':
    schemes.append(CUDAScheme(device_num=_options['devicenum']))

if _options['scheme']=='opencl':
    schemes.append(OpenCLScheme(device_num=_options['devicenum']))

tests = []

suite = unittest.TestSuite()

i = 0
for s in schemes:
    for t,otypes in types:
        for ot in otypes:
            na = 'test' + str(i)
            vars()[na] = test_maker(s, t, ot)
            suite.addTest(unittest.TestLoader().loadTestsFromTestCase(vars()[na]))
            i += 1

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite)
