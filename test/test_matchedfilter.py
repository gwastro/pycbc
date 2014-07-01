# Copyright (C) 2012  Alex Nitz, Josh Willis
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
"""
These are the unittests for the pycbc.filter.matchedfilter module
"""
import sys
import pycbc
import unittest
from pycbc.types import *
from pycbc.scheme import *
from pycbc.filter import *
from math import sqrt
import pycbc.fft
import numpy
from utils import parse_args_all_schemes, simple_exit

_scheme, _context = parse_args_all_schemes("Matched Filter")
import pycbc.fft.fftw
pycbc.fft.fftw.set_measure_level(0)

class TestMatchedFilter(unittest.TestCase):
    def setUp(self,*args):
        self.context = _context
        self.scheme = _scheme
        from math import sin
        # Use sine wave as test signal
        data = numpy.sin(numpy.arange(0,100,100/(4096.0*64)))
        self.filt = TimeSeries(data,dtype=float32,delta_t=1.0/4096)
        self.filt2 = (self.filt*1)
        self.filt2[0:len(self.filt2)/2].fill(0)
        self.filt_offset = TimeSeries(numpy.roll(data,4096*32), dtype=float32,
                                      delta_t=1.0/4096)

        self.filtD = TimeSeries(data,dtype=float64,delta_t=1.0/4096)
        self.filt2D = (self.filtD*1)
        self.filt2D[0:len(self.filt2D)/2].fill(0)
        self.filt_offsetD = TimeSeries(numpy.roll(data,4096*32), dtype=float64,
                                      delta_t=1.0/4096)

        self.filt_short =TimeSeries([0,1,2,3,4],dtype=float32,delta_t=1.0/4096)

    def test_correlate (self):
        from pycbc.filter.matchedfilter import correlate
        with self.context:
            a = Array([1j], dtype=complex64)
            b = Array([1j], dtype=complex64)
            c = zeros(1, dtype=complex64)
            correlate (a, b, c)
            self.assertEqual(1, c[0])

    def test_ave_snr_noise(self):
        with self.context:
            #Test that the average snr in noise is 2
            from numpy.random import normal

            noise = normal(0.0,2,4096*64)
            nplus= TimeSeries(noise,dtype=float32,delta_t=1.0/4096)
            ntilde = make_frequency_series(nplus) / nplus.delta_t
            # Calculate a Faux psd for normalization, replace with better algorithm
            psd = (ntilde).squared_norm()  / float(len(nplus)) * nplus.delta_t *2.0

            snr = matched_filter(self.filt, nplus, psd=psd)

            ave = snr.squared_norm().sum() /len(snr)
            self.assertAlmostEqual(2,ave,places=5)

            noise = normal(0.0,2,4096*64)
            nplus= TimeSeries(noise,dtype=float64,delta_t=1.0/4096)
            ntilde = make_frequency_series(nplus) / nplus.delta_t
            # Calculate a Faux psd for normalization, replace with better algorithm
            psd = (ntilde).squared_norm()  / float(len(nplus)) * nplus.delta_t *2.0

            snr = matched_filter(self.filtD,nplus,psd=psd)
            ave = snr.squared_norm().sum() /len(snr)
            self.assertAlmostEqual(2,ave,places=5)

    def test_perfect_match(self):
        with self.context:
            o,i = match(self.filt,self.filt)
            self.assertAlmostEqual(1,o,places=4)
            self.assertEqual(0,i)
            o,i = match(self.filtD,self.filtD)
            self.assertAlmostEqual(1,o,places=4)
            self.assertEqual(0,i)

    def test_perfect_match_offset(self):
        with self.context:
            o,i = match(self.filt,self.filt_offset)
            self.assertAlmostEqual(1,o,places=4)
            self.assertEqual(4096*32,i)

            o,i = match(self.filtD,self.filt_offsetD)
            self.assertAlmostEqual(1,o,places=4)
            self.assertEqual(4096*32,i)

    def test_imperfect_match(self):
        with self.context:
            f = make_frequency_series(self.filt)
            f2 = make_frequency_series(self.filt2)
            o,i = match(self.filt,self.filt2)
            self.assertAlmostEqual(sqrt(0.5),o,places=3)

            f = make_frequency_series(self.filtD)
            f2 = make_frequency_series(self.filt2D)
            o,i = match(self.filtD,self.filt2D)
            self.assertAlmostEqual(sqrt(0.5),o,places=3)

    def test_errors(self):
        with self.context:
            #Check that an incompatible data and filter produce an error
            self.assertRaises(ValueError,match,self.filt,self.filt[0:5])

            #Check that an incompatible psd produces an error
            self.assertRaises(TypeError,match,self.filt,self.filt,psd=self.filt)
            psd = FrequencySeries(zeros(len(self.filt)/2+1),delta_f=100000)
            self.assertRaises(TypeError,match,self.filt,self.filt,psd=psd)

            #Check that only TimeSeries or FrequencySeries are accepted
            self.assertRaises(TypeError,match,zeros(10),zeros(10))

            self.assertRaises(ValueError,match,self.filt,self.filt[0:len(self.filt)-1])


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMatchedFilter))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
