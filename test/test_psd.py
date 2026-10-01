# Copyright (C) 2012  Tito Dal Canton, Josh Willis
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
These are the unittests for the pycbc PSD module.
'''

import os
import tempfile
from types import SimpleNamespace
from unittest import mock
import pycbc
import pycbc.psd
from pycbc.psd import (psd_estimation_window, psd_estimation_data_length,
                       generate_segment_psds)
from pycbc.types import TimeSeries, FrequencySeries
from pycbc.fft import ifft
from pycbc.fft.fftw import set_measure_level
import unittest
import numpy
from utils import parse_args_all_schemes, simple_exit
set_measure_level(0)

_scheme, _context = parse_args_all_schemes("PSD")

class TestPSD(unittest.TestCase):
    def setUp(self):
        self.scheme = _scheme
        self.context = _context
        self.psd_len = 1024
        self.psd_delta_f = 0.1
        self.psd_low_freq_cutoff = 10.
        # generate 1/f noise for testing PSD estimation
        noise_size = 524288
        sample_freq = 4096.
        delta_f = sample_freq / noise_size
        numpy.random.seed(132435)
        fd_size = noise_size // 2 + 1
        noise = numpy.random.normal(loc=0, scale=1, size=fd_size) + \
            1j * numpy.random.normal(loc=0, scale=1, size=fd_size)
        noise_model = 1. / numpy.linspace(1., 100., fd_size)
        noise *= noise_model / numpy.sqrt(delta_f) / 2
        noise[0] = noise[0].real
        noise_fs = FrequencySeries(noise, delta_f=delta_f)
        self.noise = TimeSeries(numpy.zeros(noise_size), delta_t=1./sample_freq)
        ifft(noise_fs, self.noise)

    def test_analytical(self):
        """Basic test of lalsimulation's analytical noise PSDs"""
        with self.context:
            psd_list = pycbc.psd.analytical.get_lalsim_psd_list()
            self.assertTrue(psd_list)
            for psd_name in psd_list:
                psd = pycbc.psd.analytical.from_string(psd_name, self.psd_len,
                                    self.psd_delta_f, self.psd_low_freq_cutoff)
                psd_min = psd.min()
                self.assertTrue(psd_min >= 0,
                                          msg=(psd_name + ': negative values'))
                self.assertTrue(psd.min() < 1e-40,
                                msg=(psd_name + ': unreasonably high minimum'))

    def test_read(self):
        """Test reading PSDs from text files"""
        test_data = numpy.zeros((self.psd_len, 2))
        test_data[:, 0] = numpy.linspace(0.,
                           (self.psd_len - 1) * self.psd_delta_f, self.psd_len)
        test_data[:, 1] = numpy.sqrt(test_data[:, 0])
        file_desc, file_name = tempfile.mkstemp()
        os.close(file_desc)
        numpy.savetxt(file_name, test_data)
        test_data[test_data[:, 0] < self.psd_low_freq_cutoff, 1] = 0.
        with self.context:
            psd = pycbc.psd.read.from_txt(file_name, self.psd_len,
                                    self.psd_delta_f, self.psd_low_freq_cutoff, is_asd_file=True)
            self.assertAlmostEqual(abs(psd - test_data[:, 1] ** 2).max(), 0)
        os.unlink(file_name)

    def test_estimate_welch(self):
        """Test estimating PSDs from data using Welch's method"""
        for seg_len in (2048, 4096, 8192):
            noise_model = (numpy.linspace(1., 100., seg_len//2 + 1)) ** (-2)
            for seg_stride in (seg_len, seg_len//2):
                for method in ('mean', 'median', 'median-mean'):
                    with self.context:
                        psd = pycbc.psd.welch(self.noise, seg_len=seg_len, \
                            seg_stride=seg_stride, avg_method=method)
                        error = (psd.numpy() - noise_model) / noise_model
                    err_rms = numpy.sqrt(numpy.mean(error ** 2))
                    self.assertTrue(err_rms < 0.2,
                        msg='seg_len=%d seg_stride=%d method=%s -> rms=%.3f' % \
                        (seg_len, seg_stride, method, err_rms))

    def test_truncation(self):
        """Test inverse PSD truncation"""
        for seg_len in (2048, 4096, 8192):
            noise_model = (numpy.linspace(1., 100., seg_len//2 + 1)) ** (-2)
            for max_len in (1024, 512, 256):
                with self.context:
                    psd = pycbc.psd.welch(self.noise, seg_len=seg_len, \
                                          seg_stride=seg_len//2, avg_method='mean')
                    psd_trunc = pycbc.psd.inverse_spectrum_truncation(
                            psd, max_len,
                            low_frequency_cutoff=self.psd_low_freq_cutoff)
                    freq = psd.sample_frequencies.numpy()
                    error = (psd.numpy() - noise_model) / noise_model
                error = error[freq > self.psd_low_freq_cutoff]
                err_rms = numpy.sqrt(numpy.mean(error ** 2))
                self.assertTrue(err_rms < 0.1,
                                msg='seg_len=%d max_len=%d -> rms=%.3f' \
                                % (seg_len, max_len, err_rms))

class TestPSDSegmentPlacement(unittest.TestCase):
    """Tests for associating an estimated PSD with each analysis segment."""

    def setUp(self):
        self.scheme = _scheme
        self.context = _context
        self.psd_low_freq_cutoff = 10.

    def _white_noise(self, size, sample_freq):
        numpy.random.seed(132435)
        fd_size = size // 2 + 1
        delta_f = sample_freq / size
        noise = numpy.random.normal(size=fd_size) \
            + 1j * numpy.random.normal(size=fd_size)
        noise /= numpy.sqrt(delta_f) * 2
        noise[0] = noise[0].real
        out = TimeSeries(numpy.zeros(size), delta_t=1. / sample_freq)
        ifft(FrequencySeries(noise, delta_f=delta_f), out)
        return out

    def test_window_shorter_than_analysed(self):
        # PSD stretch shorter than the analysed span -> centred in it
        start, stop = psd_estimation_window(300, 1000, 2000, 1200, 1700, 5000)
        self.assertEqual((start, stop), (1300, 1600))
        self.assertEqual((start + stop) // 2, (1200 + 1700) // 2)

    def test_window_equal_to_analysed(self):
        start, stop = psd_estimation_window(500, 1000, 2000, 1200, 1700, 5000)
        self.assertEqual((start, stop), (1200, 1700))

    def test_window_between_analysed_and_segment(self):
        # covers the analysed span; spare length goes before it first
        start, stop = psd_estimation_window(700, 1000, 2000, 1200, 1700, 5000)
        self.assertEqual((start, stop), (1000, 1700))
        # once the "before" side is exhausted the rest spills after
        start, stop = psd_estimation_window(850, 1000, 2000, 1200, 1700, 5000)
        self.assertEqual((start, stop), (1000, 1850))

    def test_window_equal_to_segment(self):
        start, stop = psd_estimation_window(1000, 1000, 2000, 1200, 1700, 5000)
        self.assertEqual((start, stop), (1000, 2000))

    def test_window_longer_than_segment(self):
        # centred on the segment
        start, stop = psd_estimation_window(1400, 1000, 2000, 1200, 1700, 5000)
        self.assertEqual((start, stop), (800, 2200))
        self.assertEqual((start + stop) // 2, (1000 + 2000) // 2)

    def test_window_slides_inside_data_at_edges(self):
        # near the start: cannot begin before sample 0
        start, stop = psd_estimation_window(1400, 0, 1000, 100, 600, 5000)
        self.assertEqual((start, stop), (0, 1400))
        # near the end: cannot finish past the last sample
        start, stop = psd_estimation_window(1400, 4000, 5000, 4200, 4700, 5000)
        self.assertEqual((start, stop), (3600, 5000))

    def test_window_with_zero_padding(self):
        # zero-padded leading segment: seg_start < 0, analysed span reaches
        # into the pad; the stretch must stay within the real data [0, n]
        for pdl in (300, 600, 700, 900):
            start, stop = psd_estimation_window(pdl, -300, 700, -100, 400, 5000)
            self.assertGreaterEqual(start, 0)
            self.assertLessEqual(stop, 5000)
            self.assertEqual(stop - start, pdl)
        # zero-padded trailing segment
        for pdl in (300, 600, 700, 900):
            start, stop = psd_estimation_window(pdl, 4300, 5300, 4600, 5200,
                                                5000)
            self.assertGreaterEqual(start, 0)
            self.assertLessEqual(stop, 5000)
            self.assertEqual(stop - start, pdl)

    def test_data_length_clamped_and_warns(self):
        # --psd-num-segments omitted, data too short for the stride -> the
        # segment count is clamped to 1 so the length is never < one Welch
        # segment, and a warning is emitted.
        opt = SimpleNamespace(psd_segment_length=4, psd_segment_stride=2,
                              psd_num_segments=None)
        with mock.patch.object(pycbc.psd.logging, 'warning') as warn:
            self.assertEqual(psd_estimation_data_length(opt, 1, 3), 4)
        self.assertTrue(warn.called)
        # explicit small/zero values are also clamped to at least one segment,
        # without a warning
        opt = SimpleNamespace(psd_segment_length=4, psd_segment_stride=2,
                              psd_num_segments=0)
        with mock.patch.object(pycbc.psd.logging, 'warning') as warn:
            self.assertEqual(psd_estimation_data_length(opt, 1, 1000), 4)
        self.assertFalse(warn.called)

    def test_segment_bounds(self):
        from pycbc.psd import _segment_bounds
        # from the psd_seg_bounds set by StrainSegments.fourier_segments
        seg = SimpleNamespace(psd_seg_bounds=(10, 20, 12, 18))
        self.assertEqual(_segment_bounds(seg), (10, 20, 12, 18))
        # backward compatible: derive from seg_slice/analyze if that is all the
        # segment carries (segments built the pre-psd_seg_bounds way)
        seg = SimpleNamespace(seg_slice=slice(10, 20), analyze=slice(2, 8))
        self.assertEqual(_segment_bounds(seg), (10, 20, 12, 18))

    def test_generate_segment_psds_alignment(self):
        sr = 4096.
        noise = self._white_noise(262144, sr)   # 64 s
        n = len(noise)
        seg_len = 8192
        ana = (1024, seg_len - 1024)
        analysis_segments = [
            (0, seg_len, ana[0], ana[1]),
            (seg_len, 2 * seg_len, ana[0], ana[1]),
            (n - seg_len, n, ana[0], ana[1]),
            (0, seg_len, ana[0], ana[1]),   # repeat -> must reuse the PSD
        ]
        opt = SimpleNamespace(
            psd_estimation='median',
            psd_segment_length=4096 / sr, psd_segment_stride=2048 / sr,
            psd_num_segments=4, psd_inverse_length=None,
            psd_model=None, psd_file=None, asd_file=None,
            psd_low_frequency_cutoff=None, invpsd_trunc_method=None)
        flen = seg_len // 2 + 1
        delta_f = sr / seg_len
        with self.context:
            out = generate_segment_psds(opt, noise, analysis_segments,
                                        flen, delta_f, self.psd_low_freq_cutoff)
        self.assertEqual(len(out), len(analysis_segments))
        pdl = psd_estimation_data_length(opt, sr, n)
        for (s0, s1, a0, a1), (start, stop, _) in zip(analysis_segments, out):
            self.assertEqual((start, stop),
                             psd_estimation_window(pdl, s0, s1, a0, a1, n))
            self.assertGreaterEqual(start, 0)
            self.assertLessEqual(stop, n)
        # segments resolving to the same window share the PSD object
        self.assertIs(out[0][2], out[3][2])


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestPSD))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestPSDSegmentPlacement))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
