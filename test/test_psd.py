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
from pycbc.psd import generate_segment_psds
from pycbc.strain import StrainSegments
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

    def _white_noise(self, size, sample_freq=1.):
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

    def _psd_window(self, noise, pdl, seg_start, seg_stop, ana_start, ana_stop):
        """(start, stop) chosen by generate_segment_psds for one segment.
        psd_num_segments=1 makes the estimation stretch exactly pdl samples.
        """
        opt = SimpleNamespace(
            psd_estimation='median', psd_segment_length=pdl,
            psd_segment_stride=pdl, psd_num_segments=1,
            psd_inverse_length=None, psd_model=None, psd_file=None,
            asd_file=None, psd_low_frequency_cutoff=None,
            invpsd_trunc_method=None)
        with self.context:
            out = generate_segment_psds(
                opt, noise, [(seg_start, seg_stop, ana_start, ana_stop)],
                pdl // 2 + 1, 1. / pdl, self.psd_low_freq_cutoff)
        start, stop, _ = out[0]
        return start, stop

    def test_window_shorter_than_analysed(self):
        # PSD stretch shorter than the analysed span -> centred in it
        noise = self._white_noise(5000)
        start, stop = self._psd_window(noise, 300, 1000, 2000, 1200, 1700)
        self.assertEqual((start, stop), (1300, 1600))
        self.assertEqual((start + stop) // 2, (1200 + 1700) // 2)

    def test_window_equal_to_analysed(self):
        noise = self._white_noise(5000)
        start, stop = self._psd_window(noise, 500, 1000, 2000, 1200, 1700)
        self.assertEqual((start, stop), (1200, 1700))

    def test_window_between_analysed_and_segment(self):
        noise = self._white_noise(5000)
        # covers the analysed span; spare length goes before it first
        start, stop = self._psd_window(noise, 700, 1000, 2000, 1200, 1700)
        self.assertEqual((start, stop), (1000, 1700))
        # once the "before" side is exhausted the rest spills after
        start, stop = self._psd_window(noise, 850, 1000, 2000, 1200, 1700)
        self.assertEqual((start, stop), (1000, 1850))

    def test_window_equal_to_segment(self):
        noise = self._white_noise(5000)
        start, stop = self._psd_window(noise, 1000, 1000, 2000, 1200, 1700)
        self.assertEqual((start, stop), (1000, 2000))

    def test_window_longer_than_segment(self):
        # centred on the segment
        noise = self._white_noise(5000)
        start, stop = self._psd_window(noise, 1400, 1000, 2000, 1200, 1700)
        self.assertEqual((start, stop), (800, 2200))
        self.assertEqual((start + stop) // 2, (1000 + 2000) // 2)

    def test_window_slides_inside_data_at_edges(self):
        noise = self._white_noise(5000)
        # near the start: cannot begin before sample 0
        start, stop = self._psd_window(noise, 1400, 0, 1000, 100, 600)
        self.assertEqual((start, stop), (0, 1400))
        # near the end: cannot finish past the last sample
        start, stop = self._psd_window(noise, 1400, 4000, 5000, 4200, 4700)
        self.assertEqual((start, stop), (3600, 5000))

    def test_window_with_zero_padding(self):
        noise = self._white_noise(5000)
        # zero-padded leading segment: seg_start < 0, analysed span reaches
        # into the pad; the stretch must stay within the real data [0, n]
        for pdl in (300, 600, 700, 900):
            start, stop = self._psd_window(noise, pdl, -300, 700, -100, 400)
            self.assertGreaterEqual(start, 0)
            self.assertLessEqual(stop, 5000)
            self.assertEqual(stop - start, pdl)
        # zero-padded trailing segment
        for pdl in (300, 600, 700, 900):
            start, stop = self._psd_window(noise, pdl, 4300, 5300, 4600, 5200)
            self.assertGreaterEqual(start, 0)
            self.assertLessEqual(stop, 5000)
            self.assertEqual(stop - start, pdl)

    def test_data_length_clamped_and_warns(self):
        # --psd-num-segments omitted and data too short for the stride -> the
        # segment count is clamped to 1, so the estimation stretch is never
        # shorter than one Welch segment, and a warning is emitted. (Patch
        # logging.warning rather than use assertLogs: another test module
        # calls logging.disable().)
        noise = self._white_noise(4096)
        opt = SimpleNamespace(
            psd_estimation='median', psd_segment_length=1024,
            psd_segment_stride=4096, psd_num_segments=None,
            psd_inverse_length=None, psd_model=None, psd_file=None,
            asd_file=None, psd_low_frequency_cutoff=None,
            invpsd_trunc_method=None)
        with self.context, \
                mock.patch.object(pycbc.psd.logging, 'warning') as warn:
            out = generate_segment_psds(opt, noise, [(0, 4096, 0, 4096)],
                                        513, 1. / 1024,
                                        self.psd_low_freq_cutoff)
        self.assertTrue(
            any('Welch segment' in call.args[0] for call in warn.call_args_list))
        self.assertEqual(out[0][1] - out[0][0], 1024)
        # an explicit (even degenerate) value is also clamped, without that
        # warning (segment fully covered by the window, so no other warning
        # fires either)
        opt.psd_num_segments = 0
        with self.context, \
                mock.patch.object(pycbc.psd.logging, 'warning') as warn:
            out = generate_segment_psds(opt, noise, [(0, 1024, 0, 1024)],
                                        513, 1. / 1024,
                                        self.psd_low_freq_cutoff)
        self.assertFalse(warn.called)
        self.assertEqual(out[0][1] - out[0][0], 1024)

    def test_generate_segment_psds_alignment(self):
        # multiple real segments: in-bounds windows, repeats share the PSD
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
        for start, stop, _ in out:
            self.assertGreaterEqual(start, 0)
            self.assertLessEqual(stop, n)
            self.assertEqual(stop - start, 3 * 2048 + 4096)
        self.assertIs(out[0][2], out[3][2])

    def test_associate_psds_to_segments(self):
        # end-to-end: associate_psds_to_segments reads seg_slice/analyze
        # straight off the fourier segments and assigns every one a PSD
        noise = self._white_noise(8192)
        opt = SimpleNamespace(
            psd_estimation='median', psd_segment_length=256,
            psd_segment_stride=128, psd_num_segments=4,
            psd_inverse_length=None, psd_model=None, psd_file=None,
            asd_file=None, psd_low_frequency_cutoff=None,
            invpsd_trunc_method=None, segment_length=4096,
            segment_start_pad=128, segment_end_pad=16,
            trig_start_time=0, trig_end_time=0, filter_inj_only=False,
            injection_window=None, allow_zero_padding=False)
        strain_segments = StrainSegments.from_cli(opt, noise)
        fd_segments = strain_segments.fourier_segments()
        with self.context:
            pycbc.psd.associate_psds_to_segments(
                opt, fd_segments, noise, strain_segments.freq_len,
                strain_segments.delta_f, self.psd_low_freq_cutoff)
        self.assertTrue(fd_segments)
        self.assertTrue(all(fs.psd is not None for fs in fd_segments))


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestPSD))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestPSDSegmentPlacement))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
