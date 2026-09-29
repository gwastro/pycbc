# Copyright (C) 2026 Alex Nitz
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

"""Unit tests for gating and inpainting."""

import unittest
import numpy as np
from pycbc.types import TimeSeries
from pycbc.psd import welch, interpolate, inverse_spectrum_truncation
from pycbc.strain.gate import (gate_and_paint, gate_and_paint_matmul,
                               invert_covariance)


class TestGateAndPaint(unittest.TestCase):
    def setUp(self):
        np.random.seed(12345)
        self.sample_rate = 2048.0
        self.duration = 8.0
        self.N = int(self.duration * self.sample_rate)
        self.delta_t = 1.0 / self.sample_rate
        self.delta_f = 1.0 / self.duration

        white = np.random.randn(self.N)
        self.ts = TimeSeries(white, delta_t=self.delta_t, epoch=0)

        # Estimate PSD with Hann inverse spectrum truncation
        psd = welch(self.ts, seg_len=int(2 * self.sample_rate),
                    seg_stride=int(self.sample_rate))
        psd = interpolate(psd, self.delta_f)
        self.psd = inverse_spectrum_truncation(
            psd, max_filter_len=int(2 * self.sample_rate),
            low_frequency_cutoff=20.0, trunc_method='hann')
        self.invpsd = 1.0 / self.psd

        # Add a localized glitch at t = 4.0s
        t = self.ts.sample_times.numpy() - 4.0
        glitch = 50.0 * np.exp(- (t / 0.05)**2) * np.sin(2 * np.pi * 50.0 * t)
        self.ts_glitch = self.ts + TimeSeries(glitch, delta_t=self.delta_t, epoch=0)

    def test_inpaint_cholesky(self):
        """Test regularized Cholesky inpainting cancels overwhitened data."""
        lindex = int(3.9 * self.sample_rate)
        rindex = int(4.1 * self.sample_rate)  # 0.2s gate (410 samples)

        cleaned = gate_and_paint(self.ts_glitch, lindex, rindex, self.invpsd,
                                 copy=True, method='cholesky', ridge=1e-10)

        # Check that overwhitened data inside the gate is essentially zero
        owh_clean = (cleaned.to_frequencyseries() * self.invpsd).to_timeseries()
        owh_inside = owh_clean[lindex:rindex].numpy()
        diag = (self.invpsd.astype('complex').to_timeseries() * self.invpsd.delta_t)[0]

        # Residual should be tiny relative to diagonal
        self.assertLess(np.max(np.abs(owh_inside)) / diag, 1e-4)

    def test_inpaint_matmul(self):
        """Test matmul inpainting matches Cholesky."""
        lindex = int(3.98 * self.sample_rate)
        rindex = int(4.02 * self.sample_rate)

        invmat = invert_covariance(self.invpsd, lindex, rindex, ridge=1e-10)
        cleaned_mm = gate_and_paint_matmul(self.ts_glitch, lindex, rindex,
                                           self.invpsd, invmat=invmat, copy=True)
        cleaned_cho = gate_and_paint(self.ts_glitch, lindex, rindex, self.invpsd,
                                     copy=True, method='cholesky', ridge=1e-10)

        diff = np.max(np.abs(cleaned_mm.numpy() - cleaned_cho.numpy()))
        self.assertLess(diff, 1e-4)

    def test_compare_cholesky_to_toeplitz(self):
        """Test that Cholesky and Toeplitz solvers agree on well-behaved gate."""
        lindex = int(3.98 * self.sample_rate)
        rindex = int(4.02 * self.sample_rate)

        cleaned_cho = gate_and_paint(self.ts_glitch, lindex, rindex, self.invpsd,
                                     copy=True, method='cholesky', ridge=1e-10)
        cleaned_toep = gate_and_paint(self.ts_glitch, lindex, rindex, self.invpsd,
                                      copy=True, method='toeplitz')

        diff = np.max(np.abs(cleaned_cho.numpy() - cleaned_toep.numpy()))
        self.assertLess(diff, 1e-4)

    def test_numerical_stability_large_gate(self):
        """Test large gate where unregularized Levinson suffers precision loss."""
        lindex = int(3.8 * self.sample_rate)
        rindex = int(4.2 * self.sample_rate)

        cleaned_cho = gate_and_paint(self.ts_glitch, lindex, rindex, self.invpsd,
                                     copy=True, method='cholesky', ridge=1e-10)

        self.assertTrue(np.all(np.isfinite(cleaned_cho.numpy())))
        max_amp = np.max(np.abs(cleaned_cho[lindex:rindex].numpy()))
        self.assertLess(max_amp, 2000.0)

    def test_timeseries_gate_paint(self):
        """Test TimeSeries.gate with paint method and default Cholesky."""
        g_paint = self.ts_glitch.gate(4.0, window=0.1, method='paint',
                                      invpsd=self.invpsd, copy=True)
        self.assertIsNotNone(g_paint)
        self.assertTrue(np.all(np.isfinite(g_paint.numpy())))


if __name__ == '__main__':
    unittest.main()
