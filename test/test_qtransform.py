# Copyright (C) 2026 Alexander Harvey Nitz
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

"""
Unit tests for qtransform filter functions and bounds handling.
"""

import unittest

import numpy

from pycbc.filter.qtransform import qseries
from pycbc.types import TimeSeries


class TestQTransform(unittest.TestCase):
    def setUp(self):
        # 2 seconds at 1024 Hz -> 2048 samples
        self.data = TimeSeries(numpy.zeros(2048), delta_t=1.0 / 1024)
        self.data[1024] = 1.0  # Delta impulse at center
        self.fdata = self.data.to_frequencyseries()

    def test_qseries_standard(self):
        """Verify standard qseries execution and output length."""
        q = qseries(self.fdata, Q=8.0, f0=100.0)
        self.assertEqual(len(q), len(self.data))
        self.assertAlmostEqual(q.delta_t, self.data.delta_t)

    def test_qseries_bounds_low_frequency(self):
        """Verify qseries bounds checking when window start < 0."""
        # Low frequency with low Q pushes start < 0
        q = qseries(self.fdata, Q=2.0, f0=2.0)
        self.assertEqual(len(q), len(self.data))

    def test_qseries_bounds_high_frequency(self):
        """Verify qseries bounds checking when window end > Nyquist."""
        # High frequency near Nyquist (512 Hz) pushes end > len(fdata)
        q = qseries(self.fdata, Q=2.0, f0=500.0)
        self.assertEqual(len(q), len(self.data))

    def test_qseries_bounds_outside_nyquist(self):
        """Verify qseries clips cleanly when window lies outside data bound."""
        # Central frequency above Nyquist (512 Hz) clips to available range
        q = qseries(self.fdata, Q=8.0, f0=600.0)
        self.assertEqual(len(q), len(self.data))
        self.assertTrue(numpy.isfinite(q.numpy()).all())

    def test_qseries_return_complex(self):
        """Verify qseries return_complex returns complex TimeSeries."""
        q = qseries(self.fdata, Q=8.0, f0=100.0, return_complex=True)
        self.assertEqual(len(q), len(self.data))
        self.assertEqual(q.dtype, numpy.complex128)

    def test_qseries_odd_length_input(self):
        """Verify qseries executes cleanly on odd-length series."""
        for n in [1023, 1025, 2047, 2049]:
            ts = TimeSeries(numpy.zeros(n), delta_t=1.0 / 1024)
            ts[n // 2] = 1.0
            fs = ts.to_frequencyseries()
            expected_tlen = (len(fs) - 1) * 2

            q = qseries(fs, Q=8.0, f0=100.0)
            self.assertEqual(len(q), expected_tlen)
            self.assertTrue(numpy.isfinite(q.numpy()).all())

            qc = qseries(fs, Q=8.0, f0=100.0, return_complex=True)
            self.assertEqual(len(qc), expected_tlen)
            self.assertEqual(qc.dtype, numpy.complex128)
            self.assertTrue(numpy.isfinite(qc.numpy()).all())

    def test_qtransform_odd_length_series(self):
        """Verify TimeSeries.qtransform works on odd-length series."""
        for n in [1023, 1025, 2047, 2049]:
            ts = TimeSeries(numpy.zeros(n), delta_t=1.0 / 1024)
            ts[n // 2] = 1.0

            # Test interpolated qtransform
            times, freqs, qplane = ts.qtransform(
                delta_t=ts.delta_t, delta_f=1.0, frange=(30, 300)
            )
            self.assertEqual(len(times), n)
            self.assertEqual(qplane.shape, (len(freqs), n))
            self.assertTrue(numpy.isfinite(qplane).all())
            self.assertGreater(qplane.max(), 0.0)

            # Test uninterpolated qtransform
            times_raw, freqs_raw, qplane_raw = ts.qtransform(
                frange=(30, 300)
            )
            self.assertEqual(
                qplane_raw.shape, (len(freqs_raw), len(times_raw))
            )
            self.assertTrue(numpy.isfinite(qplane_raw).all())

            # Test complex output
            times_c, freqs_c, qplane_c = ts.qtransform(
                delta_t=ts.delta_t, delta_f=1.0, frange=(30, 300),
                return_complex=True
            )
            self.assertEqual(qplane_c.shape, (len(freqs_c), n))
            self.assertEqual(qplane_c.dtype, numpy.complex128)
            self.assertTrue(numpy.isfinite(qplane_c).all())


if __name__ == "__main__":
    unittest.main()
