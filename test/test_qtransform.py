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
from pycbc.types import TimeSeries
from pycbc.filter.qtransform import qseries


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

    def test_qseries_return_complex(self):
        """Verify qseries return_complex returns complex TimeSeries."""
        q = qseries(self.fdata, Q=8.0, f0=100.0, return_complex=True)
        self.assertEqual(len(q), len(self.data))
        self.assertEqual(q.dtype, numpy.complex128)


if __name__ == "__main__":
    unittest.main()
