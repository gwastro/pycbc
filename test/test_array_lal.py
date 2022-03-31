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
import unittest
from pycbc.types import *
from pycbc.scheme import *
from lal import LIGOTimeGPS as LTG
from utils import parse_args_all_schemes, simple_exit

_scheme, _context = parse_args_all_schemes("lal() method")

class TestUtils(unittest.TestCase):
    def setUp(self,*args):
        self.context = _context
        self.delta_t = 1.0 / 4096
        self.epoch = LTG(0,0)

        self.at = TimeSeries([1], delta_t=self.delta_t, dtype=float32,epoch=self.epoch)
        self.bt = TimeSeries([1], delta_t=self.delta_t, dtype=float64,epoch=self.epoch)
        self.ct = TimeSeries([1], delta_t=self.delta_t, dtype=complex64,epoch=self.epoch)
        self.dt = TimeSeries([1], delta_t=self.delta_t, dtype=complex128,epoch=self.epoch)

        self.a = Array([1], dtype=float32)
        self.b = Array([1], dtype=float64)
        self.c = Array([1], dtype=complex64)
        self.d = Array([1], dtype=complex128)

        self.af = FrequencySeries([1], delta_f=self.delta_t, dtype=float32,epoch=self.epoch)
        self.bf = FrequencySeries([1], delta_f=self.delta_t, dtype=float64,epoch=self.epoch)
        self.cf = FrequencySeries([1], delta_f=self.delta_t, dtype=complex64,epoch=self.epoch)
        self.df = FrequencySeries([1], delta_f=self.delta_t, dtype=complex128,epoch=self.epoch)


    if _scheme == 'cpu':
        def test_array_to_lal(self):
            al = self.a.lal()
            self.assertEqual(al.data.dtype, self.a.dtype)
            self.assertEqual(al.data[0], self.a[0])
            al = self.b.lal()
            self.assertEqual(al.data.dtype, self.b.dtype)
            self.assertEqual(al.data[0], self.b[0])
            al = self.c.lal()
            self.assertEqual(al.data.dtype, self.c.dtype)
            self.assertEqual(al.data[0], self.c[0])
            al = self.d.lal()
            self.assertEqual(al.data.dtype, self.d.dtype)
            self.assertEqual(al.data[0], self.d[0])


        def test_timeseries_to_lal(self):
            al = self.at.lal()
            self.assertEqual(al.data.data.dtype, self.at.dtype)
            self.assertEqual(al.data.data[0], self.at[0])
            self.assertEqual(al.deltaT, self.at.delta_t)
            al = self.bt.lal()
            self.assertEqual(al.data.data.dtype, self.bt.dtype)
            self.assertEqual(al.data.data[0], self.bt[0])
            self.assertEqual(al.deltaT, self.bt.delta_t)
            al = self.ct.lal()
            self.assertEqual(al.data.data.dtype, self.ct.dtype)
            self.assertEqual(al.data.data[0], self.ct[0])
            self.assertEqual(al.deltaT, self.ct.delta_t)
            al = self.dt.lal()
            self.assertEqual(al.data.data.dtype, self.dt.dtype)
            self.assertEqual(al.data.data[0], self.dt[0])
            self.assertEqual(al.deltaT, self.dt.delta_t)

        def test_frequencyseries_to_lal(self):
            al = self.af.lal()
            self.assertEqual(al.data.data.dtype, self.af.dtype)
            self.assertEqual(al.data.data[0], self.af[0])
            self.assertEqual(al.deltaF, self.af.delta_f)
            al = self.bf.lal()
            self.assertEqual(al.data.data.dtype, self.bf.dtype)
            self.assertEqual(al.data.data[0], self.bf[0])
            self.assertEqual(al.deltaF, self.bf.delta_f)
            al = self.cf.lal()
            self.assertEqual(al.data.data.dtype, self.cf.dtype)
            self.assertEqual(al.data.data[0], self.cf[0])
            self.assertEqual(al.deltaF, self.cf.delta_f)
            al = self.df.lal()
            self.assertEqual(al.data.data.dtype, self.df.dtype)
            self.assertEqual(al.data.data[0], self.df[0])
            self.assertEqual(al.deltaF, self.df.delta_f)
    else:
        def test_array_lal_errors(self):
            with self.context:
                self.assertRaises(TypeError, self.a.lal)

suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestUtils))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
