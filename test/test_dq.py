# Copyright (C) 2019  Alex Nitz
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
These are the unittests for the data quality query methods in pycbc
"""
import unittest
from utils import simple_exit
from pycbc.dq import query_flag, query_cumulative_flags, query_str


class TestDataQualityFlags(unittest.TestCase):
    def setUp(self,*args):
        pass

    def test_direct_empty_return(self):
        segs = query_flag('H1', 'DATA',  1126051217,  1126051217 + 1000, cache=True)
        self.assertTrue(len(segs) == 0)

    def test_direct_query(self):
        segs = query_flag('H1', 'DATA', 1126051217,  1126051217 + 100000, cache=True)
        self.assertTrue(len(segs) > 0)

    def test_veto_flag(self):
        d = query_flag('L1', 'data', 1126051217,  1126051217 + 100000, cache=True)

        v1n = query_flag('L1', 'CBC_HW_INJ', 1126051217,  1126051217 + 100000, cache=True)
        v1 = query_flag('L1', 'CBC_HW_INJ', 1126051217,  1126051217 + 100000, cache=True)
        self.assertTrue(abs((v1 + v1n).coalesce() - d) == 0)

        v2n = query_flag('L1', 'CBC_CAT2', 1126051217,  1126051217 + 100000, cache=True)
        v2 = query_flag('L1', 'CBC_CAT2_VETO', 1126051217,  1126051217 + 100000, cache=True)
        self.assertTrue(abs((v2 + v2n).coalesce() - d) == 0)

    def test_cumulative_query(self):
        segs1 = query_flag('H1', 'CBC_HW_INJ',
                          1126051217,  1126051217 + 100000, cache=True)
        segs2 = query_flag('H1', 'BURST_HW_INJ',
                          1126051217,  1126051217 + 100000, cache=True)
        segs = query_cumulative_flags('H1', ['CBC_HW_INJ', 'BURST_HW_INJ'],
                                      1126051217,  1126051217 + 100000, cache=True)
        csegs = (segs1 + segs2).coalesce()
        self.assertTrue(abs(csegs - segs) == 0)

    def test_bounds(self):
        segs = query_cumulative_flags('H1', ['DATA'],
                                      1126051217,  1126051217 + 100000, cache=True)

        segs_all = query_cumulative_flags('H1', ['DATA'],
                                      1126051217,  1126051217 + 100000,
                                      bounds={'DATA':(1126051217, 1126051217 + 100000)}, cache=True)
        segs_none = query_cumulative_flags('H1', ['DATA'],
                                      1126051217,  1126051217 + 100000,
                                      bounds={'DATA':(0, 10)}, cache=True)
        self.assertTrue(abs(segs) == abs(segs_all))
        self.assertTrue(abs(segs) > abs(segs_none))

    def test_padding(self):
        segs = query_cumulative_flags('H1', ['DATA'],
                                      1126051217,  1126051217 + 100000, cache=True)

        segs2 = query_cumulative_flags('H1', ['DATA'],
                                      1126051217,  1126051217 + 100000,
                                      padding={'DATA':(8, -8)}, cache=True)
        self.assertTrue(abs(segs) > abs(segs2))

    def test_query_str(self):
        d = query_flag('H1', 'data', 1126051217, 1126051217 + 100000)

        d1 = query_str('H1', '+data', 1126051217, 1126051217 + 100000)
        d2 = query_str('H1', '+H1:data', 1126051217, 1126051217 + 100000)
        d3 = query_str('H1', '+data[1126051217:1127051217]',
                       1126051217, 1126051217 + 100000)
        d4 = query_str('H1', '+data<0:0>[1126051217:1127051217]',
                       1126051217, 1126051217 + 100000)

        self.assertTrue(abs(d - d1) == 0)
        self.assertTrue(abs(d - d2) == 0)
        self.assertTrue(abs(d - d3) == 0)
        self.assertTrue(abs(d - d4) == 0)

        d5 = query_str('L1', '+data', 1126051217, 1126051217 + 100000)
        d6 = query_str('H1', '+L1:data', 1126051217, 1126051217 + 100000)
        self.assertTrue(abs(d5 - d6) == 0)


        d7 = query_flag('H1', 'CBC_HW_INJ',
                          1126051217,  1126051217 + 100000, cache=True)
        d8 = query_str('H1', '+data,-CBC_HW_INJ', 1126051217, 1126051217 + 100000)
        self.assertTrue(abs((d-d7) - d8) == 0)

suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestDataQualityFlags))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
