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
These are the unittests for the pycbc.waveform module
"""
import pycbc
import unittest
from pycbc.types import zeros, complex64
from pycbc.filter import overlap
from pycbc.waveform import get_fd_waveform, get_waveform_filter
from utils import parse_args_all_schemes, simple_exit

_scheme, _context = parse_args_all_schemes("Waveform")

class TestSPAtmplt(unittest.TestCase):
    def setUp(self,*args):
        self.context = _context
        self.scheme = _scheme

    def test_spatmplt(self):
        fl = 25
        delta_f = 1.0 / 256

        for m1 in [1, 1.4, 20]:
            for m2 in [1.4, 20]:
                for s1 in  [-2, -1, -0.5, 0, 0.5, 1, 2]:
                    for s2 in [-2, -1, -0.5, 0, 0.5, 1, 2]:
                        # Generate TaylorF2 from lalsimulation, restricting to the capabilities of spatmplt
                        hpr,_ = get_fd_waveform( mass1=m1, mass2=m2, spin1z=s1, spin2z=s2, 
                                                 delta_f=delta_f, f_lower=fl,
                                                 approximant="TaylorF2", amplitude_order=0,
                                                 spin_order=-1, phase_order=-1)
                        hpr=hpr.astype(complex64)

                        with self.context:
                            # Generate the spatmplt waveform
                            out = zeros(len(hpr), dtype=complex64)
                            hp = get_waveform_filter(out, mass1=m1, mass2=m2, spin1z=s1, spin2z=s2,
                                                     delta_f=delta_f, f_lower=fl, approximant="SPAtmplt", 
                                                     amplitude_order=0, spin_order=-1, phase_order=-1)

                            # Check the diff is sane
                            mag = abs(hpr).sum()
                            diff = abs(hp - hpr).sum() / mag
                            self.assertTrue(diff < 0.01)

                            # Point to point overlap (no phase or time maximization)
                            o =  overlap(hp, hpr)
                            self.assertAlmostEqual(1.0, o, places=4)

                            print("checked m1: %s m2:: %s s1z: %s s2z: %s] overlap = %s, diff = %s" % (m1, m2, s1, s2, o, diff))


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestSPAtmplt))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
